#!/bin/bash
# ===============================================
# GWAS ced set-lncRNA Multi-Region Association Analysis with BEDTools
# ===============================================
#
# Description:
#   This script performs genomic overlap analysis between GWAS ced set variants
#   and multiple lncRNA regions (exons and transcript boundaries) using BEDTools.
#
# Requirements:
#   - BEDTools (install: conda install -c bioconda bedtools or apt-get install bedtools)
#   - Standard Unix tools (awk, sort, join)
#
# Usage:
#   ./gwas_lncrna_analysis.sh <lncrna_file.csv> <gwas_catalog.tsv> <output_file.csv>
#
# Arguments:
#   lncrna_file.csv   - Input file with lncRNA data in CSV format
#                       Expected columns: 
#                         tl_exon1 (Chrm_Exon1): chr (10), start (11), end (12)
#                         tl_exon2 (Chrm_Exon2): chr (13), start (14), end (15)
#                         tl (Chrm_Transcript_Left): chr (16), start (17), end (18)
#                         tr (Chrm_Transcript_Right): chr (19), start (20), end (21)
#   gwas_cred_set.tsv  - The GWAS credible set data in TSV format from open_targets
#                       Expected columns: chr (7), position (8), pip (15)
#   output_file.csv   - Output CSV with original lncRNA data plus pip sum for each region:
#                       sum_pip

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# --- Functions ---
error_exit() {
    echo "Error: $1" >&2
    exit 1
}

check_command() {
    if ! command -v "$1" &> /dev/null; then
        error_exit "$1 is required but not installed. Please install it first."
    fi
}

# --- Initial Setup & Validation ---
DATE=$(date "+%A %B %d, %Y")
echo "====================================="
echo "GWAS-lncRNA Multi-Region Analysis"
echo "Execution date: ${DATE}"
echo "====================================="

# Check required tools
check_command bedtools
check_command awk

# Check arguments
if [ "$#" -ne 3 ]; then
    error_exit "Invalid number of arguments.
Usage: $0 <lncrna_file.csv> <gwas_cred_set.csv> <output_file.csv>"
fi

LNCRNA_FILE=$1
GWAS_FILE=$2
OUTPUT_FILE=$3

# Check input files
[ ! -f "$LNCRNA_FILE" ] && error_exit "lncRNA file '$LNCRNA_FILE' not found."
[ ! -f "$GWAS_FILE" ] && error_exit "GWAS credible sets file '$GWAS_FILE' not found."

# Create temporary files
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

# Define region types
REGIONS=("tl_exon1" "tl_exon2" "tl" "tr")
REGION_COLS=(10 13 16 19)  # Starting column for each region (chr column)

# Temporary files
FILTERED_GWAS="$TEMP_DIR/gwas_filtered.bed"
ALL_STATS="$TEMP_DIR/all_stats.txt"

echo "Step 1: Creating BED files for all lncRNA regions..."

# First, let's check the number of columns in the file for debugging
NCOLS=$(head -n 1 "$LNCRNA_FILE" | awk -F',' '{print NF}')
echo "  -> Input file has $NCOLS columns"

# Create separate BED files for each region type
for i in "${!REGIONS[@]}"; do
    REGION="${REGIONS[$i]}"
    START_COL="${REGION_COLS[$i]}"
    
    BED_FILE="$TEMP_DIR/${REGION}.bed"
    
    echo "  -> Extracting ${REGION} from columns ${START_COL}, $((START_COL+1)), $((START_COL+2))..."
    
    # Extract coordinates for this region
    awk -F',' -v start_col="$START_COL" -v region="$REGION" 'NR>1 {
        # Get chromosome, start, end from appropriate columns
        chr = $(start_col)
        start = $(start_col + 1)
        end = $(start_col + 2)
        
        # Remove any whitespace, carriage returns, and newlines from all fields
        gsub(/[[:space:]]+$/, "", chr)
        gsub(/^[[:space:]]+/, "", chr)
        gsub(/[[:space:]]+$/, "", start)
        gsub(/^[[:space:]]+/, "", start)
        gsub(/[[:space:]]+$/, "", end)
        gsub(/^[[:space:]]+/, "", end)
        gsub(/\r/, "", chr)
        gsub(/\r/, "", start)
        gsub(/\r/, "", end)
        
        # Debug: print first lines values for tr region
        if (NR == 2 && region == "tr") {
            print "     Debug " region " after cleaning:" > "/dev/stderr"
            print "       chr: [" chr "]" > "/dev/stderr"
            print "       start: [" start "]" > "/dev/stderr"
            print "       end: [" end "]" > "/dev/stderr"
        }
        
        # Create unique ID: line number and region type
        id = (NR - 1) "_" region
        
        # Normalize chromosome names
        gsub(/^chr/, "", chr)
        
        # Validate and convert coordinates
        if (chr != "" && start ~ /^[0-9]+$/ && end ~ /^[0-9]+$/) {
            # Convert to integers and 0-based
            start_int = int(start) - 1
            end_int = int(end)
            
            if (start_int >= 0 && end_int > start_int) {
                printf "chr%s\t%d\t%d\t%s\n", chr, start_int, end_int, id
            }
        }
    }' "$LNCRNA_FILE" | sort -k1,1 -k2,2n > "$BED_FILE" 2>&1
    
    COUNT=$(wc -l < "$BED_FILE")
    echo "  -> Processed $COUNT ${REGION} regions"
done

# Combine all regions into one BED file for summary
cat "$TEMP_DIR"/tl_exon1.bed "$TEMP_DIR"/tl_exon2.bed "$TEMP_DIR"/tl.bed "$TEMP_DIR"/tr.bed | sort -k1,1 -k2,2n > "$TEMP_DIR/all_regions.bed"
TOTAL_REGIONS=$(wc -l < "$TEMP_DIR/all_regions.bed")
echo "  -> Total regions to analyze: $TOTAL_REGIONS"

echo "Step 2: Filtering and converting GWAS data to BED format..."

# First, determine the correct field separator based on file extension
if [[ "$GWAS_FILE" == *.tsv ]]; then
    FIELD_SEP='\t'
    echo "  -> Detected TSV format for GWAS file"
elif [[ "$GWAS_FILE" == *.csv ]]; then
    FIELD_SEP=','
    echo "  -> Detected CSV format for GWAS file"
else
    # Try to auto-detect by checking first line
    if head -n 1 "$GWAS_FILE" | grep -q $'\t'; then
        FIELD_SEP='\t'
        echo "  -> Auto-detected TSV format for GWAS file"
    else
        FIELD_SEP=','
        echo "  -> Auto-detected CSV format for GWAS file"
    fi
fi

# Filter GWAS by p-value and convert to BED format
awk -F"$FIELD_SEP" 'NR > 1 {
    # Get fields
    chr = $7
    pos = $8
    pip = $15  # pip of variant

    # Remove any whitespace and carriage returns
    gsub(/[[:space:]]+$/, "", chr)
    gsub(/^[[:space:]]+/, "", chr)
    gsub(/[[:space:]]+$/, "", pos)
    gsub(/^[[:space:]]+/, "", pos)
    gsub(/[[:space:]]+$/, "", pip)
    gsub(/^[[:space:]]+/, "", pip)
    gsub(/\r/, "", chr)
    gsub(/\r/, "", pos)
    gsub(/\r/, "", pip)
    
    # Debug: print first few lines
    if (NR <= 3) {
        print "     Debug GWAS line " NR ":" > "/dev/stderr"
        print "       chr: [" chr "]" > "/dev/stderr"
        print "       pos: [" pos "]" > "/dev/stderr"
        print "       pip: [" pip "]" > "/dev/stderr"
    }
    
    # Normalize chromosome (handle various formats)
    gsub(/^chr/, "", chr)
    gsub(/^CHR/, "", chr)
    gsub(/^Chr/, "", chr)
    
    # Convert position to integer
    pos_int = int(pos)
    
    
    # Validate PIP (can be decimal or NA)
    if (pip == "" || pip == "NA" || pip == "NULL") {
        pip = "0"  # Default to 0 if missing
    }
    
    # Create BED entry (SNP position as single base)
    # Format: chr start end name pip
    printf "chr%s\t%d\t%d\tSNP_%d\t%s\n", chr, pos_int-1, pos_int, NR, pip
}' "$GWAS_FILE" | grep -v "^     Debug" | sort -k1,1 -k2,2n > "$FILTERED_GWAS"

# Check if file is empty or has content
if [ -s "$FILTERED_GWAS" ]; then
    GWAS_COUNT=$(wc -l < "$FILTERED_GWAS")
    echo "  -> Found $GWAS_COUNT significant GWAS variants"
    
    # Debug: Show first few lines
    echo "  -> Sample GWAS BED entries:"
    head -n 3 "$FILTERED_GWAS" | sed 's/^/     /'
else
    GWAS_COUNT=0
    echo "  -> WARNING: No valid GWAS variants found!"
    echo "  -> Checking raw GWAS file format..."
    echo "     First 2 lines of GWAS file:"
    head -n 2 "$GWAS_FILE" | sed 's/^/     /'
    echo "     Total lines in GWAS file: $(wc -l < "$GWAS_FILE")"
    
    # Additional debugging - check specific columns
    echo "  -> Sampling column 7 (chr), 8 (pos), 15 (pip) from first 3 data rows:"
    awk -F"$FIELD_SEP" 'NR > 1 && NR <= 4 {
        print "     Row " (NR-1) ": chr=[" $7 "], pos=[" $8 "], pip=[" $15 "]"
    }' "$GWAS_FILE"
fi

echo "Step 3: Finding overlaps for each region type..."
# Initialize stats file
> "$ALL_STATS"

# Process each region type
for REGION in "${REGIONS[@]}"; do
    BED_FILE="$TEMP_DIR/${REGION}.bed"
    INTERSECT_FILE="$TEMP_DIR/${REGION}_intersect.bed"
    
    echo "  -> Processing ${REGION} regions..."
    
    # Find overlaps for this region type
    if [ -s "$BED_FILE" ]; then
        if bedtools intersect -wa -wb -a "$BED_FILE" -b "$FILTERED_GWAS" > "$INTERSECT_FILE" 2>/dev/null; then
            if [ -s "$INTERSECT_FILE" ]; then
                OVERLAP_COUNT=$(wc -l < "$INTERSECT_FILE")
                echo "     Found $OVERLAP_COUNT overlaps"
                
                # Calculate statistics for this region type
                awk -F'\t' -v region="$REGION" '{
                    # Extract ID and statistics
                    id = $4  # Format: linenum_regiontype
                    split(id, parts, "_")
                    line_num = parts[1]
                    pip = $9
                    
                    # Track counts
                    count[line_num]++
                    
                    # Track beta statistics if valid
                    if (pip != "NA") {
                        pip_val = pip + 0  # Convert to number    
                    
                        # Track sum of pips
                        if (!(line_num in sum_pip)) {
                            sum_pip[line_num] = 0
                        }
                        sum_pip[line_num] += pip_val
                    }
                }
                END {
                    for (ln in count) {
                        pip_sum = (ln in sum_pip) ? sum_pip[ln] : "NA"
                        printf "%s\t%s\t%d\t%s\n", ln, region, count[ln], pip_sum
                    }
                }' "$INTERSECT_FILE" >> "$ALL_STATS"
            else
                echo "     No overlaps found"
            fi
        fi
    fi
done

echo "Step 4: Merging results with original lncRNA file..."
# Create the final output by merging all statistics
awk -F',' -v stats_file="$ALL_STATS" 'BEGIN {
    OFS=","
    
    # Read all statistics
    while ((getline line < stats_file) > 0) {
        split(line, arr, "\t")
        line_num = arr[1]
        region = arr[2]
        
        # Store stats indexed by line number and region
        count[line_num, region] = arr[3]
        sum_pip[line_num, region] = arr[4]
    }
    close(stats_file)
}
NR==1 {
    # Print header with new columns for each region
    sub(/[[:space:]]*$/, "")
    printf "%s", $0
    printf ",tl_exon1_cred_set_count,tl_exon1_sum_pip"
    printf ",tl_exon2_cred_set_count,tl_exon2_sum_pip"
    printf ",tl_cred_set_count,tl_sum_pip"
    printf ",tr_cred_set_count,tr_sum_pip\n"
    next
}
{
    # Process each lncRNA line
    current_line = NR - 1
    
    # Remove trailing whitespace from original line
    sub(/[[:space:]]*$/, "")
    printf "%s", $0
    
    # Add statistics for each region
    regions[1] = "tl_exon1"
    regions[2] = "tl_exon2"
    regions[3] = "tl"
    regions[4] = "tr"
    
    for (i=1; i<=4; i++) {
        region = regions[i]
        
        # Get values or set defaults
        cred_set_count = ((current_line, region) in count) ? count[current_line, region] : "NA"
        pip_sum = ((current_line, region) in sum_pip) ? sum_pip[current_line, region] : "NA"
        
        printf ",%s,%s", cred_set_count, pip_sum
    }
    printf "\n"
}' "$LNCRNA_FILE" > "$OUTPUT_FILE"

# --- Summary Statistics ---
echo "====================================="
echo "Analysis Complete!"
echo "====================================="
echo "Summary:"
echo "  - Total coordinate sets analyzed: $TOTAL_REGIONS"
echo "  - Significant GWAS variants: $GWAS_COUNT"

# Count overlaps per region type
for REGION in "${REGIONS[@]}"; do
    INTERSECT_FILE="$TEMP_DIR/${REGION}_intersect.bed"
    if [ -s "$INTERSECT_FILE" ]; then
        OVERLAPS=$(wc -l < "$INTERSECT_FILE")
        UNIQUE=$(cut -f4 "$INTERSECT_FILE" | cut -d'_' -f1 | sort -u | wc -l)
        echo "  - ${REGION}: $OVERLAPS overlaps in $UNIQUE lncRNAs"
    else
        echo "  - ${REGION}: No overlaps found"
    fi
done

echo "  - Output written to: $OUTPUT_FILE"
echo ""
echo "First few lines of output (truncated for display):"
head -n 3 "$OUTPUT_FILE" | cut -c1-150
echo "..."
echo "Done."