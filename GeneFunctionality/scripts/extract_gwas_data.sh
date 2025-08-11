#!/bin/bash
# ===============================================
# GWAS-lncRNA Multi-Region Association Analysis with BEDTools
# ===============================================
#
# Description:
#   This script performs genomic overlap analysis between GWAS variants
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
#   gwas_catalog.tsv  - The full GWAS catalog data in TSV format
#                       Expected columns: chr (12), position (13), p-value (28), beta (31), beta_desc (32)
#   output_file.csv   - Output CSV with original lncRNA data plus statistics for each region:
#                       SNP_count, min_p_value, max_abs_beta, and sum_beta

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
Usage: $0 <lncrna_file.csv> <gwas_catalog.tsv> <output_file.csv>"
fi

LNCRNA_FILE=$1
GWAS_FILE=$2
OUTPUT_FILE=$3

# Check input files
[ ! -f "$LNCRNA_FILE" ] && error_exit "lncRNA file '$LNCRNA_FILE' not found."
[ ! -f "$GWAS_FILE" ] && error_exit "GWAS catalog file '$GWAS_FILE' not found."

# --- Configuration ---
THRESHOLD=5e-8

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
# Filter GWAS by p-value and convert to BED
awk -F'\t' -v threshold="$THRESHOLD" '
BEGIN {
    # Convert scientific notation threshold to number
    threshold_num = threshold + 0
}
NR > 1 {
    # Get fields
    chr = $12
    pos = $13
    pval = $28
    pmlog = $29 # Use the maximum of this field to get minimum p-value
    beta = $31  # Beta coefficient in column 31
    beta_desc = $32  # Beta description in column 32
    
    # Skip if essential fields are empty
    if (chr == "" || pos == "" || pval == "") next
    
    # Convert p-value to number for comparison
    pval_num = pval + 0
    pmlog_num = pmlog + 0
    
    # Check p-value threshold
    if (pval_num >= threshold_num) next
    
    # Normalize chromosome
    gsub(/^chr/, "", chr)
    
    # Validate position is numeric
    if (pos !~ /^[0-9]+$/) next
    
    pos_int = int(pos)
    if (pos_int <= 0) next
    
    # Check if beta is valid based on description in column 32
    valid_beta = 0
    if (beta_desc ~ /unit decrease/ || beta_desc ~ /unit increase/) {
        valid_beta = 1
    }
    
    # Handle beta value - use "NA" if not valid or not numeric
    if (!valid_beta || beta == "" || beta !~ /^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$/) {
        beta = "NA"
    }
    
    # Create BED entry (SNP position as single base)
    # Format: chr start end name pvalue pmlog beta
    printf "chr%s\t%d\t%d\tSNP_%d\t%s\t%s\t%s\n", chr, pos_int-1, pos_int, NR, pval, pmlog, beta
}' "$GWAS_FILE" | sort -k1,1 -k2,2n > "$FILTERED_GWAS"

# Check if file is empty or has content
if [ -s "$FILTERED_GWAS" ]; then
    GWAS_COUNT=$(wc -l < "$FILTERED_GWAS")
    echo "  -> Found $GWAS_COUNT significant GWAS variants (p < $THRESHOLD)"
    
    # Debug: Show first few lines
    echo "  -> Sample GWAS BED entries:"
    head -n 3 "$FILTERED_GWAS" | sed 's/^/     /'
else
    GWAS_COUNT=0
    echo "  -> Found $GWAS_COUNT significant GWAS variants (p < $THRESHOLD)"
fi

if [ "$GWAS_COUNT" -eq 0 ]; then
    echo "Warning: No GWAS variants passed the p-value threshold."
    echo "Creating output with 0 counts and NA values..."
    
    # Add columns for all regions with 0 and NA
    awk -F',' 'NR==1 {
        # Remove any trailing whitespace and add new columns for each region
        sub(/[[:space:]]*$/, "")
        printf "%s", $0
        printf ",tl_exon1_SNP_count,tl_exon1_min_p_value,tl_exon1_max_abs_beta,tl_exon1_sum_beta"
        printf ",tl_exon2_SNP_count,tl_exon2_min_p_value,tl_exon2_max_abs_beta,tl_exon2_sum_beta"
        printf ",tl_SNP_count,tl_min_p_value,tl_max_abs_beta,tl_sum_beta"
        printf ",tr_SNP_count,tr_min_p_value,tr_max_abs_beta,tr_sum_beta\n"
    } 
    NR>1 {
        sub(/[[:space:]]*$/, "")
        printf "%s", $0
        # Add zeros and NAs for all four regions (4 columns each)
        for (i=0; i<4; i++) {
            printf ",0,NA,NA,NA"
        }
        printf "\n"
    }' "$LNCRNA_FILE" > "$OUTPUT_FILE"
    
    echo "Done."
    exit 0
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
                    pval = $9
                    pmlog = $10
                    beta = $11
                    
                    # Track counts
                    count[line_num]++
                    
                    # Track minimum p-value
                    if (!(line_num in min_pval) || pmlog > max_pmlog[line_num]) {
                        min_pval[line_num] = pval
                        max_pmlog[line_num] = pmlog
                    }
                    
                    # Track beta statistics if valid
                    if (beta != "NA") {
                        beta_val = beta + 0  # Convert to number
                        abs_beta = (beta_val < 0) ? -beta_val : beta_val
                        
                        # Track maximum absolute beta value
                        if (!(line_num in max_abs_beta) || abs_beta > max_abs_beta[line_num]) {
                            max_abs_beta[line_num] = abs_beta
                        }
                        
                        # Track sum of betas
                        if (!(line_num in sum_beta)) {
                            sum_beta[line_num] = 0
                        }
                        sum_beta[line_num] += abs_beta
                    }
                }
                END {
                    for (ln in count) {
                        max_beta = (ln in max_abs_beta) ? max_abs_beta[ln] : "NA"
                        beta_sum = (ln in sum_beta) ? sum_beta[ln] : "NA"
                        printf "%s\t%s\t%d\t%s\t%s\t%s\n", ln, region, count[ln], min_pval[ln], max_beta, beta_sum
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
        snp_count[line_num, region] = arr[3]
        min_p[line_num, region] = arr[4]
        max_beta[line_num, region] = arr[5]
        sum_beta[line_num, region] = arr[6]
    }
    close(stats_file)
}
NR==1 {
    # Print header with new columns for each region
    sub(/[[:space:]]*$/, "")
    printf "%s", $0
    printf ",tl_exon1_SNP_count_ct,tl_exon1_min_p_value_ct,tl_exon1_max_abs_beta_ct,tl_exon1_sum_beta_ct"
    printf ",tl_exon2_SNP_count_ct,tl_exon2_min_p_value_ct,tl_exon2_max_abs_beta_ct,tl_exon2_sum_beta_ct"
    printf ",tl_SNP_count_ct,tl_min_p_value_ct,tl_max_abs_beta_ct,tl_sum_beta_ct"
    printf ",tr_SNP_count_ct,tr_min_p_value_ct,tr_max_abs_beta_ct,tr_sum_beta_ct\n"
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
        count = ((current_line, region) in snp_count) ? snp_count[current_line, region] : 0
        pval = ((current_line, region) in min_p) ? min_p[current_line, region] : "NA"
        beta_max = ((current_line, region) in max_beta) ? max_beta[current_line, region] : "NA"
        beta_sum = ((current_line, region) in sum_beta) ? sum_beta[current_line, region] : "NA"
        
        printf ",%s,%s,%s,%s", count, pval, beta_max, beta_sum
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