#!/opt/local/bin/bash

# CORRECTED SCRIPT (v3)
# This version matches the gRNA_ID from the CSV (col 3) with the
# gRNA_ID extracted from the PSL (col 10).
# It correctly reads the strand from CSV col 6.
# It also now REPLACES columns 8 and 9 instead of appending.
#
# Usage: ./correct_match_script.sh input.csv input.psl output.csv

# Check bash version (associative arrays require bash 4+)
if [ "${BASH_VERSION%%.*}" -lt 4 ]; then
    echo "Error: This script requires bash version 4 or higher"
    echo "Current version: $BASH_VERSION"
    exit 1
fi

# Check if correct number of arguments provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_csv> <input_psl> <output_csv>"
    echo "Example: $0 file1.csv blat_output.psl file1_with_match.csv"
    exit 1
fi

INPUT_CSV="$1"
INPUT_PSL="$2"
OUTPUT_CSV="$3"

# Check if input files exist
if [ ! -f "$INPUT_CSV" ]; then
    echo "Error: Input CSV file '$INPUT_CSV' not found!"
    exit 1
fi

if [ ! -f "$INPUT_PSL" ]; then
    echo "Error: Input PSL file '$INPUT_PSL' not found!"
    exit 1
fi

# Declare two associative arrays
# psl_strands stores the extracted strand (+ or -)
# psl_col14 stores the full value of the 14th column
# *** KEY IS gRNA_ID ***
declare -A psl_strands
declare -A psl_col14

echo "Processing PSL file to extract gRNA_IDs, strands, and column 14..."

# Process PSL file
{
    # Skip headers (lines starting with 'psLayout', 'match', dashes, or empty lines)
    while IFS=$'\t' read -r line; do
        # Skip header lines
        if [[ "$line" =~ ^psLayout ]] || [[ "$line" =~ ^match ]] || [[ "$line" =~ ^-+ ]] || [[ -z "$line" ]]; then
            continue
        fi
        
        # Parse the line
        IFS=$'\t' read -ra fields <<< "$line"
        
        # Check for at least 14 columns
        if [ ${#fields[@]} -ge 14 ]; then
            query_name="${fields[9]}"  # Column 10 (e.g., gL_000182_Hum_XLOC_000265_lncRNA)
            target_col14="${fields[13]}" # Column 14 (e.g., sp|...(+))
            
            # Extract strand from column 14
            psl_strand=""
            if [[ "$target_col14" == *"(+)" ]]; then
                psl_strand="+"
            elif [[ "$target_col14" == *"(-)" ]]; then
                psl_strand="-"
            fi
            
            # Extract gRNA_ID from query_name
            if [[ "$query_name" =~ ^(gL_[0-9]+) ]]; then
                gRNA_id_key="${BASH_REMATCH[1]}"
                
                # Store strand and col14 value using gRNA_ID as the key
                # Note: If a gRNA_ID appears multiple times, this stores the *last* one found.
                psl_strands["$gRNA_id_key"]="$psl_strand"
                psl_col14["$gRNA_id_key"]="$target_col14"
            fi
        fi
    done
} < "$INPUT_PSL"

NUM_IDS=${#psl_strands[@]}
echo "Found $NUM_IDS unique gRNA_IDs in PSL file"

# Debug: Show first few IDs found (optional)
if [ $NUM_IDS -gt 0 ]; then
    echo "Sample gRNA_IDs found:"
    count=0
    for id in "${!psl_strands[@]}"; do
        echo "  - $id (Strand: ${psl_strands[$id]})"
        ((count++))
        [ $count -ge 3 ] && break
    done
fi

echo "Processing CSV file and adding/replacing match_pc and psl_target_col14 columns..."

# Process CSV file
{
    # Read and modify header
    IFS=$'\t' read -r header
    
    # *** CORRECTION: Read header into array to replace cols 8/9 ***
    IFS=$'\t' read -ra header_fields <<< "$header"
    
    # Assume input has at least 7 columns.
    # We will print the first 7 header columns, then our new 8th and 9th.
    if [ ${#header_fields[@]} -ge 7 ]; then
        col1_h="${header_fields[0]}"
        col2_h="${header_fields[1]}"
        col3_h="${header_fields[2]}"
        col4_h="${header_fields[3]}"
        col5_h="${header_fields[4]}"
        col6_h="${header_fields[5]}"
        col7_h="${header_fields[6]}"
        # Print first 7 headers + new 8th and 9th headers
        echo -e "${col1_h}\t${col2_h}\t${col3_h}\t${col4_h}\t${col5_h}\t${col6_h}\t${col7_h}\tmatch_pc\tpsl_target_col14"
    else
        # Fallback for 6-column file (original logic)
        echo -e "${header}\tmatch_pc\tpsl_target_col14"
    fi

    # Process data rows
    row_count=0
    matched_count=0
    unmatched_count=0
    
    while IFS=$'\t' read -r line; do
        IFS=$'\t' read -ra csv_fields <<< "$line"
        
        # *** CORRECTION: Check for at least 6 columns (for gRNA_ID and Strand) ***
        # We need col 3 and col 6, so at least 6 columns.
        if [ ${#csv_fields[@]} -lt 6 ]; then
            echo "Skipping malformed CSV line: $line" >&2
            continue
        fi
        
        # *** CORRECTION 1: Get gRNA_ID (col 3) ***
        gRNA_id="${csv_fields[2]}"
        
        # Get lncRNA_Strand (col 6)
        csv_strand_raw="${csv_fields[5]}" # e.g., (+) or (-)
        
        # Normalize CSV strand
        csv_strand_norm=""
        if [[ "$csv_strand_raw" == "(+)" ]]; then
            csv_strand_norm="+"
        elif [[ "$csv_strand_raw" == "(-)" ]]; then
            csv_strand_norm="-"
        fi
        
        # New logic with strand check
        match_value="NO"
        col14_value="NA" # Default value if no match
        
        # Check if gRNA_id exists in the PSL data (using the correct key)
        if [[ -n "${psl_strands[$gRNA_id]}" ]]; then
            # ID exists, get stored PSL strand and col14 value
            psl_strand_val="${psl_strands[$gRNA_id]}"
            col14_value="${psl_col14[$gRNA_id]}"
            
            # Check for strand coincidence (and that CSV strand was valid)
            if [[ -n "$csv_strand_norm" && "$csv_strand_norm" == "$psl_strand_val" ]]; then
                match_value="YES"
                ((matched_count++))
            else
                # ID matched, but strand did not (or CSV/PSL strand was blank/invalid)
                ((unmatched_count++))
            fi
        else
            # ID did not match
            ((unmatched_count++))
        fi
        
        # *** CORRECTION 2: Rebuild output line to replace cols 8/9 ***
        # This logic assumes the input file has at least 7 columns
        # and we want to output 9 columns total.
        if [ ${#csv_fields[@]} -ge 7 ]; then
            col1="${csv_fields[0]}"
            col2="${csv_fields[1]}"
            col3="${csv_fields[2]}" # This is gRNA_id
            col4="${csv_fields[3]}"
            col5="${csv_fields[4]}"
            col6="${csv_fields[5]}" # This is csv_strand_raw
            col7="${csv_fields[6]}" # e.g., the 0.912 column
            
            # Print first 7 cols + new 8th and 9th cols
            echo -e "${col1}\t${col2}\t${col3}\t${col4}\t${col5}\t${col6}\t${col7}\t${match_value}\t${col14_value}"
        else
            # Fallback for 6-column file (original logic)
            echo -e "${line}\t${match_value}\t${col14_value}"
        fi

        ((row_count++))
        
        # Progress indicator for large files
        if [ $((row_count % 1000)) -eq 0 ]; then
            echo "  Processed $row_count rows..." >&2
        fi
    done
    
    echo "Processing complete!" >&2
    echo "Total data rows processed: $row_count" >&2
    echo "Rows with strand-coincident matches: $matched_count" >&2
    echo "Rows without matches (or strand mismatch): $unmatched_count" >&2
    
} < "$INPUT_CSV" > "$OUTPUT_CSV".tmp

# Finally, change the names of the groups
awk 'BEGIN{
    FS="\t"; OFS="\t";
    # Define the replacement map
    map["Shared"] = "Core";
    map["Partially shared"] = "Common";
    map["Cell-type specific"] = "Specific";
}
{
    # Check if the value of column 2 is a key in our map
    if ($2 in map) {
        $2 = map[$2]; # If yes, replace it with the maps value
    }
    # Print the entire line (either modified or as-is)
    print
}' "$OUTPUT_CSV".tmp > "$OUTPUT_CSV"

rm -rf "$OUTPUT_CSV".tmp
echo "Output written to: $OUTPUT_CSV"

# Validation check
# Note: This check might fail if the input was 6 columns and output is 8
# A better check is that output lines >= input lines (due to header)
output_lines=$(wc -l < "$OUTPUT_CSV")
input_lines=$(wc -l < "$INPUT_CSV")

if [ "$output_lines" -eq "$input_lines" ]; then
    echo "✓ Validation passed: Output has same number of lines as input"
else
    echo "⚠ Warning: Line count mismatch (input: $input_lines, output: $output_lines)"
fi