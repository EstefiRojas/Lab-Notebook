#!/bin/bash

# Script to add a match_pc column to a CSV file based on Target_Gene_ID presence in a PSL file
# Usage: ./add_match_column.sh input.csv input.psl output.csv

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

# Create a temporary file to store Target_Gene_IDs from PSL
TEMP_IDS=$(mktemp)

echo "Processing PSL file to extract Target_Gene_IDs..."

# Extract Target_Gene_IDs from PSL file
# The query name format appears to be: gRNA_ID_Target_Gene_ID_lncRNA
# We need to extract the Target_Gene_ID part (e.g., Hum_XLOC_000060 from gL_000047_Hum_XLOC_000060_lncRNA)

# Skip header lines in PSL (first 5 lines typically) and extract Target_Gene_IDs
awk 'NR > 5 {
    # Column 10 is the query name
    query = $10
    # Split by underscore
    n = split(query, parts, "_")
    if (n >= 4) {
        # Reconstruct Target_Gene_ID (typically Hum_XLOC_XXXXXX)
        # This assumes format: gRNA_ID_Hum_XLOC_NUMBER_lncRNA
        target_id = ""
        # Find where "Hum" starts and concatenate until before "lncRNA"
        start_found = 0
        for (i = 1; i <= n; i++) {
            if (parts[i] == "Hum") {
                start_found = 1
            }
            if (start_found && parts[i] != "lncRNA") {
                if (target_id == "") {
                    target_id = parts[i]
                } else {
                    target_id = target_id "_" parts[i]
                }
            }
            if (parts[i] == "lncRNA") {
                break
            }
        }
        if (target_id != "") {
            print target_id
        }
    }
}' "$INPUT_PSL" | sort -u > "$TEMP_IDS"

# Count how many unique Target_Gene_IDs found
NUM_IDS=$(wc -l < "$TEMP_IDS")
echo "Found $NUM_IDS unique Target_Gene_IDs in PSL file"

echo "Processing CSV file and adding match_pc column..."

# Process the CSV file
# Read header and add new column
head -n 1 "$INPUT_CSV" | sed 's/$/\tmatch_pc/' > "$OUTPUT_CSV"

# Process data rows
tail -n +2 "$INPUT_CSV" | while IFS=$'\t' read -r line; do
    # Extract Target_Gene_ID (first column)
    target_gene_id=$(echo "$line" | cut -f1)
    
    # Check if this ID exists in our PSL-derived list
    if grep -q "^${target_gene_id}$" "$TEMP_IDS"; then
        match_value="YES"
    else
        match_value="NO"
    fi
    
    # Append the match value to the line
    echo -e "${line}\t${match_value}"
done >> "$OUTPUT_CSV"

# Clean up temporary file
rm -f "$TEMP_IDS"

# Count results
TOTAL_ROWS=$(tail -n +2 "$OUTPUT_CSV" | wc -l)
MATCHED_ROWS=$(tail -n +2 "$OUTPUT_CSV" | grep -c $'\tYES$')
UNMATCHED_ROWS=$(tail -n +2 "$OUTPUT_CSV" | grep -c $'\tNO$')

echo "Processing complete!"
echo "Output written to: $OUTPUT_CSV"
echo "Total data rows: $TOTAL_ROWS"
echo "Rows with matches: $MATCHED_ROWS"
echo "Rows without matches: $UNMATCHED_ROWS"
