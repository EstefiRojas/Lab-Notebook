#!/bin/bash

# ==============================================================================
# Script: filter_na_prob.sh
# Description: This script filters the output from 'filter_unique_hits.sh'.
#              It identifies records that have an assigned lncRNA_ENSG_ID
#              (i.e., not 'NA') but a 'highest_prob' value of 'NA'. From this
#              filtered set, it keeps only the first record for each unique
#              lncRNA_ENSG_ID.
# Date: September 23, 2025
#
# Prerequisites:
#    - awk
#
# Usage:
#    1. Run this script from the project's root directory.
#    2. Make the script executable: chmod +x scripts/filter_na_prob.sh
#    3. Run the script with the input file as an argument:
#       ./scripts/filter_na_prob.sh results/gRNA_lncRNA_matches_unique_sorted.tsv
#
# ==============================================================================

# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error when substituting.
set -u

# --- STEP 0: Check for input file and set up directories ---
if [ -z "${1-}" ]; then
    echo "ERROR: An input file is required."
    echo "Usage: ./scripts/filter_na_prob.sh <path_to_input_report.tsv>"
    exit 1
fi
INPUT_FILE="$1"

if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file '$INPUT_FILE' not found."
    exit 1
fi

# Define the output file path
RESULTS_DIR="../results"
mkdir -p "$RESULTS_DIR"
OUTPUT_FILE="$RESULTS_DIR/gRNA_lncRNA_matches_unique_ensg_na_prob.tsv"

echo "--- Filtering for records with a unique ENSG ID but NA probability ---"
echo "Input file: $INPUT_FILE"

# Use awk to filter the report.
# We keep the header and any line where the lncRNA_ENSG_ID (col 9) is not 'NA',
# the highest_prob (col 13) is 'NA', AND we have not seen the ENSG ID before.
awk '
BEGIN {
    FS="\t";
    OFS="\t";
}
# Print the header line
FNR == 1 {
    print $0;
    next;
}
# Check conditions
{
    # Condition 1: ENSG ID is not "NA"
    # Condition 2: highest_prob is "NA"
    # Condition 3: We have not seen this ENSG ID before
    if ($5 != "NA" && $7 == "NA" && !($5 in seen_ensg)) {
        print $0;
        seen_ensg[$5] = 1; # Mark this ENSG ID as seen
    }
}
' "$INPUT_FILE" > "$OUTPUT_FILE"

echo ""
echo "--- Filtering Complete! ---"
echo "Filtered report is saved in: '$OUTPUT_FILE'"

