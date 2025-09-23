#!/bin/bash

# ==============================================================================
# Script: filter_unique_hits.sh
# Description: This script processes the output from the main find_lncRNA_guides.sh
#              script to create a filtered and sorted report. For each unique
#              combination of a Target_Gene_ID and a lncRNA_ENSG_ID, it keeps
#              only a single entry, prioritizing "gene" hits. The final output
#              is then sorted by the Target_Gene_ID.
# Date: September 23, 2025
#
# Prerequisites:
#   - awk, head, tail, sort
#
# Usage:
#   1. Run this script from the project's root directory.
#   2. Make the script executable: chmod +x scripts/filter_unique_hits.sh
#   3. Run the script with the input file as an argument:
#      ./scripts/filter_unique_hits.sh results/gRNA_lncRNA_matches.tsv
#
# ==============================================================================

# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error when substituting.
set -u

# --- STEP 0: Check for input file and set up directories ---
if [ -z "${1-}" ]; then
    echo "ERROR: An input file is required."
    echo "Usage: ./scripts/filter_unique_hits.sh <path_to_input_report.tsv>"
    exit 1
fi
INPUT_FILE="$1"

if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file '$INPUT_FILE' not found."
    exit 1
fi

# Define directory paths and a temporary file
RESULTS_DIR="../results"
TMP_DIR="tmp"
mkdir -p "$RESULTS_DIR" "$TMP_DIR"

TEMP_OUTPUT_FILE="$TMP_DIR/filtered_unsorted.tsv"
OUTPUT_FILE="$RESULTS_DIR/gRNA_lncRNA_matches_unique_sorted.tsv"

echo "--- Filtering Report: $INPUT_FILE ---"

# Use awk to process the report and save the unique hits to a temporary file.
# For each unique key (Target_Gene_ID + lncRNA_ENSG_ID), we store the "best" line found so far.
# A "gene" line will always replace a "transcript" line for the same key.
awk '
BEGIN {
    FS="\t";
    OFS="\t";
}
# Skip the header line, but save it for later
FNR == 1 {
    header = $0;
    next;
}
# Process all other lines
{
    # The key is the combination of Target_Gene_ID (col 3) and lncRNA_ENSG_ID (col 9)
    key = $3 FS $9;
    current_type = $8;
    
    # If we have not seen this key before, store the current line and its type.
    if (!(key in best_line)) {
        best_line[key] = $0;
        line_type[key] = current_type;
    } else {
        # If we have seen this key, we check if the new line is better.
        stored_type = line_type[key];
        # A "gene" line is always better than a "transcript" line.
        if (stored_type == "transcript" && current_type == "gene") {
            best_line[key] = $0;
            line_type[key] = current_type;
        }
    }
}
END {
    # Print the stored header
    print header;
    
    # Print the final, filtered results
    for (k in best_line) {
        print best_line[k];
    }
}
' "$INPUT_FILE" > "$TEMP_OUTPUT_FILE"

echo "--- Sorting filtered report ---"

# Sort the temporary file by the 3rd column (Target_Gene_ID), keeping the header on top.
(head -n 1 "$TEMP_OUTPUT_FILE" && tail -n +2 "$TEMP_OUTPUT_FILE" | sort -t $'\t' -k3,3) > "$OUTPUT_FILE"

# Clean up the temporary file
rm "$TEMP_OUTPUT_FILE"

echo ""
echo "--- Filtering and Sorting Complete! ---"
echo "Unique, sorted report is saved in: '$OUTPUT_FILE'"

