#!/bin/bash

# ==============================================================================
# Script: add_probability.sh
# Description: This script joins the main gRNA analysis report with a second
#              file containing functional probability data. It matches records
#              by the lncRNA Ensembl Gene ID, appends the "highest_prob"
#              value, and reorders columns to prioritize the Target_Gene_ID.
# Author: Gemini
# Date: September 23, 2025
#
# Prerequisites:
#   - awk
#
# Usage:
#   1. Run this script from the project's root directory.
#   2. Make the script executable: chmod +x scripts/add_probability.sh
#   3. Run the script with the two required files as arguments:
#      ./scripts/add_probability.sh results/gRNA_lncRNA_matches.tsv path/to/probability_file.csv
#
# ==============================================================================

# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error when substituting.
set -u

# --- STEP 0: Check for input files and set up paths ---
if [ -z "${2-}" ]; then
    echo "ERROR: Two input files are required."
    echo "Usage: ./scripts/add_probability.sh <main_report.tsv> <probability_data.csv>"
    exit 1
fi
MAIN_REPORT_FILE="$1"
PROBABILITY_FILE="$2"

if [ ! -f "$MAIN_REPORT_FILE" ]; then
    echo "ERROR: Main report file '$MAIN_REPORT_FILE' not found."
    exit 1
fi
if [ ! -f "$PROBABILITY_FILE" ]; then
    echo "ERROR: Probability data file '$PROBABILITY_FILE' not found."
    exit 1
fi

# Define the output file path
RESULTS_DIR="../results"
OUTPUT_FILE="$RESULTS_DIR/gRNA_lncRNA_matches_with_prob.tsv"

echo "--- Joining reports ---"
echo "Input 1 (Main Report): $MAIN_REPORT_FILE"
echo "Input 2 (Probability Data): $PROBABILITY_FILE"

# Use awk to perform the join operation.
# It first reads the probability file into memory, creating a map.
# Then, it reads the main report, swaps the first and third columns, and
# appends the probability for each matching gene ID.
awk '
BEGIN {
    OFS="\t";
}
# Process the first file (probability data, which is a CSV)
FILENAME == ARGV[1] {
    FS=",";
    # Skip header
    if (FNR == 1) next;
    
    gene_id = $1;
    highest_prob = $6;
    
    # Strip version number (e.g., .14) from gene ID for a robust match
    sub(/\..*/, "", gene_id);
    
    # Store the probability in a map with the gene ID as the key
    prob_map[gene_id] = highest_prob;
    
    next; # Move to the next line of the first file
}
# Process the second file (main report, which is a TSV)
FILENAME == ARGV[2] {
    FS="\t";
    # Handle the header line
    if (FNR == 1) {
        # Swap columns 1 and 3 in the header
        temp = $1; $1 = $3; $3 = temp;
        print $0, "highest_prob";
        next;
    }
    
    # Swap columns 1 (gRNA_ID) and 3 (Target_Gene_ID) for the data line
    temp = $1; $1 = $3; $3 = temp;
    
    # The lncRNA ENSG ID is in the 9th column of the original file
    lncRNA_ensg_id = $9;
    
    # If there is no lncRNA hit, just print the modified line with NA for the new column
    if (lncRNA_ensg_id == "NA") {
        print $0, "NA";
        next;
    }
    
    # Strip version number from the ENSG ID to match the key format
    key = lncRNA_ensg_id;
    sub(/\..*/, "", key);
    
    # Check if a probability exists for this gene ID
    if (key in prob_map) {
        # If found, print the modified line plus the probability
        print $0, prob_map[key];
    } else {
        # If not found, print the modified line plus NA
        print $0, "NA";
    }
}
' "$PROBABILITY_FILE" "$MAIN_REPORT_FILE" > "$OUTPUT_FILE"

echo ""
echo "--- Join Complete! ---"
echo "Final enriched report is saved in: '$OUTPUT_FILE'"

