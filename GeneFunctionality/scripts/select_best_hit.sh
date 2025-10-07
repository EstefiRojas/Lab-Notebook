#!/bin/bash
# ==============================================================================
# Script: select_best_hit.sh
# Description: This script filters an enriched gRNA report to find the best
#              lncRNA hit for each unique Target_Gene_ID. If a gene has
#              multiple lncRNA matches, it determines the best one by finding
#              the highest probability. It then keeps ALL gRNA records that
#              target that single best lncRNA. If a gene only has non-matches
#              (NA), it keeps one of those records. The final output is sorted.
#
# Date: September 24, 2025
#
# Prerequisites:
#   - awk, head, tail, sort
#
# Usage:
#   1. Run this script from the project's root directory.
#   2. Make the script executable: chmod +x scripts/select_best_hit.sh
#   3. Run the script with the enriched report file as an argument:
#      ./scripts/select_best_hit.sh results/gRNA_lncRNA_matches_with_prob.tsv
#
# ==============================================================================
# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error when substituting.
set -u
# --- STEP 0: Check for input file and set up paths ---
if [ -z "${1-}" ]; then
    echo "ERROR: An input file is required."
    echo "Usage: ./scripts/select_best_hit.sh <path_to_enriched_report.tsv>"
    exit 1
fi
INPUT_FILE="$1"
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file '$INPUT_FILE' not found."
    exit 1
fi
# Define directory paths and a temporary file
RESULTS_DIR="../results"
TMP_DIR="../tmp"
mkdir -p "$RESULTS_DIR" "$TMP_DIR"
TEMP_OUTPUT_FILE="$TMP_DIR/best_hits_unsorted.tsv"
OUTPUT_FILE="$RESULTS_DIR/gRNA_lncRNA_matches_best_hit_multi.tsv"
echo "--- Selecting best hit(s) for each Target_Gene_ID from: $INPUT_FILE ---"
# Use a two-pass awk approach. The input file is passed twice.
awk '
BEGIN {
    FS="\t";
    OFS="\t";
}
# PASS 1: Iterate through the file to find the max probability and corresponding
# best ENSG ID for each Target_Gene_ID.
FNR==NR {
    if (FNR == 1) next; # Skip header on first pass
    
    target_gene = $1;  # Only use Target_Gene_ID as the key
    ensg_id = $9;
    prob = $13;
    
    # Only consider non-NA entries
    if (ensg_id != "NA" && prob != "NA") {
        # Convert prob to number for comparison
        prob_num = prob + 0;
        
        # If we haven'\''t stored a probability for this target gene yet, or if the
        # current line has a higher probability, update the best hit.
        if (! (target_gene in max_prob) || prob_num > max_prob[target_gene]) {
            max_prob[target_gene] = prob_num;
            best_ensg[target_gene] = ensg_id;
        }
    }
    next; # Move to the next line of the first pass
}
# PASS 2: Iterate through the file again. Print lines that match our criteria.
{
    if (FNR == 1) { # Print header on second pass
        print $0;
        next;
    }
    
    target_gene = $1;
    current_ensg = $9;
    
    # Case A: This Target Gene had successful matches.
    # Print all lines where the ENSG ID matches the best one we found.
    if (target_gene in best_ensg) {
        if (current_ensg == best_ensg[target_gene]) {
            print $0;
        }
    }
    # Case B: This Target Gene had NO successful matches (all were NA).
    # We need to print one representative NA line.
    else {
        if (! (target_gene in printed_na)) {
            print $0;
            printed_na[target_gene] = 1; # Mark this target gene as printed
        }
    }
}
' "$INPUT_FILE" "$INPUT_FILE" > "$TEMP_OUTPUT_FILE"
echo "--- Sorting final report ---"
# Sort the temporary file by the 1st column (Target_Gene_ID), keeping the header on top.
(head -n 1 "$TEMP_OUTPUT_FILE" && tail -n +2 "$TEMP_OUTPUT_FILE" | sort -t $'\t' -k1,1) > "$OUTPUT_FILE"
# Clean up the temporary file
rm "$TEMP_OUTPUT_FILE"
echo ""
echo "--- Selection Complete! ---"
echo "Final report of best hits is saved in: '$OUTPUT_FILE'"