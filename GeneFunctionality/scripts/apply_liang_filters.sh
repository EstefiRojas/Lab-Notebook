#!/bin/bash
set -ueo pipefail

# ==============================================================================
# Script: apply_liang_filters.sh
# Description: Applies filtering logic from the Liang study to the unified dataset.
#              1. Deduplicates hits (Unique Hits).
#              2. Selects the best lncRNA target for each gene based on probability.
#
# Usage: ./apply_liang_filters.sh <input_csv> <output_csv>
# ==============================================================================

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_csv> <output_csv>"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"
TMP_DIR="../tmp"
mkdir -p "$TMP_DIR"

TEMP_UNIQUE="$TMP_DIR/temp_unique.csv"
TEMP_BEST="$TMP_DIR/temp_best.csv"

echo "Applying Liang filters..."

# ---------------------------------------------------------
# Step 1: Unique Hits
# Deduplicate based on gRNA_ID, Target_Gene_ID, and ENSG_ID.
# ---------------------------------------------------------
echo "Step 1: Filtering for unique hits..."

awk -F',' '
BEGIN { OFS="," }
NR==1 { print $0; next }
{
    # Columns: 
    # 1:Study, 2:gRNA_ID, 3:Target_Gene_ID, 4:gRNA_Type, 5:Essentiality, 
    # 6:ENSG_ID, 7:Strand, 8:Probability_Functional, 9:Protein_Off_Target
    
    # Key = gRNA_ID + Target_Gene_ID + ENSG_ID
    key = $2 SUBSEP $3 SUBSEP $6
    
    if (!(key in seen)) {
        seen[key] = 1
        print $0
    }
}' "$INPUT_FILE" > "$TEMP_UNIQUE"

echo "Unique hits saved to temporary file."

# ---------------------------------------------------------
# Step 2: Best Hit Selection
# For each Target_Gene_ID, find the ENSG_ID with the highest Probability_Functional.
# Then keep only gRNAs targeting that best ENSG_ID.
# ---------------------------------------------------------
echo "Step 2: Selecting best hit based on probability..."

awk -F',' '
BEGIN { OFS="," }
# Pass 1: Find max probability and best ENSG for each Target Gene
NR==FNR {
    if (NR==1) next
    
    target_gene = $3
    ensg_id = $6
    prob = $8
    
    # Skip if prob is NA
    if (prob == "NA") next
    
    # Convert to number
    prob_num = prob + 0
    
    if (!(target_gene in max_prob) || prob_num > max_prob[target_gene]) {
        max_prob[target_gene] = prob_num
        best_ensg[target_gene] = ensg_id
    }
    next
}
# Pass 2: Filter lines
{
    if (FNR==1) { print $0; next }
    
    target_gene = $3
    ensg_id = $6
    
    # If we found a best match for this target gene
    if (target_gene in best_ensg) {
        # Only print if this line matches the best ENSG ID
        if (ensg_id == best_ensg[target_gene]) {
            print $0
        }
    } else {
        # If no probability info was found for this target gene at all,
        # we keep the entries (or should we filter them? Liang script keeps one NA representative, 
        # but here we might want to keep all if no info is available, or follow strict logic.
        # The Liang script "select_best_hit.sh" keeps one NA line if no match found.
        # However, in our unified context, we might have multiple potential targets.
        # Let"s keep all if no probability is available to avoid data loss, 
        # or we can refine this later.
        print $0
    }
}
' "$TEMP_UNIQUE" "$TEMP_UNIQUE" > "$TEMP_BEST"

# Move final result
mv "$TEMP_BEST" "$OUTPUT_FILE"
rm "$TEMP_UNIQUE"

echo "Filtering complete. Result saved to $OUTPUT_FILE"
