#!/bin/bash
set -ueo pipefail

# Define paths
DATA_DIR="../data"
RESULTS_DIR="../results"
INPUT_FILE="$RESULTS_DIR/annotated_unified_genome_alignments.csv"
OUTPUT_FILE="$RESULTS_DIR/essentiality_matrix.csv"

echo "Generating essentiality matrix..."
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found at $INPUT_FILE"
    exit 1
fi

awk -F',' '
BEGIN {
    OFS=","
    # Define studies for consistent column order
    studies[1] = "Huang"
    studies[2] = "Liang"
    studies[3] = "Liu"
    studies[4] = "Montero"
}
NR > 1 {
    # Clean carriage returns
    gsub(/\r/, "", $0)

    # Columns based on annotated_unified_genome_alignments.csv structure:
    # 1: Study
    # 2: gRNA_ID
    # 3: Target_Gene_ID
    # 4: gRNA_Type
    # 5: Essentiality
    # 6: ENSG_ID
    # 7: Strand
    # 8: Probability_Functional
    # 9: Protein_Off_Target
    # 10: Gene_Name (Added by map_lncRNA_genes.sh)
    # 11: gRNA_Sequence (Added by run_unified_analysis.sh)
    
    study = $1
    essentiality = $5
    ensg_id = $6
    prob = $8
    gene_name = $11
    
    # Filter:
    # 1. Essentiality is Rare, Common, or Core
    # 2. Probability_Functional is NOT NA
    if (prob != "NA" && (essentiality == "Rare" || essentiality == "Common" || essentiality == "Core")) {
        
        # Store data
        # If multiple gRNAs map to the same gene in the same study with different essentialities,
        # we need a policy. For now, lets assume we overwrite or keep the "strongest".
        # However, the prompt implies a single status per gene/study.
        # If duplicates exist, this simple assignment takes the last one seen.
        
        # Check if we already have a value for this gene/study
        if ((ensg_id, study) in matrix) {
            current_val = matrix[ensg_id, study]
            # Priority logic (optional): Core > Common > Rare
            # If current is Rare and new is Common, update.
            # If current is Common and new is Core, update.
            if (essentiality == "Core") {
                matrix[ensg_id, study] = essentiality
            } else if (essentiality == "Common" && current_val == "Rare") {
                matrix[ensg_id, study] = essentiality
            }
            # Else keep current
        } else {
            matrix[ensg_id, study] = essentiality
        }
        
        unique_genes[ensg_id] = 1
        gene_probs[ensg_id] = prob
        gene_names[ensg_id] = gene_name
    }
}
END {
    # Print Header
    printf "ENSG_ID,Gene_Name,Probability_Functional"
    for (i=1; i<=4; i++) {
        printf ",%s", studies[i]
    }
    printf "\n"
    
    # Print Rows
    for (gene in unique_genes) {
        printf "%s,%s,%s", gene, gene_names[gene], gene_probs[gene]
        for (i=1; i<=4; i++) {
            study = studies[i]
            val = "-"
            if ((gene, study) in matrix) {
                val = matrix[gene, study]
            }
            printf ",%s", val
        }
        printf "\n"
    }
}
' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Done."
