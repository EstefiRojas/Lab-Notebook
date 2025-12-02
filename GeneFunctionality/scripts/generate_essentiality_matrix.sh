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

    # Columns based on NEW annotated_unified_genome_alignments.csv structure:
    # 1: Study
    # 2: gRNA_ID
    # 3: gRNA_Sequence
    # 4: gRNA_Type
    # 5: Target_Gene_ID
    # 6: ENSG_ID
    # 7: Gene_Name
    # 8: Strand
    # 9: Essentiality
    # 10: Probability_Functional
    # 11: Protein_Off_Target
    # 12: Antisense_to_CDS
    # 13: Antisense_Target_Name
    
    study = $1
    essentiality = $9
    ensg_id = $6
    prob = $10
    gene_name = $7
    
    # Filter:
    # 1. Essentiality is Rare, Common, Core, or Non-essential
    # 2. Probability_Functional is NOT NA
    if (prob != "NA" && (essentiality == "Rare" || essentiality == "Common" || essentiality == "Core" || essentiality == "Non-essential")) {
        
        # Store data
        # If multiple gRNAs map to the same gene in the same study with different essentialities,
        # we need a policy. For now, lets assume we overwrite or keep the "strongest".
        # However, the prompt implies a single status per gene/study.
        # If duplicates exist, this simple assignment takes the last one seen.
        
        # Check if we already have a value for this gene/study
        if ((ensg_id, study) in matrix) {
            current_val = matrix[ensg_id, study]
            # Priority logic: Core > Common > Rare > Non-essential
            
            # Helper function logic (simulated):
            # 3 = Core, 2 = Common, 1 = Rare, 0 = Non-essential
            
            new_score = 0
            if (essentiality == "Rare") new_score = 1
            if (essentiality == "Common") new_score = 2
            if (essentiality == "Core") new_score = 3
            
            curr_score = 0
            if (current_val == "Rare") curr_score = 1
            if (current_val == "Common") curr_score = 2
            if (current_val == "Core") curr_score = 3
            
            if (new_score > curr_score) {
                matrix[ensg_id, study] = essentiality
            }
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
