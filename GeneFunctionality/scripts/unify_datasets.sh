#!/bin/bash
set -ueo pipefail

# Define paths
DATA_DIR="../data"
UNIFIED_DIR="$DATA_DIR/unified"
UNIFIED_FILE="$UNIFIED_DIR/unified_gRNAs.csv"

# Create unified directory
mkdir -p "$UNIFIED_DIR"

# Initialize unified file with header
echo "Study,gRNA_ID,Target_Gene_ID,Sequence,gRNA_Type,Essentiality" > "$UNIFIED_FILE"

echo "Processing Huang et al..."
# Huang: SuppT6b.csv
# ID: Col 1 (crRNA_A_location) / Col 2 (crRNA_B_location)
# Gene: Col 7 (gene ID)
# Seq: Col 3 (crRNA_A) / Col 4 (crRNA_B)
# Type: Col 5 (gene type)
# Filter: Keep only if Type == "lncRNA"
# Essentiality: Check data/Huang/column_matches_ensgids.csv (Cols 2,3,4 count)
awk -F',' '
BEGIN {
    # Load Essentiality Map
    while ((getline < "'"$DATA_DIR/Huang/column_matches_ensgids.csv"'") > 0) {
        # Skip header if needed, but logic below handles empty checks
        # Col 1 is GeneID
        # Col 2,3,4 are cell lines
        count = 0
        if ($2 != "") count++
        if ($3 != "") count++
        if ($4 != "") count++
        
        if (count >= 2) {
            ess_map[$1] = "Common"
        } else if (count == 1) {
            ess_map[$1] = "Specific"
        } else {
            ess_map[$1] = "Non-essential"
        }
    }
    close("'"$DATA_DIR/Huang/column_matches_ensgids.csv"'")
}
NR>1 {
    # Clean carriage returns
    gsub(/\r/, "", $0)
    
    if ($5 == "lncRNA") {
        gene_id = $7
        essentiality = "Non-essential"
        if (gene_id in ess_map) {
            essentiality = ess_map[gene_id]
        }
        
        # gRNA A
        if ($3 != "") {
            print "Huang," $1 "," $7 "," $3 "," $5 "," essentiality >> "'"$UNIFIED_FILE"'"
        }
        # gRNA B
        if ($4 != "") {
            print "Huang," $2 "," $7 "," $4 "," $5 "," essentiality >> "'"$UNIFIED_FILE"'"
        }
    }
}' "$DATA_DIR/Huang/SuppT6b.csv"

echo "Processing Liang et al..."
# Liang: LiangMuller_May2025_Table1_Filtered-gRNAs.csv
# ID: Col 1 (gRNA_ID)
# Gene: Col 2 (Target_Gene_ID)
# Seq: Col 4 (Sequence)
# Type: Col 3 (Type)
# Filter: Keep only if Type == "lncRNA"
# Essentiality: Map Group from LiangMuller_May2025-Table3-Gene-RRA-Ranking.csv
awk -F',' '
BEGIN {
    # Load Essentiality Map
    while ((getline < "'"$DATA_DIR/Liang/LiangMuller_May2025-Table3-Gene-RRA-Ranking.csv"'") > 0) {
        # Col 1: Gene
        # Col 2: Group
        group = $2
        gsub(/\r/, "", group)
        if (group == "Shared") {
            ess_map[$1] = "Core"
        } else if (group == "Partially shared") {
            ess_map[$1] = "Common"
        } else if (group == "Cell-type specific") {
            ess_map[$1] = "Specific"
        } else {
            ess_map[$1] = "Non-essential"
        }
    }
    close("'"$DATA_DIR/Liang/LiangMuller_May2025-Table3-Gene-RRA-Ranking.csv"'")
}
NR>1 {
    gsub(/\r/, "", $0)
    if ($3 == "lncRNA" && $4 != "") {
        gene_id = $2
        essentiality = "Non-essential"
        if (gene_id in ess_map) {
            essentiality = ess_map[gene_id]
        }
        print "Liang," $1 "," $2 "," $4 "," $3 "," essentiality >> "'"$UNIFIED_FILE"'"
    }
}' "$DATA_DIR/Liang/LiangMuller_May2025_Table1_Filtered-gRNAs.csv"

echo "Processing Liu et al..."
# Liu: TableS2_gRNA_sequences.csv
# ID: Col 1 (Gene ID)
# Gene: Col 1 (Gene ID)
# Seq: Col 2 (Protospacer sequence)
# Type: "NA" -> "lncRNA"
# Filter: Keep only if Gene ID starts with "LH"
# Essentiality: Count occurrences in Essential_lncrnas_cell_lines_IDs.csv
awk -F',' '
BEGIN {
    # Load Essentiality Map
    while ((getline < "'"$DATA_DIR/Liu/Essential_lncrnas_cell_lines_IDs.csv"'") > 0) {
        # Iterate over all 7 columns
        for (i=1; i<=7; i++) {
            id = $i
            gsub(/\r/, "", id)
            if (id != "") {
                counts[id]++
            }
        }
    }
    close("'"$DATA_DIR/Liu/Essential_lncrnas_cell_lines_IDs.csv"'")
}
NR>1 {
    gsub(/\r/, "", $0)
    if ($1 ~ /^LH/ && $2 != "") {
        gene_id = $1
        count = 0
        if (gene_id in counts) {
            count = counts[gene_id]
        }
        
        if (count == 7) {
            essentiality = "Core"
        } else if (count >= 2) {
            essentiality = "Common"
        } else if (count == 1) {
            essentiality = "Specific"
        } else {
            essentiality = "Non-essential"
        }
        
        print "Liu," $1 "," $1 "," $2 ",lncRNA," essentiality >> "'"$UNIFIED_FILE"'"
    }
}' "$DATA_DIR/Liu/TableS2_gRNA_sequences.csv"

echo "Processing Montero et al..."
# Montero: Supplementary_table9_useThis.csv
# ID: Col 2 (gRNA_ID)
# Gene: Col 1 (Target_Gene_ID)
# Seq: Col 4 (gRNA1_Seq) / Col 5 (gRNA2_Seq)
# Type: Col 3 (Type)
# Filter: Keep only if Type contains "lncRNA"
# Essentiality: Map group from Supplementary_table13_new_labels.csv (Col 8)
awk -F',' '
BEGIN {
    # Load Essentiality Map
    while ((getline < "'"$DATA_DIR/Montero/processed/Supplementary_table13_new_labels.csv"'") > 0) {
        # Col 1: target (Gene ID)
        # Col 8: group
        # Need to handle potential CSV parsing issues if there are quotes, but assuming simple CSV for now.
        gene = $1
        group = $8
        gsub(/\r/, "", group)
        
        # Normalize group case if needed? User said "use the group column".
        # Assuming values like "common", "specific".
        # Let"s capitalize first letter to match others: Common, Specific.
        if (tolower(group) == "common") {
            ess_map[gene] = "Common"
        } else if (tolower(group) == "specific") {
            ess_map[gene] = "Specific"
        } else if (tolower(group) == "core") {
             ess_map[gene] = "Core"
        } else {
            ess_map[gene] = "Non-essential"
        }
    }
    close("'"$DATA_DIR/Montero/processed/Supplementary_table13_new_labels.csv"'")
}
NR>1 {
    gsub(/\r/, "", $0)
    if ($3 ~ /lncRNA/) {
        gene_id = $1
        essentiality = "Non-essential"
        if (gene_id in ess_map) {
            essentiality = ess_map[gene_id]
        }
        
        # gRNA 1
        if ($4 != "") {
            print "Montero," $2 "_1," $1 "," $4 "," $3 "," essentiality >> "'"$UNIFIED_FILE"'"
        }
        # gRNA 2
        if ($5 != "") {
            print "Montero," $2 "_2," $1 "," $5 "," $3 "," essentiality >> "'"$UNIFIED_FILE"'"
        }
    }
}' "$DATA_DIR/Montero/Supplementary_table9_useThis.csv"

echo "Processing Zhu et al..."
# Zhu: Supp_Table10.csv
# ID: Col 1 (ID)
# Gene: Col 2 (Target_gene)
# Seq: Col 8 (sgRNA1_sequence) / Col 13 (sgRNA2_sequence)
# Type: Col 15 (Type)
# Filter: Exclude if Type == "HOTAIR"
# Essentiality: All "Non-essential"
awk -F',' 'NR>1 {
    gsub(/\r/, "", $0)
    if ($15 != "HOTAIR") {
        essentiality = "Non-essential"
        # gRNA 1
        if ($8 != "") {
            print "Zhu," $1 "_1," $2 "," $8 "," $15 "," essentiality >> "'"$UNIFIED_FILE"'"
        }
        # gRNA 2
        if ($13 != "") {
            print "Zhu," $1 "_2," $2 "," $13 "," $15 "," essentiality >> "'"$UNIFIED_FILE"'"
        }
    }
}' "$DATA_DIR/Zhu/Supp_Table10.csv"

echo "Unified dataset created at $UNIFIED_FILE"
# Count rows
wc -l "$UNIFIED_FILE"
