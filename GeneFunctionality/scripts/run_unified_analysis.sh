#!/bin/bash
set -ueo pipefail

# Activate environment
# Note: This assumes micromamba is initialized in the shell. 
# If running from a non-interactive shell, we might need to source the profile.
# Trying direct activation first as requested.
eval "$(micromamba shell hook --shell bash)"
micromamba activate base

# Define paths
DATA_DIR="../data"
UNIFIED_DIR="$DATA_DIR/unified"
UNIFIED_CSV="$UNIFIED_DIR/unified_gRNAs.csv"
UNIFIED_FASTA="$UNIFIED_DIR/unified_gRNAs.fasta"
MODEL_DIR="$DATA_DIR/model_predictions"
EXON1_FASTA="$MODEL_DIR/lncrna_exon1.fasta"
EXON2_FASTA="$MODEL_DIR/lncrna_exon2.fasta"
MODEL_PREDICTIONS="$MODEL_DIR/gencode-lncrna-ranking.csv"
RESULTS_DIR="../results"
REF_DIR="$DATA_DIR/references"
GENOME_FASTA="$REF_DIR/gencode.v49.lncRNA_transcripts.fa"

mkdir -p "$RESULTS_DIR"

# Step 1: Convert unified CSV to FASTA
echo "Converting unified CSV to FASTA..."
# Header format: Study|gRNA_ID|Target_Gene_ID|gRNA_Type|Essentiality
awk -F',' 'NR>1 {
    # Sanitize fields to remove pipes if any
    gsub(/\|/, "_", $1)
    gsub(/\|/, "_", $2)
    gsub(/\|/, "_", $3)
    gsub(/\|/, "_", $5)
    gsub(/\|/, "_", $6)
    
    header = $1 "|" $2 "|" $3 "|" $5 "|" $6
    print ">" header
    print $4
}' "$UNIFIED_CSV" > "$UNIFIED_FASTA"

# Step 2: Run Genome-wide BLAT
echo "Running Genome-wide BLAT (this may take a while)..."
# Check if genome file exists
if [ ! -f "$GENOME_FASTA" ]; then
    echo "Error: Genome FASTA not found at $GENOME_FASTA. Cannot proceed."
    exit 1
fi

GENOME_PSL="$UNIFIED_DIR/unified_vs_genome.psl"
blat -t=dna -q=dna "$GENOME_FASTA" "$UNIFIED_FASTA" -minScore=15 -minIdentity=100 "$GENOME_PSL".tmp

# Keep only antisense blat matches
awk '$9 ~ /^-/' "$GENOME_PSL".tmp > "$GENOME_PSL"
rm "$GENOME_PSL".tmp

# Step 3: Run Protein BLAT and Off-Target Check
echo "Running Protein BLAT and Off-Target Check..."
PROTEIN_FASTA="$REF_DIR/human_proteins.fasta"
GTF_FILE="$REF_DIR/gencode.v49.long_noncoding_RNAs.gtf"
PROTEIN_PSL="$UNIFIED_DIR/unified_vs_proteins.psl"
PROTEIN_MATCHES="$UNIFIED_DIR/protein_matches.csv"

if [ ! -f "$PROTEIN_FASTA" ]; then
    echo "Warning: Protein FASTA not found at $PROTEIN_FASTA. Skipping protein check."
    # Create empty matches file to avoid errors downstream
    touch "$PROTEIN_MATCHES"
else
    blat -t=dna -q=dna "$PROTEIN_FASTA" "$UNIFIED_FASTA" -minScore=15 -minIdentity=100 "$PROTEIN_PSL".tmp
    awk '$9 ~ /^-/' "$PROTEIN_PSL".tmp > "$PROTEIN_PSL"
    rm "$PROTEIN_PSL".tmp
    
    echo "Checking protein matches using genome-wide BLAT results..."
    #./check_protein_matches.sh "$PROTEIN_PSL" "$GTF_FILE" > "$PROTEIN_MATCHES"
    ./join_psl.sh "$PROTEIN_PSL" "$GENOME_PSL" "$PROTEIN_MATCHES"
fi

# Step 4: Process genome-wide BLAT results and add protein off-target and model predictions
echo "Processing genome-wide BLAT results, adding protein off-target annotation and model predictions..."


awk -F'\t' -v protein_file="$PROTEIN_MATCHES" -v predictions="$MODEL_PREDICTIONS" -v OFS=',' '
BEGIN {
    # ---------------------------------------------------------
    # 1. Load Protein Matches
    # Key = Column 1 (Full ID string)
    # Value = Column 5 (YES/NO Result)
    # ---------------------------------------------------------
    FS=","
    while ((getline < protein_file) > 0) {
        if ($1 == "Key") continue
        split($0, arr, ",")
        protein_matches[arr[1]] = arr[5]
    }
    close(protein_file)
    
    # ---------------------------------------------------------
    # 2. Load Model Predictions
    # Key = Base ENSG ID
    # Value = Column 6 (highest_prob)
    # ---------------------------------------------------------
    while ((getline < predictions) > 0) {
        if ($1 == "GeneID") continue
        
        gene_id = $1
        
        # [MODIFIED] Now reading from Column 6 (highest_prob)
        # Previous logic used $4 (Probability_Functional)
        prob_val = $6
        
        # Store prediction (strip version from gene_id key)
        split(gene_id, id_parts, ".")
        base_gene_id = id_parts[1]
        
        if (base_gene_id != "" && prob_val != "") {
            predictions_map[base_gene_id] = prob_val
        }
    }
    close(predictions)
    
    # Print Output Header
    print "Study,gRNA_ID,Target_Gene_ID,gRNA_Type,Essentiality,ENSG_ID,Strand,Probability_Functional,Protein_Off_Target"
    
    # Reset Field Separator to Tab for the PSL file
    FS="\t"
}

# ---------------------------------------------------------
# 3. Process Genome PSL (Main Loop)
# ---------------------------------------------------------
FNR > 5 {
    q_name = $10
    
    # Parse Q_name
    split(q_name, q_parts, "|")
    study = q_parts[1]
    grna_id = q_parts[2]
    target_gene = q_parts[3]
    type = q_parts[4]
    essentiality = q_parts[5]
    
    # Parse T_name (Col 14)
    t_name = $14
    split(t_name, t_parts, "|")
    ensg_with_version = t_parts[2]
    
    # Extract strand
    strand_field = t_parts[length(t_parts)]
    if (strand_field ~ /\(\+\)/) {
        strand = "+"
    } else if (strand_field ~ /\(-\)/) {
        strand = "-"
    } else {
        strand = "NA"
    }
    
    # Lookup 1: Protein Match Result
    prot_match = "NO"
    if (q_name in protein_matches) {
        prot_match = protein_matches[q_name]
    }
    
    # Lookup 2: Model Prediction
    split(ensg_with_version, ensg_parts, ".")
    base_ensg = ensg_parts[1]
    
    prob_func = "NA"
    if (base_ensg in predictions_map) {
        prob_func = predictions_map[base_ensg]
    }
    
    print study, grna_id, target_gene, type, essentiality, ensg_with_version, strand, prob_func, prot_match
}
' "$GENOME_PSL" > "$RESULTS_DIR/unified_genome_alignments_raw.csv"

# Step 4.5: Apply Liang Filters (Unique Hits & Best Hit Selection)
echo "Applying Liang filters (Unique Hits & Best Hit Selection)..."
./apply_liang_filters.sh "$RESULTS_DIR/unified_genome_alignments_raw.csv" "$RESULTS_DIR/unified_genome_alignments.csv"

# Step 4.6: Append Unmatched gRNAs
echo "Appending unmatched gRNAs..."
awk -F',' '
BEGIN { OFS="," }
# 1. Read the filtered alignments (File 1)
# Store seen gRNAs
NR==FNR {
    # Key: Study|gRNA_ID|Target_Gene_ID
    # Columns: 1:Study, 2:gRNA_ID, 3:Target_Gene_ID
    key = $1 "|" $2 "|" $3
    seen[key] = 1
    next
}
# 2. Read the original unified CSV (File 2)
# Columns: 1:Study, 2:gRNA_ID, 3:Target_Gene_ID, 4:Sequence, 5:gRNA_Type, 6:Essentiality
FNR > 1 {
    key = $1 "|" $2 "|" $3
    
    if (!(key in seen)) {
        # Construct the line with NAs
        # Output format: Study,gRNA_ID,Target_Gene_ID,gRNA_Type,Essentiality,ENSG_ID,Strand,Probability_Functional,Protein_Off_Target
        # Input cols: 1, 2, 3, 5, 6
        print $1, $2, $3, $5, $6, "NA", "NA", "NA", "NA"
    }
}
' "$RESULTS_DIR/unified_genome_alignments.csv" "$UNIFIED_CSV" >> "$RESULTS_DIR/unified_genome_alignments.csv"

# Step 5: Add gRNA sequences from FASTA
echo "Adding gRNA sequences to results..."

# Extract gRNA sequences from FASTA and create a lookup map
awk '
/^>/ {
    # Parse header: >Study|gRNA_ID|Target_Gene_ID|gRNA_Type|Essentiality
    sub(/^>/, "")
    split($0, hdr, "|")
    grna_id = hdr[2]
    getline
    seq = $0
    print grna_id "," seq
}
' "$UNIFIED_FASTA" > "$RESULTS_DIR/grna_seq_map.csv"

# Merge sequences into the results
awk -F',' 'BEGIN {OFS=","}
NR==FNR {
    # Load sequence map
    seq_map[$1] = $2
    next
}
FNR==1 {
    # Print header with new column
    print $0, "gRNA_Sequence"
    next
}
{
    # Add sequence for each row
    grna_id = $2
    seq = (grna_id in seq_map) ? seq_map[grna_id] : "NA"
    print $0, seq
}
' "$RESULTS_DIR/grna_seq_map.csv" "$RESULTS_DIR/unified_genome_alignments.csv" > "$RESULTS_DIR/unified_genome_alignments_with_seq.csv"

# Replace original with new version
mv "$RESULTS_DIR/unified_genome_alignments_with_seq.csv" "$RESULTS_DIR/unified_genome_alignments.csv"
rm "$RESULTS_DIR/grna_seq_map.csv" "$RESULTS_DIR/unified_genome_alignments_raw.csv"

echo "Analysis complete. Results in $RESULTS_DIR/unified_genome_alignments.csv"
