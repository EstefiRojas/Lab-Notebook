#!/bin/bash
# ==============================================================================
# Script: find_lncRNA_guides.sh
# Description: This script aligns gRNAs against the full human genome (GRCh38)
#              and then intersects the results with a dedicated GENCODE lncRNA
#              annotation file. The final report includes ALL BLAT hits; if a
#              hit does not overlap a lncRNA, "NA" is reported in lncRNA fields.
#              The report also includes the gRNA sequence and the specific
#              feature type (gene, transcript) of any lncRNA that was hit.
# Date: September 23, 2025
#
# Prerequisites:
#   - BLAT: (https://genome.ucsc.edu/goldenPath/help/blatSpec.html)
#   - bedtools: (https://bedtools.readthedocs.io/en/latest/)
#   - wget, gunzip, awk
#
# Usage:
#   1. Run this script from the project's root directory.
#   2. Make the script executable: chmod +x scripts/find_lncRNA_guides.sh
#   3. Run the script with your two CSV files as arguments:
#      ./scripts/find_lncRNA_guides.sh path/to/gRNA_info.csv path/to/essentiality.csv
#
# ==============================================================================

# --- Script Configuration ---

# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error when substituting.
set -u

# --- Project Directory Structure ---
echo "--- Setting up project directories ---"
DATA_DIR="../data"
REF_DIR="$DATA_DIR/references"
RESULTS_DIR="../results"
TMP_DIR="../tmp"

mkdir -p "$REF_DIR" "$RESULTS_DIR" "$TMP_DIR"
echo "Directories are ready."
echo ""

# --- URLs for Genome and lncRNA Annotation ---
GENOME_FASTA_URL="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
LNCRNA_GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.long_noncoding_RNAs.gtf.gz"

# --- Filenames (with directory paths) ---
GENOME_FASTA_NAME="GRCh38.primary_assembly.genome.fa"
LNCRNA_GTF_NAME="gencode.v49.long_noncoding_RNAs.gtf"

GENOME_FASTA="$REF_DIR/$GENOME_FASTA_NAME"
LNCRNA_GTF="$REF_DIR/$LNCRNA_GTF_NAME"

LNCRNA_BED="$TMP_DIR/lncrna_features.bed"
FILTERED_GUIDES_FASTA="$TMP_DIR/filtered_guides.fa"
GRNA_METADATA_MAP="$TMP_DIR/gRNA_metadata_map.tsv"
BLAT_OUTPUT_PSL="$TMP_DIR/blat_output.psl"
BLAT_OUTPUT_BED="$TMP_DIR/blat_output.bed"
TEMP_INTERSECT="$TMP_DIR/temp_intersect.txt"
FINAL_REPORT="$RESULTS_DIR/gRNA_lncRNA_matches.tsv"

# ==============================================================================
# --- STEP 0: Check for prerequisites and input files ---
# ==============================================================================
echo "--- Step 0: Checking prerequisites ---"

if [ -z "${2-}" ]; then
    echo "ERROR: Two input CSV files are required."
    echo "Usage: ./scripts/find_lncRNA_guides.sh <gRNA_info.csv> <essentiality_info.csv>"
    exit 1
fi
INPUT_GRNA_CSV="$1"
INPUT_ESSENTIALITY_CSV="$2"

if ! command -v blat &> /dev/null; then
    echo "ERROR: BLAT could not be found. Please install it and ensure it's in your PATH."
    exit 1
fi

if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools could not be found. Please install it."
    exit 1
fi

if [ ! -f "$INPUT_GRNA_CSV" ]; then
    echo "ERROR: gRNA input file '$INPUT_GRNA_CSV' not found."
    exit 1
fi
if [ ! -f "$INPUT_ESSENTIALITY_CSV" ]; then
    echo "ERROR: Essentiality input file '$INPUT_ESSENTIALITY_CSV' not found."
    exit 1
fi

echo "Prerequisites met. Starting workflow..."
echo ""


# ==============================================================================
# --- STEP 1: Download and Prepare Reference Files ---
# ==============================================================================
echo "--- Step 1: Downloading and Preparing Reference Files ---"

if [ ! -f "$GENOME_FASTA" ]; then
    echo "Downloading human genome GRCh38 to $REF_DIR..."
    wget -O ${GENOME_FASTA}.gz "$GENOME_FASTA_URL"
    echo "Decompressing genome..."
    gunzip ${GENOME_FASTA}.gz
else
    echo "Genome FASTA file '$GENOME_FASTA' already exists. Skipping download."
fi

if [ ! -f "$LNCRNA_GTF" ]; then
    echo "Downloading lncRNA annotations to $REF_DIR..."
    wget -O ${LNCRNA_GTF}.gz "$LNCRNA_GTF_URL"
    echo "Decompressing annotations..."
    gunzip ${LNCRNA_GTF}.gz
else
    echo "lncRNA GTF file '$LNCRNA_GTF' already exists. Skipping download."
fi
echo "File preparation complete."
echo ""


# ==============================================================================
# --- STEP 2: Create a BED file of lncRNA Features ---
# ==============================================================================
echo "--- Step 2: Creating a BED file of lncRNA features ---"

# Parse only gene or transcript features from the GTF.
awk 'BEGIN{OFS="\t"} $3 == "gene" {
    attributes = $0;
    sub(/.*gene_id "/, "", attributes);
    sub(/".*/, "", attributes);
    # chrom, start, end, gene_id, score, strand, feature_type
    print $1, $4-1, $5, attributes, "0", $7, $3
}' "$LNCRNA_GTF" > "$LNCRNA_BED"

echo "Created BED file for lncRNA features: '$LNCRNA_BED'"
LNCRNA_GENE_COUNT=`cut -f 4 "$LNCRNA_BED" | sort | uniq | wc -l`
echo "$LNCRNA_GENE_COUNT unique lncRNA genes found in annotation."
echo ""


# ==============================================================================
# --- STEP 3: Filter gRNAs and Create FASTA File ---
# ==============================================================================
echo "--- Step 3: Filtering gRNAs based on essentiality and creating FASTA file ---"

rm -f "$FILTERED_GUIDES_FASTA" "$GRNA_METADATA_MAP"
touch "$FILTERED_GUIDES_FASTA"
touch "$GRNA_METADATA_MAP"

awk 'BEGIN {FS=","; OFS="\t"}
FNR==NR {
    sub(/\r$/, "", $2);
    if (FNR > 1 && ($2 == "Cell-type specific" || $2 == "Partially shared" || $2 == "Shared")) {
        essential_genes[$1] = $2
    }
    next
}
{
    sub(/\r$/, "", $4);
    if (FNR > 1 && ($2 in essential_genes)) {
        print ">" $1 "\n" $4 >> "'"$FILTERED_GUIDES_FASTA"'"
        # Add the sequence ($4) as the last column to the map
        print $1, essential_genes[$2], $2, $4 >> "'"$GRNA_METADATA_MAP"'"
    }
}' "$INPUT_ESSENTIALITY_CSV" "$INPUT_GRNA_CSV"

FILTERED_COUNT=`grep -c "^>" "$FILTERED_GUIDES_FASTA"`
echo "Created FASTA file '$FILTERED_GUIDES_FASTA' with $FILTERED_COUNT filtered gRNAs."
echo ""

# ==============================================================================
# --- STEP 4: Run BLAT against the Whole Genome ---
# ==============================================================================
echo "--- Step 4: Running BLAT against the Whole Genome ---"

if [ "$FILTERED_COUNT" -eq 0 ]; then
    echo "No gRNAs passed the essentiality filter. Exiting."
    touch "$FINAL_REPORT"
    exit 0
fi

echo "Aligning gRNAs to the genome with BLAT (this may take a while)..."
blat -t=dna -q=dna "$GENOME_FASTA" "$FILTERED_GUIDES_FASTA" -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 "$BLAT_OUTPUT_PSL" -noHead
echo "BLAT alignment complete. Output saved to '$BLAT_OUTPUT_PSL'"
echo ""


# ==============================================================================
# --- STEP 5: Process BLAT results and find overlaps ---
# ==============================================================================
echo "--- Step 5: Processing results and finding overlaps ---"

echo "Converting BLAT output to BED format..."
awk 'BEGIN{OFS="\t"} {print $14, $16, $17, $10, $1, $9}' "$BLAT_OUTPUT_PSL" > "$BLAT_OUTPUT_BED"

echo "Intersecting gRNA alignments with lncRNA regions (reporting all alignments)..."
# Use -loj for a "left outer join" to keep all BLAT hits, even if they dont overlap a lncRNA.
bedtools intersect -a "$BLAT_OUTPUT_BED" -b "$LNCRNA_BED" -loj > "$TEMP_INTERSECT"

echo "Formatting final report..."
echo -e "gRNA_ID\tEssentiality\tTarget_Gene_ID\tgRNA_Sequence\tgRNA_Chr\tgRNA_Start\tgRNA_End\tlncRNA_Feature_Type\tlncRNA_ENSG_ID\tlncRNA_Chr\tlncRNA_Start\tlncRNA_End" > "$FINAL_REPORT"

awk 'BEGIN{FS="\t"; OFS="\t"}
# 1. Read gRNA metadata map
FNR==NR {
    essential_map[$1] = $2;
    target_gene_map[$1] = $3;
    sequence_map[$1] = $4;
    next;
}
# 2. Process the intersected results from the left outer join
{
    gRNA_id = $4;
    gRNA_chr = $1;
    gRNA_start = $2;
    gRNA_end = $3;

    # Check the 7th field of the bedtools output. If it is ".", there was no match.
    if ($7 == ".") {
        lncRNA_feature = "NA";
        lncRNA_ensg = "NA";
        lncRNA_chr = "NA";
        lncRNA_start = "NA";
        lncRNA_end = "NA";
    } else {
        # A match was found, extract the lncRNA info
        lncRNA_feature = $13;
        lncRNA_ensg = $10;
        lncRNA_chr = $7;
        lncRNA_start = $8;
        lncRNA_end = $9;
    }
    
    print gRNA_id, essential_map[gRNA_id], target_gene_map[gRNA_id], sequence_map[gRNA_id], gRNA_chr, gRNA_start, gRNA_end, lncRNA_feature, lncRNA_ensg, lncRNA_chr, lncRNA_start, lncRNA_end
}' "$GRNA_METADATA_MAP" "$TEMP_INTERSECT" | sort | uniq >> "$FINAL_REPORT"

# Clean up intermediate files from the tmp directory
#rm "$TMP_DIR"/*

echo ""
echo "--- Workflow Complete! ---"
echo "Final report is saved in: '$FINAL_REPORT'"

