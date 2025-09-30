#!/bin/bash

# ==============================================================================
# Script: get_exon_sequences.sh
# Description: This script takes a list of gene IDs and extracts the DNA
#              sequences for their first two exons. It downloads the human
#              genome and gene annotations, identifies the coordinates for
#              exons 1 and 2 for each specified gene, and then uses bedtools
#              to extract the corresponding sequences into a FASTA file.
# Date: September 23, 2025
#
# Prerequisites:
#    - bedtools: For extracting sequences.
#    - samtools: For indexing the genome FASTA file.
#    - Standard tools: wget, gunzip, awk, grep, sort, cut, tail.
#
# Usage:
#    1. Make the script executable: chmod +x scripts/get_exon_sequences.sh
#    2. Run the script with the filtered report as input:
#       ./scripts/get_exon_sequences.sh results/gRNA_lncRNA_matches_unique_ensg_na_prob.tsv
#
# ==============================================================================

# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error when substituting.
set -u

# --- STEP 0: Check for input file and dependencies ---
if [ -z "${1-}" ]; then
    echo "ERROR: An input file is required."
    echo "Usage: ./scripts/get_exon_sequences.sh <path_to_filtered_report.tsv>"
    exit 1
fi
INPUT_FILE="$1"

if ! command -v bedtools &> /dev/null || ! command -v samtools &> /dev/null; then
    echo "ERROR: 'bedtools' and 'samtools' are required but not found in your PATH."
    echo "Please install them before running this script."
    exit 1
fi

# --- STEP 1: Define paths and create directories ---
DATA_DIR="../data"
RESULTS_DIR="../results"
TMP_DIR="../tmp"
mkdir -p "$DATA_DIR" "$RESULTS_DIR" "$TMP_DIR"

# Define URLs and local file paths
GENOME_URL="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz"

GENOME_FASTA_GZ="$DATA_DIR/references/GRCh38.primary_assembly.genome.fa.gz"
GENOME_FASTA="$DATA_DIR/references/GRCh38.primary_assembly.genome.fa"
GTF_FILE="$DATA_DIR/references/gencode.v49.primary_assembly.annotation.gtf.gz"
GTF_UNCOMPRESSED="$DATA_DIR/references/gencode.v49.primary_assembly.annotation.gtf"

TARGET_IDS_FILE="$TMP_DIR/target_ensg_ids.txt"
ALL_EXONS_GTF="$TMP_DIR/all_target_exons.gtf"
EXONS_BED_FILE="$TMP_DIR/exons_to_extract.bed"
OUTPUT_FASTA="$TMP_DIR/exon_1_2_sequences.fa"
OUTPUT_CSV_EXON1="$RESULTS_DIR/exon_1_sequences.csv"
OUTPUT_CSV_EXON2="$RESULTS_DIR/exon_2_sequences.csv"

echo "--- STEP 2: Downloading and preparing reference files (if necessary) ---"

# Download and unzip genome FASTA
if [ ! -f "$GENOME_FASTA" ]; then
    echo "Genome FASTA not found. Downloading..."
    wget -O "$GENOME_FASTA_GZ" "$GENOME_URL"
    echo "Decompressing genome..."
    gunzip -c "$GENOME_FASTA_GZ" > "$GENOME_FASTA"
    rm "$GENOME_FASTA_GZ"
fi

# Index genome FASTA with samtools
if [ ! -f "${GENOME_FASTA}.fai" ]; then
    echo "Indexing genome FASTA with samtools..."
    samtools faidx "$GENOME_FASTA"
fi

# Download GTF annotation
if [ ! -f "$GTF_FILE" ]; then
    echo "Annotation GTF not found. Downloading..."
    wget -O "$GTF_FILE" "$GTF_URL"
fi

# Decompress GTF for faster access if not already done
GTF_UNCOMPRESSED="$DATA_DIR/references/gencode.v49.primary_assembly.annotation.gtf"
if [ ! -f "$GTF_UNCOMPRESSED" ]; then
    echo "Decompressing GTF file for faster access (one-time operation, may take 1-2 minutes)..."
    gunzip -c "$GTF_FILE" > "$GTF_UNCOMPRESSED"
fi

echo "Reference files are ready."

# --- STEP 3: Extract unique ENSG IDs from input file ---
echo "--- Extracting target ENSG IDs from '$INPUT_FILE' ---"
tail -n +2 "$INPUT_FILE" | cut -f9 | grep -v "^NA$" | sort -u > "$TARGET_IDS_FILE"

if [ ! -s "$TARGET_IDS_FILE" ]; then
    echo "ERROR: No valid ENSG IDs found in the input file."
    exit 1
fi
NUM_GENES=$(wc -l < "$TARGET_IDS_FILE")
echo "Found $NUM_GENES unique ENSG IDs to process."

# --- STEP 4: Generate BED file of first two exons for each gene ---
echo "--- Identifying coordinates for exons 1 and 2 ---"
echo "Step 4a: Searching GTF file for target genes (should complete in under 1 minute)..."
echo "Processing $NUM_GENES gene IDs against GTF file using optimized hash lookup..."

# Debug: Check first few target IDs
echo "Sample target IDs (first 3):"
head -3 "$TARGET_IDS_FILE"

# Use awk with associative array for MUCH faster lookup (O(1) vs O(n) per line)
# This loads all gene IDs into memory once, then does a single hash lookup per GTF line
awk 'BEGIN{FS="\t"; count=0}
     NR==FNR {ids[$1]=1; id_count++; next} 
     FNR==1 {print "Loaded " id_count " gene IDs into hash table" > "/dev/stderr"}
     {
         if (match($9, /gene_id "([^"]+)"/)) {
             gene_id = substr($9, RSTART+9, RLENGTH-10);
             if (gene_id in ids && $3 == "exon" && ($9 ~ /tag "Ensembl_canonical"/ || $9 ~ /tag "MANE_Select"/)) {
                 count++;
                 if (count <= 3) print "Found match: " gene_id > "/dev/stderr";
                 print
             }
         }
     }
     END {print "Total exon lines matched: " count > "/dev/stderr"}' "$TARGET_IDS_FILE" "$GTF_UNCOMPRESSED" > "$ALL_EXONS_GTF"

NUM_EXONS=$(wc -l < "$ALL_EXONS_GTF")
echo "Step 4a: Complete. Found $NUM_EXONS total exon records for $NUM_GENES genes."

echo "Step 4b: Extracting gene IDs, exon numbers, and coordinates..."
# Process the GTF records to get exons 1 and 2 based on exon_number annotation
# Extract: chr, start, end, strand, gene_id, exon_number
awk 'BEGIN{FS="\t"; OFS="\t"} {
    if (match($9, /gene_id "([^"]+)"/)) {
        gene_id = substr($9, RSTART+9, RLENGTH-10);
        if (match($9, /exon_number ([0-9]+)/)) {
            exon_num = substr($9, RSTART+12, RLENGTH-12);
            print $1, $4, $5, $7, gene_id, exon_num
        }
    }
}' "$ALL_EXONS_GTF" > "$TMP_DIR/exons_parsed.txt"

echo "Step 4c: Filtering for exon 1 and exon 2 only..."
# Filter to keep only exon_number 1 and 2, then convert to BED format
awk 'BEGIN{OFS="\t"} ($6 == 1 || $6 == 2) {
    print $1, $2-1, $3, $5"_exon"$6, "0", $4
}' "$TMP_DIR/exons_parsed.txt" > "$EXONS_BED_FILE"

NUM_BED_ENTRIES=$(wc -l < "$EXONS_BED_FILE")
echo "Generated BED file with $NUM_BED_ENTRIES exon coordinate entries."

# --- STEP 5: Extract sequences using bedtools ---
echo "--- Extracting exon sequences into FASTA format ---"
echo "Running bedtools getfasta (this may take a minute)..."
bedtools getfasta -fi "$GENOME_FASTA" -bed "$EXONS_BED_FILE" -s -name -fo "$OUTPUT_FASTA"

NUM_SEQUENCES=$(grep -c "^>" "$OUTPUT_FASTA" || echo "0")
echo "Extracted $NUM_SEQUENCES sequences."

# --- STEP 6: Convert FASTA to CSV format (split by exon) ---
echo "--- Converting FASTA to CSV format (splitting exon 1 and exon 2) ---"
echo "ENSG_ID,start,end,sequence" > "$OUTPUT_CSV_EXON1"
echo "ENSG_ID,start,end,sequence" > "$OUTPUT_CSV_EXON2"

awk '
BEGIN {
    seq = ""
    ensg_id = ""
    start_pos = ""
    end_pos = ""
    exon_num = ""
}
/^>/ {
    # Output previous sequence if exists
    if (ensg_id != "" && seq != "") {
        output_line = ensg_id "," start_pos "," end_pos "," seq
        if (exon_num == "1") {
            print output_line >> "'"$OUTPUT_CSV_EXON1"'"
        } else if (exon_num == "2") {
            print output_line >> "'"$OUTPUT_CSV_EXON2"'"
        }
    }
    
    # Parse new header: >ENSG00000310526.1_exon1::chr1:29533-29757(-)
    header = substr($0, 2)  # Remove ">"
    
    # Split by "::"
    split(header, parts, "::")
    
    # Part 1: ENSG00000310526.1_exon1
    name_part = parts[1]
    # Remove version number (.1, .2, etc) but keep exon suffix
    gsub(/\.[0-9]+_/, "_", name_part)
    ensg_id = name_part
    
    # Extract exon number
    if (match(name_part, /_exon([0-9]+)$/)) {
        exon_str = substr(name_part, RSTART, RLENGTH)
        gsub(/_exon/, "", exon_str)
        exon_num = exon_str
    }
    
    # Part 2: chr1:29533-29757(-)
    if (length(parts) >= 2) {
        coords_part = parts[2]
        # Extract start-end
        if (match(coords_part, /[^:]+:([0-9]+)-([0-9]+)/)) {
            coord_section = substr(coords_part, RSTART, RLENGTH)
            split(coord_section, c, ":")
            split(c[2], positions, "-")
            start_pos = positions[1]
            end_pos = positions[2]
            gsub(/\(.*/, "", end_pos)  # Remove strand marker
        }
    }
    
    seq = ""
    next
}
{
    # Accumulate sequence lines
    seq = seq $0
}
END {
    # Output last sequence
    if (ensg_id != "" && seq != "") {
        output_line = ensg_id "," start_pos "," end_pos "," seq
        if (exon_num == "1") {
            print output_line >> "'"$OUTPUT_CSV_EXON1"'"
        } else if (exon_num == "2") {
            print output_line >> "'"$OUTPUT_CSV_EXON2"'"
        }
    }
}
' "$OUTPUT_FASTA"

NUM_EXON1=$(tail -n +2 "$OUTPUT_CSV_EXON1" | wc -l)
NUM_EXON2=$(tail -n +2 "$OUTPUT_CSV_EXON2" | wc -l)
echo "Created exon 1 CSV with $NUM_EXON1 sequences."
echo "Created exon 2 CSV with $NUM_EXON2 sequences."

# --- STEP 7: Clean up temporary files ---
echo "--- Cleaning up temporary files ---"
rm -r "$TMP_DIR"

echo ""
echo "=========================================="
echo "--- Sequence Extraction Complete! ---"
echo "=========================================="
echo "Exon 1 CSV file saved in: '$OUTPUT_CSV_EXON1'"
echo "Exon 2 CSV file saved in: '$OUTPUT_CSV_EXON2'"
echo "Total exon 1 sequences: $NUM_EXON1"
echo "Total exon 2 sequences: $NUM_EXON2"
echo "Expected approximately: $NUM_GENES sequences per file"