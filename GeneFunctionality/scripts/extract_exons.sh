#!/bin/bash

# This script takes a list of ENSG gene IDs and fetches data for exon 1 and 2
# of the canonical transcript for each gene.
# It produces three files:
# 1. A FASTA file with the exon sequences.
# 2. A CSV file with details for exon 1.
# 3. A CSV file with details for exon 2.
# Dependencies: curl, jq
# Usage: ./fetch_gene_exons.sh <input_file_with_ids> <output_fasta_file_base>

# --- Argument and Dependency Check ---
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file_with_ids> <output_fasta_file>"
    exit 1
fi

if ! command -v curl &> /dev/null; then
    echo "Error: 'curl' is not installed. Please install it to use this script."
    exit 1
fi

if ! command -v jq &> /dev/null; then
    echo "Error: 'jq' is not installed. Please install it to use this script (e.g., sudo apt-get install jq)."
    exit 1
fi

INPUT_FILE=$1
OUTPUT_FASTA=$2
# Derive CSV filenames from the FASTA output filename
OUTPUT_BASE="${OUTPUT_FASTA%.*}"
OUTPUT_CSV_E1="${OUTPUT_BASE}_exon1.csv"
OUTPUT_CSV_E2="${OUTPUT_BASE}_exon2.csv"

ENSEMBL_API="https://rest.ensembl.org"

# --- Initialize Output Files ---
# Clear the FASTA file to start fresh
> "$OUTPUT_FASTA"
# Create CSV files and write headers.
echo "ID,Functional,Chromosome,Start,End,Sequence,GeneID" > "$OUTPUT_CSV_E1"
echo "ID,Functional,Chromosome,Start,End,Sequence,GeneID" > "$OUTPUT_CSV_E2"

echo "Starting download of exon data..."
echo "FASTA output will be saved to: $OUTPUT_FASTA"
echo "Exon 1 CSV output will be saved to: $OUTPUT_CSV_E1"
echo "Exon 2 CSV output will be saved to: $OUTPUT_CSV_E2"
echo "----------------------------------------"

# --- Main Loop ---
# Read each gene ID, trimming whitespace automatically.
while read -r gene_id || [[ -n "$gene_id" ]]; do
    # Skip empty lines
    if [ -z "$gene_id" ]; then
        continue
    fi

    echo "Processing Gene ID: $gene_id"

    # 1. Look up the gene ID to get its versioned ID and canonical transcript ID.
    gene_lookup_json=$(curl -s -X GET "${ENSEMBL_API}/lookup/id/${gene_id}" -H "Content-Type:application/json")
    
    versioned_gene_id=$(echo "$gene_lookup_json" | jq -r '.id')
    canonical_transcript_id=$(echo "$gene_lookup_json" | jq -r '.canonical_transcript')

    if [ -z "$canonical_transcript_id" ] || [ "$canonical_transcript_id" == "null" ]; then
        echo "  -> Warning: No canonical transcript found for $gene_id. Skipping."
        continue
    fi
    
    ct_id_no_version=${canonical_transcript_id%%.*}
    echo "  -> Found Canonical Transcript: $ct_id_no_version"

    # 2. Look up the canonical transcript to get a list of its exon data.
    # We will capture ID, chromosome, start, and end for each exon.
    exon_data=()
    while IFS= read -r line; do
        exon_data+=("$line")
    done < <(curl -s -X GET "${ENSEMBL_API}/lookup/id/${ct_id_no_version}?expand=exons" -H "Content-Type:application/json" | jq -r '.Exon[] | "\(.id),\(.seq_region_name),\(.start),\(.end)"')

    if [ ${#exon_data[@]} -eq 0 ]; then
        echo "  -> Warning: No exons found for transcript $ct_id_no_version. Skipping."
        continue
    fi

    # 3. Process Exon 1
    IFS=',' read -r exon1_id exon1_chr exon1_start exon1_end <<< "${exon_data[0]}"
    if [ -n "$exon1_id" ]; then
        echo "  -> Fetching Exon 1: $exon1_id"
        sequence_data=$(curl -sf -X GET "${ENSEMBL_API}/sequence/id/${exon1_id}?content-type=text/plain")
        if [ -n "$sequence_data" ]; then
            # Write to FASTA
            echo ">${versioned_gene_id}_${ct_id_no_version}_exon_1" >> "$OUTPUT_FASTA"
            echo "$sequence_data" >> "$OUTPUT_FASTA"
            echo "" >> "$OUTPUT_FASTA"
            # Write to CSV
            csv_id="${ct_id_no_version}_exon_1"
            gene_id_col="${versioned_gene_id}/${canonical_transcript_id}"
            echo "${csv_id},NA,${exon1_chr},${exon1_start},${exon1_end},${sequence_data},${gene_id_col}" >> "$OUTPUT_CSV_E1"
        fi
        sleep 0.2
    fi

    # 4. Process Exon 2, if it exists
    if [ ${#exon_data[@]} -gt 1 ]; then
        IFS=',' read -r exon2_id exon2_chr exon2_start exon2_end <<< "${exon_data[1]}"
        if [ -n "$exon2_id" ]; then
            echo "  -> Fetching Exon 2: $exon2_id"
            sequence_data=$(curl -sf -X GET "${ENSEMBL_API}/sequence/id/${exon2_id}?content-type=text/plain")
            if [ -n "$sequence_data" ]; then
                # Write to FASTA
                echo ">${versioned_gene_id}_${ct_id_no_version}_exon_2" >> "$OUTPUT_FASTA"
                echo "$sequence_data" >> "$OUTPUT_FASTA"
                echo "" >> "$OUTPUT_FASTA"
                # Write to CSV
                csv_id="${ct_id_no_version}_exon_2"
                gene_id_col="${versioned_gene_id}/${canonical_transcript_id}"
                echo "${csv_id},NA,${exon2_chr},${exon2_start},${exon2_end},${sequence_data},${gene_id_col}" >> "$OUTPUT_CSV_E2"
            fi
            sleep 0.2
        fi
    fi

done < "$INPUT_FILE"

echo "----------------------------------------"
echo "Script finished. All files have been created."

