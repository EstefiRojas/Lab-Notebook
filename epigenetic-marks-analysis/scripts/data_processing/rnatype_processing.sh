#!/bin/bash
#
# Dependencies: 
# - jq:
#       https://jqlang.github.io/jq/
#
# Script Name: rnatype_processing.sh
#
# Author: Estefania Rojas
#
# Description: This script retrieves the rna type annotation from RNAcentral API
# 			   for each sequence in the dataset.
#
# Input: 1. Csv file containing the sequeces dataset to be analysed. 
#           Ex. data/functional-short-ncrna-dataset.csv
#        2. Name of the output csv file to save mean and max features for each 
#           sequence in input file 1.
#           Ex. short-ncrna-rnatype-feature
#
###############################################################################
# Bash strict mode. -e report all errors, -u exit on error, -x print every command executed.
#set -uex

# Get the date of execution. Format example: Tuesday January 17, 2024
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}
##############################################################

# Check for correct usage
if [ $# -ne 2 ]; then
    echo "Usage: $0 regions_to_extract.csv outname.csv"
    exit 1
fi

# Store input files
REGIONS_FILE=$1
OUTPUT_NAME=$2

echo "Processing file $REGIONS_FILE"

# Define path to store output
output_path=../data/datasets/rnatype_feature

# Create output directory if needed
mkdir -p "$output_path"

echo "Output file: $output_path/$OUTPUT_NAME.csv"

# Print output header
echo "rna_type,gene_name,rnacentral_id,hgnc_id" > "$output_path"/"$OUTPUT_NAME".csv

while read -r line; do
	echo "Processing $line"
    # Query RNAcentral API with ID, retry 2 times, and ask for JSON format
    json_response=$(curl -s --retry 2 "https://rnacentral.org/api/v1/rna/$(echo "${line}" | tr -d '"')?format=json")

    # Parse JSON response with jq, get rna_type and gene_name fields if present, NA otherwise
    fields=$(echo "$json_response" | jq -r 'try ([.rna_type // "NA", .genes[0] // "NA"] | @csv) // "NA,NA"')

    # Query HGNC API with RNAcentral ID, request json format, parse HGNC ID:
    hgnc_id=$(curl -H "Accept: application/json" -s --retry 2 \
              "https://rest.genenames.org/search/rna_central_id/$(echo "${line}" | tr -d '"' | sed 's/_.*//')" \
              | jq -r 'try ([.response.docs[0].hgnc_id // "NA"] | @csv) // "NA"')

    # Send to output file
    echo "${fields},${line},${hgnc_id}" >> "$output_path/$OUTPUT_NAME.csv"


done < <(awk -F, 'NR > 1 {print $7}' "$REGIONS_FILE")

