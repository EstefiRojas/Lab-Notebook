#!/bin/bash
# This script retrieves exon 1 and 2 coords from ENSEMBL API using ENSEMBL gene ids.
#
# Dependencies: 
# - jq:
#       https://jqlang.github.io/jq/
#
# Script Name: fetch_ensembl_data.sh
#
# Author: Estefania Rojas
#
# Input: 1. Csv file: containing the sequeces dataset to be analysed. 
#           Ex. data/functional-lncrna-dataset.csv
#        2. Name of the output csv file: coords features for each 
#           sequence in input file 1.
#           Ex. short-ncrna-ensembl-feature
#
###############################################################################
# Bash strict mode. -e report all errors, -u exit on error, -x print every command executed.
set -uex

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

# Store input variables
REGIONS_FILE=$1
OUTPUT_NAME=$2

echo "Processing file $REGIONS_FILE"

# Define path to store output
output_path=../data/datasets/ensembl_feature

# Create output directory if needed
mkdir -p "$output_path"

echo "Output file: $output_path/$OUTPUT_NAME.csv"

# Add header to output file
echo "rnacentral_id,gene_name,chromosome,start,end" > "$output_path"/"$OUTPUT_NAME".csv

# Loop through all ensembl transcript ids
while read -r line; do
	echo "Processing $line"
	# Split line into variables
	IFS='|' read -r t_id1 temp prob2 <<< "$line"
	IFS=' ' read -r t_id2 prob1 <<< "$temp"

	# Remove decimal point for comparison (multiply by 100)
	prob1_int=$(echo $prob1 | awk '{print $1 * 100}')
	prob2_int=$(echo $prob2 | awk '{print $1 * 100}')

	# Compare and select ID
	if (( $(echo "$prob1_int > $prob2_int" | bc -l) )); then
    	selected_id=$t_id1
	else
	    selected_id=$t_id2
	fi

    # Query ENSEMBL API, retry 2 times, and ask for JSON format, parse gene_name, chromosome, start, and end.
    ensembl_response=$(curl -s --retry 2 "https://rest.ensembl.org/lookup/id/$(echo "${selected_id%.*}")?content-type=application/json")

    # Parse JSON response with jq, get rna_type and gene_name fields if present
    fields=$(echo "$ensembl_response" | jq -r '[.display_name // "NA", .seq_region_name // "NA", .start // "NA", .end // "NA"] | @csv')

    # Query RNACentral API, parce RNACentral id of each transcript
    rnacentral_response1=$(curl -s --retry 2 "https://rnacentral.org/api/v1/rna/?external_id=$(echo "${t_id1%.*}")&format=json")
    rnacentral_response2=$(curl -s --retry 2 "https://rnacentral.org/api/v1/rna/?external_id=$(echo "${t_id2%.*}")&format=json")

    rnacentral_id1=$(echo "$rnacentral_response1" | jq -r '[.results[0].rnacentral_id // "NA"] | @csv')
    rnacentral_id2=$(echo "$rnacentral_response2" | jq -r '[.results[0].rnacentral_id // "NA"] | @csv')

    rnacentral_id=$(echo "${rnacentral_id1}|${rnacentral_id2}" | tr -d '"')

    # Send to output file
    echo "${rnacentral_id},${fields}" >> "$output_path/$OUTPUT_NAME.csv"

done < <(awk -F, 'NR > 1 {print $2,$4}' "$REGIONS_FILE")

