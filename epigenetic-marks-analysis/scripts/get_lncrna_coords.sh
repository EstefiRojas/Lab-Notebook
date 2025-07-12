#!/bin/bash
# This script joins exon 1 and 2 coords for lncRNAs from ENSEMBL API using 
# ENSEMBL gene and transcript ids.
#
# Dependencies: 
# - jq:
#       https://jqlang.github.io/jq/
#
# Script Name: join_ensembl_exon_coords.sh
#
# Author: Estefania Rojas
#
# Input: 1. Csv file: containing the dataset of sequences to be analysed. This 
# 			file should have the gene id in the first column; and exon 1 and 2 
#			transcript ids in the second column
# 			in the following format: ENST00000667305.1|ENST00000667305.1".
#           Ex. ../data/functional-lncrna-dataset.csv
#		 2. Coordinates file
#        3. Name of the output csv file: coords features for each sequence in 
# 			input file 1.
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
if [ $# -ne 3 ]; then
    echo "Usage: $0 regions_to_extract.csv coords_file.csv outname.csv"
    exit 1
fi

# Store input variables
REGIONS_FILE=$1
COORDS_FILE=$2
OUTPUT_NAME=$3
echo "Processing file $REGIONS_FILE"

# Define path to store output
output_path=../data/datasets/ensembl_coords_feature
# Define sorted file paths
sorted_regions_file="$output_path/${OUTPUT_NAME}_sorted_regions.tmp"
sorted_coords_file="$output_path/${OUTPUT_NAME}_sorted_coords.tmp"

# Create output directory if needed
mkdir -p "$output_path"

# Only run this part if file do not already exists.
if [ ! -e "$output_path/$OUTPUT_NAME"_coordsId.csv ]; then
	echo "CoordsID" > "$output_path/$OUTPUT_NAME"_coordsId.csv
	# Process the regions file first, skipping header
	awk -F, 'NR > 1 { print $1,$2 }' "$REGIONS_FILE" | while read -r gene_id trans_ids; do
	    # Split line into variables
	    IFS='|' read -r t_id1 t_id2 <<< "$trans_ids"
	    echo "${gene_id}/${t_id1},${gene_id}"
	done >> "$output_path/$OUTPUT_NAME"_coordsId.csv
fi

# Same for paste command
if [ ! -e "$output_path/$OUTPUT_NAME"_GeneID.csv ]; then
	# Now paste files (both should have same number of lines)
	paste -d',' "$output_path/$OUTPUT_NAME"_coordsId.csv "$REGIONS_FILE" > "$output_path/$OUTPUT_NAME"_GeneID.csv
fi


# Normalize and preprocess the regions file
awk -F, 'NR > 1 { 
    gsub(/\.[0-9]+$/, "", $1); 
    print $1 "," $2
}' "$output_path/$OUTPUT_NAME"_GeneID.csv | sort -t',' -k1,1 > "$sorted_regions_file"

# Normalize and preprocess the coords file, keeping only specific columns (e.g., 3rd, 4th, 5th)
awk -F, 'NR > 1 {
    gsub(/\.[0-9]+$/, "", $7); 
    print $7 "," $3 "," $4 "," $5
}' "$COORDS_FILE" | sort -t',' -k1,1 > "$sorted_coords_file"


# join two files on their 7th fields
join -t',' -1 1 -2 1 -a 1 "$sorted_regions_file" "$sorted_coords_file" > "$output_path/$OUTPUT_NAME".csv

# Debugging
if [ ! -s "$output_path/$OUTPUT_NAME".csv ]; then
    echo "Join operation returned an empty file. Check input files for issues."
    wc -l "$sorted_regions_file" "$sorted_coords_file"
    cut -d',' -f1 "$sorted_regions_file" | sort | uniq | wc -l
    cut -d',' -f1 "$sorted_coords_file" | sort | uniq | wc -l
fi

# Clean up temporary files
rm "$sorted_regions_file" "$sorted_coords_file"
