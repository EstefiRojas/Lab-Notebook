#!/bin/bash
#
# Dependencies: 
# - bedtools:
#       https://bedtools.readthedocs.io/en/latest/content/installation.html
#
# Script Name: chrm_acc_processing.sh
#
# Author: Estefania Rojas
#
# Description: This script calculates the weighted sum of enrichment for sequences of functional long and short-ncRNA 
# databases and protein-coding-RNA databases for chromatin accessibility. The script computes the intersection between the sequences
# and chromatin accessibility marks reported in ENCODE using bedtools intersect.
#
# Input: 1. Csv file containing the sequeces dataset to be analysed. Ex. data/functional-lncrna-exon1-dataset.csv
#        2. Name of the output csv file to save weighted sum of enrichment feature for each sequence in input file 1.
#           Ex. lncrna-exon1-histone-feature
#
############################################################################################################################

# Bash strict mode. -e report all errors, -u exit on error, -x print every command executed.
set -uex

# Get the date of execution. Format example: Tuesday January 17, 2024
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}
##############################################################


# Check for correct usage
if [ $# -ne 2 ]; then
    echo "Usage: $0 regions_to_extract.csv outname"
    exit 1
fi

# Store input files with descriptive names
REGIONS_FILE=$1
#EXPERIMENT_TYPE=$2
#BIOSAMPLE_CLASS=$3
OUTPUT_NAME=$2


# Temporary files
#temp_chrm_acc_marks_bed=$(mktemp)  # Temporary file for BED conversion
#sorted_chrm_acc_marks_bed=$(mktemp)
sorted_regions_file_bed=$(mktemp)
#overlaps_bed=$(mktemp)
augmented_regions_bed=$(mktemp)

# Variables
data_path=../data/chrm_acc_feature/by_biosample_class
marks_path=../data/epigenetic_data/chromatin_accessibility/

# Create a directory if not present already
if [ ! -d "$data_path" ]; then
    mkdir -p "$data_path"
fi


# Remove header from regions file and sort regions
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"


# Heavy processing ahead. Check for previously created files and skip if present.
if [ ! -f "$data_path"/chrm_acc_sorted_marks.bed ]; then
    for file in $(find "$marks_path" -type f -print | grep .bed.gz); do
        gunzip "$file"
        #filename=$(basename "$file" .bed)  # Extract filename without .bed extension
        awk -v OFS="\t" '{print $0}' "${file%.gz}" >> "$marks_path"/temp_chrm_acc_marks.bed
    done

    awk -F '\t' '{
             print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9
         }' OFS=, "$marks_path"/temp_chrm_acc_marks.bed | \
    sort -k1,1 -k2,2n -k3,3n > "$data_path"/chrm_acc_sorted_marks.bed
fi


# Perform interval intersection using BEDTools
bedtools intersect -a "$data_path"/chrm_acc_sorted_marks.bed -b "$sorted_regions_file_bed" -wo > overlaps_chrm_acc_"$OUTPUT_NAME".bed

# AvgSignal = sum(wi * ei)
# Calculate the weighted sum of enrichment, grouped by Chromosome, StartExon, EndExon
awk -F '\t' '
    BEGIN { OFS="\t" }
    {
        key = $1 FS $8 FS $9;
        sum[key] += $4 * ($12 / ($9 - $8));
    }
    END {
        for (k in sum) {
            print k, sum[k];
        }
    }
' overlaps_chrm_acc_"$OUTPUT_NAME".bed | sort -k1,1 -k2,2n -k3,3n > "$OUTPUT_NAME".bed


awk -F '\t' 'NR==FNR { seen[$1,$2,$3]=$4; next } 
             { if (($1,$2,$3) in seen) print $0"\t"seen[$1,$2,$3]; else print $0"\t0" }' \
             "$OUTPUT_NAME".bed "$sorted_regions_file_bed" > "$augmented_regions_bed"

# Add header row
echo Chromosome'\t'Start'\t'End'\t'ID'\t'Functional'\t'AvgSignal > "$data_path"/chrm_acc_"$OUTPUT_NAME".csv

# Convert BEDTools output to desired CSV format
awk -F '\t' '{print $0}' OFS=, "$augmented_regions_bed" >> "$data_path"/chrm_acc_"$OUTPUT_NAME".csv


# Remove the temporary files
# rm "$temp_chrm_acc_marks_bed"
# rm "$sorted_chrm_acc_marks_bed"
rm "$sorted_regions_file_bed"
#rm "$overlaps_bed"
rm "$augmented_regions_bed"
rm "$OUTPUT_NAME".bed
