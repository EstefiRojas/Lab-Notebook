#!/bin/bash
#
# Dependencies: 
# - bedtools:
#       https://bedtools.readthedocs.io/en/latest/content/installation.html
#
# Script Name: histone_processing.sh
#
# Author: Estefania Rojas
#
# Description: This script calculates the weighted sum of enrichment for sequences of functional long and short-ncRNA 
# databases and protein-coding-RNA databases for a given histone. The script computes the intersection between the sequences
# and histone marks reported in ENCODE using bedtools intersect.
#
# Input: 1. Csv file containing the sequeces dataset to be analysed. Ex. data/functional-lncrna-exon1-dataset.csv
#        2. Name of the histone to analyse. Ex. H3K36me3
#        3. Name of the output csv file to save weighted sum of enrichment feature for each sequence in input file 1.
#           Ex. lncrna-exon1-histone-feature
#
############################################################################################################################

# Bash strict mode. -e exit on error, -u error on unset variables, -x print every command executed.
set -uex

# Get the date of execution. Format example: Tuesday January 17, 2024
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}
##############################################################

# Check for correct usage
if [ $# -ne 3 ]; then
    echo "Usage: $0 regions_to_extract.csv histone_name outname"
    exit 1
fi

# Store input files with descriptive names
REGIONS_FILE=$1
HISTONE_NAME=$2
OUTPUT_NAME=$3

# Temporary files
sorted_regions_file_bed=$(mktemp)
augmented_regions_bed=$(mktemp)

# Variables
data_path=../data/datasets/histone_feature/"$HISTONE_NAME"
marks_path=../data/datasets/epigenetic_data/histone_marks/"$HISTONE_NAME"

# Create a directory if not present already
if [ ! -d "$data_path" ]; then
    mkdir -p "$data_path"
fi


# Remove header from regions file and sort regions by chromosome, start, and end descending
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"


# Heavy processing ahead. Check for previously created files and skip if present.
# Unzip histone marks and join them in one .bed file for sorting
if [ ! -f "$data_path"/"$HISTONE_NAME"_sorted_histone_marks.bed ]; then
    for file in $(find "$marks_path" -type f -print | grep .bed.gz); do
        gunzip "$file"
        #filename=$(basename "$file" .bed)  # Extract filename without .bed extension
        awk -v OFS="\t" '{print $0}' "${file%.gz}" >> "$marks_path"/temp_histone_marks.bed
    done

    # Keep columns Chromosome, Start, End, Enrichment Value, p, and q-value. Sort by chromosome, start, and end descending
    awk -F '\t' '{
             print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9
         }' OFS=, "$marks_path"/temp_histone_marks.bed | \
    sort -k1,1 -k2,2n -k3,3n > "$data_path"/"$HISTONE_NAME"_sorted_histone_marks.bed
fi


# Perform interval intersection using BEDTools between histone marks positions and regions file
bedtools intersect -a "$data_path"/"$HISTONE_NAME"_sorted_histone_marks.bed -b "$sorted_regions_file_bed" -wo > overlaps_"$HISTONE_NAME"_"$OUTPUT_NAME".bed

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
' overlaps_"$HISTONE_NAME"_"$OUTPUT_NAME".bed | sort -k1,1 -k2,2n -k3,3n > "$OUTPUT_NAME".bed


awk -F '\t' 'NR==FNR { seen[$1,$2,$3]=$4; next } 
             { if (($1,$2,$3) in seen) print $0"\t"seen[$1,$2,$3]; else print $0"\t0" }' \
             "$OUTPUT_NAME".bed "$sorted_regions_file_bed" > "$augmented_regions_bed"

# Add header row to output file
echo Chromosome','Start','End','ID','Functional','"$HISTONE_NAME"_AvgSignal > "$data_path"/"$HISTONE_NAME"_"$OUTPUT_NAME".csv

# Convert BED output to desired CSV format
sed 's/\t/,/g' "$augmented_regions_bed" >> "$data_path"/"$HISTONE_NAME"_"$OUTPUT_NAME".csv


# Remove the temporary files
rm "$sorted_regions_file_bed"
rm overlaps_"$HISTONE_NAME"_"$OUTPUT_NAME".bed
rm "$augmented_regions_bed"
rm "$OUTPUT_NAME".bed
