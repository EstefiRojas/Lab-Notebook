#!/bin/bash
#
# Dependencies: 
# - bedtools:
#       https://bedtools.readthedocs.io/en/latest/content/installation.html
#
# Script Name: methylome_processing.sh
#
# Author: Estefania Rojas
#
# Description: This script calculates the average methylation percentage for sequences of functional long and short-ncRNA 
# databases and protein-coding-RNA databases. The script computes the intersection between the sequences
# and methylation percentage reported in ENCODE, using bedtools intersect.
#
# Input: 1. Csv file containing the sequeces dataset to be analysed. Ex. data/functional-lncrna-exon1-dataset.csv
#        2. Name of the output csv file to save methylation percentage feature for each sequence in input file 1.
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
if [ $# -ne 2 ]; then
    echo "Usage: $0 regions_to_extract.csv outname"
    exit 1
fi

# Store input files with descriptive names
REGIONS_FILE=$1
OUTPUT_NAME=$2

# Temporary files
sorted_regions_file_bed=$(mktemp)
augmented_regions_bed=$(mktemp)

# Variables
data_path=../data/datasets/methylome_feature
beds_path=../data/datasets/epigenetic_data/methylome

# Create a directory if not present already
if [ ! -d "$data_path"/sorted_beds/ ]; then
    mkdir -p "$data_path"/sorted_beds/
fi

# Remove header from regions file and sort regions and convert to bed format
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"


for file in $(find "$beds_path" -type f -print | grep .bed.gz); do
    filename=$(basename "$file" .gz)
    if [ ! -f "$data_path"/sorted_beds/"${filename%.bed}"_sorted.bed ]; then
        gunzip "$file"
        #Get columns chrom, start, end, methylation level, and percentage of methylation
        awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$10"\t"$11}' "${file%.gz}" | \
        sort -k1,1 -k2,2n -k3,3n > "$data_path"/sorted_beds/"${filename%.bed}"_sorted.bed
    fi

    # Perform interval intersection using BEDTools
    bedtools intersect -a "$data_path"/sorted_beds/"${filename%.bed}"_sorted.bed -b "$sorted_regions_file_bed" -wa -wb >> overlaps_"$OUTPUT_NAME".bed

    if [ ! -f "$file" ]; then
        gzip "${file%.gz}" &
    fi
done

# Group hits by chromosome, start, and end of reads to obtain average methilation percentage
awk -F '\t' '
    BEGIN { OFS="\t" }
    {
        key = $1 FS $7 FS $8;
        sum[key]+=$5;
        count[key]++
    }
    END {
        for (k in sum) {
            print k, sum[k]/count[k];
        }
    }
' overlaps_"$OUTPUT_NAME".bed | sort -k1,1 -k2,2n -k3,3n > "$OUTPUT_NAME".bed


awk -F '\t' 'NR==FNR { seen[$1,$2,$3]=$4; next } 
             { if (($1,$2,$3) in seen) print $0"\t"seen[$1,$2,$3]; else print $0"\t0" }' \
             "$OUTPUT_NAME".bed "$sorted_regions_file_bed" > "$augmented_regions_bed"


# Add header row
echo Chromosome','Start','End','ID','Functional','methylome > "$data_path"/"$OUTPUT_NAME".csv

# Convert BEDTools output to desired CSV format
sed 's/\t/,/g' "$augmented_regions_bed" >> "$data_path"/"$OUTPUT_NAME".csv


# Remove the temporary files
rm "$sorted_regions_file_bed"
rm "$augmented_regions_bed"
rm "$OUTPUT_NAME".bed
rm overlaps_"$OUTPUT_NAME".bed
