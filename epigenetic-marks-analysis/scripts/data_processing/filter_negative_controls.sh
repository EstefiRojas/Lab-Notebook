#!/bin/bash

set -uex

# Check for correct usage
if [ $# -ne 2 ]; then
    echo "Usage: $0 query outname"
    exit 1
fi

# Store input files with descriptive names
REGIONS_FILE=$1
OUTPUT_NAME=$2

sorted_regions_file_bed=$(mktemp)

# Remove header from regions file and sort regions
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"

# Return all rows that dont intersect.
bedtools intersect -wa \
     -a "$sorted_regions_file_bed" \
     -b uniprot_sorted.bed gencode_v45_sorted.bed \
     -sorted \
     -v > "$OUTPUT_NAME".bed
#\-b ../data/F5.hg38.enhancers-sorted.bed ../data/hg38_fair+new_CAGE_peaks_phase1and2-promoter-FANTOM5-sorted.bed \


# Add header row
echo ID,Functional,Chromosome,Start,End > "$OUTPUT_NAME".csv

awk -F '\t' '{print $4,$5,$1,$2,$3}' OFS=, "$OUTPUT_NAME".bed >> "$OUTPUT_NAME".csv

rm "$sorted_regions_file_bed"
rm "$OUTPUT_NAME".bed
