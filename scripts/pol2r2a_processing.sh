#!/bin/bash

set -uex

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
output_path=../data/pol2r2a_feature/
beds_path=../data/POL2R2A/

# Create a directory if not present already
if [ ! -d "$output_path"/sorted_beds/ ]; then
    mkdir -p "$output_path"/sorted_beds/
fi

# Remove header from regions file and sort regions
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"

# Find all bed.gz files in the beds_path directory
for file in $(find "$beds_path" -type f -print | grep .bed.gz); do
    
    filename=$(basename "$file" .gz)
    # Only sort the file if it has not been sorted previously
    if [ ! -f "$output_path"/sorted_beds/"${filename%.bed}"_sorted.bed ]; then
        # Unzip each archive
        gunzip "$file"
        #Get columns chrom, start, end, and enrichment value, and sort descending by Chromosome name, start, and end position.
        awk -F '\t' '{ print $1"\t"$2"\t"$3"\t"$7 }' "${file%.gz}" | \
        sort -k1,1 -k2,2n -k3,3n > "$output_path"/sorted_beds/"${filename%.bed}"_sorted.bed
    fi

    # Perform interval intersection using BEDTools
    bedtools intersect -a "$output_path"/sorted_beds/"${filename%.bed}"_sorted.bed -b "$sorted_regions_file_bed" -wo >> overlaps_"$OUTPUT_NAME".bed

    if [ ! -f "$file" ]; then
        # Compress bed file again to regain disk space
        gzip "${file%.gz}"
    fi
done

# Group hits by chromosome, start, and end of reads to obtain the weighted sum of enrichment
awk -F '\t' '
    BEGIN { OFS="\t" }
    {
        key = $1 FS $6 FS $7;
        sum[key] += $4 * ($10 / ($7 - $6));
    }
    END {
        for (k in sum) {
            print k, sum[k];
        }
    }
' overlaps_"$OUTPUT_NAME".bed | sort -k1,1 -k2,2n -k3,3n > "$OUTPUT_NAME".bed

# Keep wheighted sum of enrichment for all hits, put 0 to the rest
awk -F '\t' 'NR==FNR { seen[$1,$2,$3]=$4; next } 
             { if (($1,$2,$3) in seen) print $0"\t"seen[$1,$2,$3]; else print $0"\t0" }' \
             "$OUTPUT_NAME".bed "$sorted_regions_file_bed" > "$augmented_regions_bed"


# Add header row
echo Chromosome'\t'Start'\t'End'\t'ID'\t'Functional'\t'AvgSignal > "$output_path"/"$OUTPUT_NAME".csv

# Convert BEDTools output to desired CSV format
awk -F '\t' '{print $0}' OFS=, "$augmented_regions_bed" >> "$output_path"/"$OUTPUT_NAME".csv


# Remove the temporary files
rm "$sorted_regions_file_bed"
rm "$augmented_regions_bed"
rm "$OUTPUT_NAME".bed

