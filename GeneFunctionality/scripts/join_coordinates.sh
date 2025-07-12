#!/bin/bash
#
###############################################################################
# Get the date of execution. Format example: Tuesday January 17, 2024
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}

# Check for correct usage
if [ $# -ne 2 ]; then
    echo "Usage: $0 file1.csv file2.csv"
    exit 1
fi
##############################################################

INPUT1=$1
INPUT2=$2


# Ensure both files have consistent line endings
tr -d '\r' < ${INPUT1} > first_clean.csv
tr -d '\r' < ${INPUT2} > second_clean.csv

# 1. Extract headers
head -1 first_clean.csv > header1.csv
#head -1 ${INPUT2} | cut -d, -f2- > header2.csv
echo "chromosome,strand,gene_start,gene_end" > header2.csv

# 3. Join the headers
paste -d, header1.csv header2.csv > merged_header.csv

# 4. Remove headers from input files
tail -n +2 first_clean.csv > first_noheader.csv

awk '{print $4","$1","$5","$2","$3}' second_clean.csv > fixed_input2.csv
tail -n +1 fixed_input2.csv > second_noheader.csv

# 1. First, ensure both files are sorted by the join field
sort -t, -k1,1 first_noheader.csv > sorted_first.csv
sort -t, -k1,1 second_noheader.csv > sorted_second.csv


# 5. Join the data
join -t, -1 1 -2 1 sorted_first.csv sorted_second.csv > joined_data.csv

# 6. Combine header and joined data
cat merged_header.csv joined_data.csv > ../data/merged_coordinates_hg38.csv

# 7. Clean up temporary files
rm first_clean.csv second_clean.csv first_noheader.csv second_noheader.csv fixed_input2.csv sorted_first.csv sorted_second.csv header1.csv header2.csv merged_header.csv joined_data.csv
