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

# 1. First, ensure both files are sorted by the join field
sort -t, -k1,1 ${INPUT1} > sorted_first.csv
sort -t, -k1,1 ${INPUT2} > sorted_second.csv

# Ensure both files have consistent line endings
tr -d '\r' < sorted_first.csv > first_file_clean.csv
tr -d '\r' < sorted_second.csv > second_file_clean.csv

# 2. Extract headers
head -1 ${INPUT1} > header1.csv
head -1 ${INPUT2} | cut -d, -f2- > header2.csv

# 3. Join the headers
paste -d, header1.csv header2.csv > merged_header.csv

# 4. Remove headers from sorted files
tail -n +2 first_file_clean.csv > sorted_first_noheader.csv
tail -n +2 second_file_clean.csv > sorted_second_noheader.csv

# 5. Join the data
join -t, -1 1 -2 1 sorted_first_noheader.csv sorted_second_noheader.csv > joined_data.csv

# 6. Combine header and joined data
cat merged_header.csv joined_data.csv > merged_ids_file.csv

# 7. Clean up temporary files
rm first_file_clean.csv second_file_clean.csv sorted_first.csv sorted_second.csv header1.csv header2.csv merged_header.csv sorted_first_noheader.csv sorted_second_noheader.csv joined_data.csv
