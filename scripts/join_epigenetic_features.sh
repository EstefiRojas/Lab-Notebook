#!/bin/bash
#
# Script Name: join_epigenetic_features.sh
#
# Author: Estefania Rojas
#
# Description: This script joins weighted sum of enrichment features from histone marks of interest for further analysis
# in R.
#
############################################################################################################################

# Bash strict mode. -e report all errors, -u exit on error, -x print every command executed.
set -uex

# Get the date of execution. Format example: Tuesday January 17, 2023
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}
##############################################################


# Check for correct usage
if [ $# -ne 3 ]; then
    echo "Usage: $0 initial_matrix.csv list_of_files.txt output_matrix_name"
    exit 1
fi

# Store input files with descriptive names
first_file=$1
list_of_files=$2
output=$3

# Extract common columns
cut -d $'\t' -f 5-6 "$first_file" > temp2_file.csv

for file in $(cat "$list_of_files"); do
    cut -d $'\t' -f6 $file > temp_file.csv
    paste -d $'\t' temp2_file.csv temp_file.csv > temp3_file.csv
    mv temp3_file.csv temp2_file.csv     # Update temp2 file
done


echo Functionality$'\t'H2KAFZ$'\t'H2AK5ac$'\t'H2AK9ac$'\t'H2BK120ac$'\t'H2BK12ac$'\t'H2BK15ac$'\t' \
     H2BK20ac$'\t'H2BK5ac$'\t'H3KF3A$'\t'H3K14ac$'\t'H3K18ac$'\t'H3K20me1$'\t'H3K23ac$'\t'H3K23me2$'\t' \
     H3K27ac$'\t'H3K27me3$'\t'H3K36me3$'\t'H3K4ac$'\t'H3K4me1$'\t'H3K4me3$'\t'H3K56ac$'\t'H3K79me1$'\t' \
     H3K79me2$'\t'H3K9ac$'\t'H3K9me1$'\t'H3K9me2$'\t'H3K9me3$'\t'H3T11ph$'\t'H4K12ac$'\t'H4K4me2$'\t' \
     H4K5ac$'\t'H4K8ac$'\t'H4K91ac > "$output"

grep -v "Functional" temp2_file.csv | awk -F'\t' '{print $0 }' >> "$output"



# Remove temps
rm -rf temp_file.csv
rm -rf temp2_file.csv
rm -rf temp3_file.csv
