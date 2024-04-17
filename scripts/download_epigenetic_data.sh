#!/bin/bash
#
# Dependencies: 
# - markovProperties.pl:
#
# Script Name: download_epigenetic_data.sh
#
# Author: Estefania Rojas
#
# Description: This script downloads epigenetic data in bed format from ENCODE database by finding all links provided in
# the folder data/epigenetic_data/. If some or all of the files are already present in the corresponding directory,
# the script skips that file.
#
############################################################################################################################

# Bash strict mode. -e report all errors, -u exit on error, -x print every command executed.
set -uex

# Get the date of execution. Format example: Tuesday January 17, 2024
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}
##############################################################

echo "Starting Epigenetic marks download..."

# Find all .txt files and loop through each
for file in $(find ../data/epigenetic_data/ -type f -print | grep .txt); do
    mkdir -p "${file%.txt}"  # Create directory if not present

    # Get all links in the txt file
    for link in $(cat "$file" | grep bed.gz); do
        code=$(basename "$link" .bed.gz) # Extract the code name from the link
        # Only download if not present already
        if [ ! -f "${file%.txt}"/"$code".bed.gz ]; then
            wget -q "$link" -O "${file%.txt}"/"$code".bed.gz
            echo "Download successful: " "$code"
        else
            echo "File already present: $code"
        fi
    done
    
done
echo "Finished downloading Epigenetic marks."

