#!/bin/bash

# This script finds ENSG IDs that are present in the first file but not in the second.
# Usage: ./find_missing_ensg.sh <path_to_file1> <path_to_file2>

# Check if exactly two arguments are provided.
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <file1.csv> <file2.csv>"
    exit 1
fi

FILE1=$1
FILE2=$2

# Ensure both input files exist before proceeding.
if [ ! -f "$FILE1" ]; then
    echo "Error: File '$FILE1' not found."
    exit 1
fi

if [ ! -f "$FILE2" ]; then
    echo "Error: File '$FILE2' not found."
    exit 1
fi

#echo "The following ENSG IDs are in $FILE1 but not in $FILE2:"

# Use awk for a robust solution.
# The gsub() function is added to remove potential trailing carriage returns (\r)
# from the fields, which can cause comparison issues with files from different OS.
#
# 1. BEGIN { FS = "," }: Sets the field separator to a comma.
# 2. NR==FNR { ... }: This block runs only for the first input file (FILE2).
#    - It reads column 11, cleans it, and stores it as a key in the 'seen' array.
# 3. Second block runs for the second file (FILE1).
#    - It cleans column 5.
#    - It checks if the cleaned ID is NOT in the 'seen' array.
#    - '!printed[$5]++' ensures that each unique ID is printed only once.
awk '
    BEGIN { FS = "," }
    NR==FNR {
        if (FNR > 1) {
            id = $11;
            gsub(/\r/,"", id);
            seen[id] = 1;
        }
        next
    }
    (FNR > 1) {
        id = $5;
        gsub(/\r/,"", id);
        if (!(id in seen) && !printed[id]++) {
            print id;
        }
    }
' "$FILE2" "$FILE1"

