#!/bin/bash

# -----------------------------------------------------------------------------
# fasta_filter.sh
#
# Description:
# A script to filter a FASTA file, keeping only records that have a specific
# keyword (e.g., "lncrna") in their header line.
#
# -----------------------------------------------------------------------------

# --- 1. Check for Correct Usage ---
# The script requires exactly two arguments: the input file and the output file.
# If they are not provided, print a usage message and exit.
if [ "$#" -ne 3 ]; then
    echo "Filters a FASTA file based on a keyword in the header."
    echo "Usage: $0 <input.fasta> <output.fasta> <keyword>"
    exit 1
fi

# --- 2. Assign Arguments to Variables ---
# This makes the script more readable.
INPUT_FILE="$1"
OUTPUT_FILE="$2"
KEYWORD="$3"

# --- 3. Check if Input File Exists ---
# It's good practice to ensure the input file is actually there before running.
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found."
    exit 1
fi

echo "Filtering '$INPUT_FILE' for records with '$KEYWORD' in the header..."

# --- 4. Core Logic using awk ---
# The result is that once a matching header is found, it and all its
# following sequence lines are printed, until a new header is found
# that *doesn't* match, which turns printing off again.

awk '/^>/ { p = /'"$KEYWORD"'/ ? 1 : 0 } p' "$INPUT_FILE" > "$OUTPUT_FILE"


echo "Filtering complete. Output saved to '$OUTPUT_FILE'."
exit 0
