#!/bin/bash

# Check if input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 input_file.csv"
    exit 1
fi

input_file=$1

# Create a temporary directory
temp_dir=$(mktemp -d)
trap 'rm -rf "$temp_dir"' EXIT

# Extract header line and write to output
head -1 "$input_file" > "$temp_dir/header.txt"
cat "$temp_dir/header.txt"

# Read the header to get column names
IFS=',' read -ra columns < "$temp_dir/header.txt"
num_columns=${#columns[@]}

# Extract all unique IDs from all columns (excluding header and "pseudo" IDs)
for ((i=1; i<=num_columns; i++)); do
    cut -d',' -f$i "$input_file" | tail -n +2 | grep -v "pseudo" >> "$temp_dir/all_ids.txt"
done

# Sort and get unique IDs (removing empty lines)
sort "$temp_dir/all_ids.txt" | uniq | grep -v "^$" > "$temp_dir/unique_ids.txt"

# For each unique ID, check which columns it appears in
while IFS= read -r id; do
    line=""
    for ((col=1; col<=num_columns; col++)); do
        # Check if this ID exists in this column
        if grep -q "^$id$" <(cut -d',' -f$col "$input_file" | tail -n +2); then
            if [ $col -eq $num_columns ]; then
                line+="$id"
            else
                line+="$id,"
            fi
        else
            if [ $col -eq $num_columns ]; then
                line+=""
            else
                line+=","
            fi
        fi
    done
    echo "$line"
done < "$temp_dir/unique_ids.txt"