#!/bin/bash

# File names (You can change these or pass them as arguments)
INPUT_FILE=$1      # File 1 (The large file to update)
MAPPING_FILE=$2    # File 2 (The list of new names)
OUTPUT_FILE="updated_gRNA_file.txt"

# Check if files exist
if [[ ! -f "$INPUT_FILE" ]] || [[ ! -f "$MAPPING_FILE" ]]; then
    echo "Error: Input files not found. Please ensure $INPUT_FILE and $MAPPING_FILE exist."
    exit 1
fi

echo "Processing files..."

# Use awk to process the files
# We use a temporary file strategy to ensure robust joining
awk '
BEGIN {
    # Set Output Field Separator to Tab (standard for this data type)
    OFS = "\t"
}

# 1. PROCESSING THE MAPPING FILE (NR == FNR)
NR == FNR {
    # We create a "clean" ID by stripping version numbers (everything after a dot)
    # This handles cases where one file has "ENSG.1" and the other has "ENSG"
    clean_id = $1
    sub(/\..*/, "", clean_id)
    
    # Store the name in an array: map[ENSG_ID] = Gene_Name
    map[clean_id] = $2
    next
}

# 2. PROCESSING THE INPUT FILE (NR != FNR)
{
    # Print the header line exactly as is
    if (FNR == 1) {
        print
        next
    }

    # Isolate the ID from Column 5 (lncRNA_ENSG_ID)
    current_id = $5
    
    # Clean the ID (remove version number) for matching purposes
    clean_id = current_id
    sub(/\..*/, "", clean_id)

    # Check if this ID exists in our mapping array
    if (clean_id in map) {
        # Update Column 10 (Gene_Name) with the value from the map
        $10 = map[clean_id]
    }

    # Print the (potentially modified) line
    print
}
' "$MAPPING_FILE" "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Done! Updated data saved to $OUTPUT_FILE"