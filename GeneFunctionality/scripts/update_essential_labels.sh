#!/usr/bin/env bash

# Usage: ./update_group_column.sh input.csv output.csv
# Rules for the "group" column (here set as column 8 via COL):
# - "Common essential" -> "common"
# - "Core essential"   -> "core"
# - empty              -> "Non-essential"
# - anything else      -> "specific"

if [ $# -ne 2 ]; then
    echo "Usage: $0 <input.csv> <output.csv>"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"

# Set this to the correct column index (1-based).
# From your example header, 'group' is the 8th column.
COL=8

awk -v COL="$COL" -F',' -v OFS=',' '
NR == 1 {
    # Header: print as-is
    print
    next
}
{
    # Normalize the target column according to the rules
    if ($COL == "Common essential") {
        $COL = "common"
    } else if ($COL == "Core essential") {
        $COL = "core"
    } else if ($COL == "" ) {
        $COL = "Non-essential"
    } else {
        $COL = "specific"
    }

    print
}
' "$INPUT" > "$OUTPUT"
