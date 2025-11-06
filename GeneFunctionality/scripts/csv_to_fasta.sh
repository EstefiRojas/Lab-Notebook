#!/bin/bash

# Bash strict mode. -e exit on error, -u error on unset variables, -o pipefail ensures pipeline errors are caught
set -ueo pipefail

# Optionally uncomment the line below if you want to see all commands executed
# set -x

# --- Configuration ---
# Define the separator used to join multiple name columns in the FASTA header
NAME_SEPARATOR="_"

# --- Argument Parsing and Validation ---

# Check for correct number of arguments
if [ "$#" -ne 4 ]; then
    # Print usage instructions to standard error (stderr)
    echo "Usage: $0 <input.csv> <name_column_spec> <sequence_column_number> <suffix_column_number>" >&2
    echo "  <name_column_spec> can be:" >&2
    echo "    - A single column number (e.g., 7)" >&2
    echo "    - A range of columns (e.g., 1-3)" >&2
    echo "    - A comma-separated list of columns (e.g., 1,3,2)" >&2
    echo "  <suffix_column_number> is the column number to use as suffix" >&2
    echo "  Example (range): $0 data/regions.csv 1-2 6 3 > output.fasta" >&2
    echo "  Example (list):  $0 data/regions.csv 1,3 6 3 > output.fasta" >&2
    exit 1
fi

# Assign arguments to variables
INPUT_FILE=$1
NAME_COLUMN_SPEC=$2 # e.g., "7" or "1-2" or "1,3"
SEQUENCE_COLUMN=$3
SUFFIX_COLUMN=$4  # Now this is a column number instead of a fixed string

# Validate that SUFFIX_COLUMN is a number
if ! [[ "$SUFFIX_COLUMN" =~ ^[0-9]+$ ]]; then
    echo "Error: Suffix column must be a number: $SUFFIX_COLUMN" >&2
    exit 1
fi

# --- Process Name Column Specification ---

NAME_COLS_LIST="" # Will hold space-separated list of column numbers

if [[ "$NAME_COLUMN_SPEC" == *-* ]]; then
    # Handle range (e.g., 1-2)
    START_COL=$(echo "$NAME_COLUMN_SPEC" | cut -d'-' -f1)
    END_COL=$(echo "$NAME_COLUMN_SPEC" | cut -d'-' -f2)
    # Validate that START_COL and END_COL are numbers
    if ! [[ "$START_COL" =~ ^[0-9]+$ ]] || ! [[ "$END_COL" =~ ^[0-9]+$ ]]; then
        echo "Error: Invalid range format in name column specification: $NAME_COLUMN_SPEC" >&2
        exit 1
    fi
    # Generate space-separated list using seq
    NAME_COLS_LIST=$(seq -s' ' "$START_COL" "$END_COL")
elif [[ "$NAME_COLUMN_SPEC" == *,* ]]; then
    # Handle comma-separated list (e.g., 1,3)
    # Validate that all parts are numbers
    if ! echo "$NAME_COLUMN_SPEC" | grep -E '^[0-9]+(,[0-9]+)*$' > /dev/null; then
         echo "Error: Invalid list format in name column specification: $NAME_COLUMN_SPEC" >&2
         exit 1
    fi
    # Replace commas with spaces
    NAME_COLS_LIST=$(echo "$NAME_COLUMN_SPEC" | tr ',' ' ')
elif [[ "$NAME_COLUMN_SPEC" =~ ^[0-9]+$ ]]; then
    # Handle single column number
    NAME_COLS_LIST="$NAME_COLUMN_SPEC"
else
    echo "Error: Invalid format for name column specification: $NAME_COLUMN_SPEC" >&2
    exit 1
fi

# Optional: Print processing info to stderr
# echo "Input File: ${INPUT_FILE}" >&2
# echo "Name Column Spec: ${NAME_COLUMN_SPEC}" >&2
# echo "Resolved Name Columns: ${NAME_COLS_LIST}" >&2
# echo "Sequence Column: ${SEQUENCE_COLUMN}" >&2
# echo "Suffix Column: ${SUFFIX_COLUMN}" >&2
# echo "Name Separator: ${NAME_SEPARATOR}" >&2
# echo "--------------------------------------------------" >&2

# --- AWK Processing ---

# Process the CSV file using awk to generate FASTA output to standard output (stdout)
# -F',' sets the field separator to a comma
# 'NR > 1' skips the header row
# -v passes shell variables into awk variables
#   NAME_COLS_STR: Space-separated list of columns for the name
#   SEQ_COL: Column number for the sequence
#   SUFFIX_COL: Column number for the suffix
#   SEP: Separator to use between name parts
awk -F',' \
    -v NAME_COLS_STR="$NAME_COLS_LIST" \
    -v SEQ_COL="$SEQUENCE_COLUMN" \
    -v SUFFIX_COL="$SUFFIX_COLUMN" \
    -v SEP="$NAME_SEPARATOR" \
'
NR > 1 {
    # Split the space-separated column list into an array
    split(NAME_COLS_STR, name_cols_array, " ");

    # Construct the name from the specified columns
    constructed_name = "";
    num_name_cols = length(name_cols_array);
    for (i = 1; i <= num_name_cols; i++) {
        col_num = name_cols_array[i]; # Get the actual column number
        # Basic check if column number is valid (greater than 0)
        if (col_num > 0 && col_num <= NF) { # NF is number of fields in current record
             constructed_name = constructed_name $col_num; # Append field content
             # Add separator if not the last part
             if (i < num_name_cols) {
                 constructed_name = constructed_name SEP;
             }
        } else {
            # Handle invalid column number - print error to stderr for this record
             print "Warning: Invalid column number " col_num " for record " NR > "/dev/stderr";
        }
    }

    # Get the suffix from the specified column
    suffix_value = "";
    if (SUFFIX_COL > 0 && SUFFIX_COL <= NF) {
        suffix_value = $SUFFIX_COL;
    } else {
        print "Warning: Invalid suffix column number " SUFFIX_COL " for record " NR > "/dev/stderr";
    }

    # Print the FASTA formatted output for this record
    # Check if constructed_name is not empty before printing
    if (constructed_name != "" && suffix_value != "") {
        print ">" constructed_name "_" suffix_value "\n" $SEQ_COL;
    } else if (constructed_name != "" && suffix_value == "") {
        # Print without suffix if suffix is empty
        print ">" constructed_name "\n" $SEQ_COL;
    } else {
        print "Warning: Could not construct name for record " NR > "/dev/stderr";
    }
}
' \
"$INPUT_FILE"

# Optional completion message to stderr
# echo "Processing complete. Output sent to standard output." >&2