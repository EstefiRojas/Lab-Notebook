#!/opt/local/bin/bash

# Script to horizontally combine columns 5-8 from CSV files based on filename patterns.
# It iterates through histone subdirectories within a parent directory.
# For each CSV file pattern (e.g., lncrna-exon1, protein-exon2-NC), it finds all
# matching files (e.g., H2AFZ/..._lncrna-exon1_..., H2AK5ac/..._lncrna-exon1_...).
# It extracts columns 5-8 from each matching file (skipping headers) into temporary files,
# and then uses 'paste' to combine these temporary files side-by-side into an
# output file named combined_PATTERN.csv.
# A header row is generated indicating the source histone mark for each block of columns.
# Compatible with bash versions < 4.0 (avoids associative arrays).

# Exit script immediately if a command exits with a non-zero status,
# if an unset variable is used, or if any command in a pipeline fails.
set -euo pipefail

# --- Configuration ---
OUTPUT_DIR="pasted_pattern_outputs" # Directory to store the output files
FIELD_SEPARATOR=','                # Set the delimiter for your CSV files
# Create a temporary directory for intermediate files
TEMP_DIR=$(mktemp -d -t paste_csv_temp.XXXXXX) || { echo "Failed to create temp dir"; exit 1; }

# --- Cleanup Function ---
# Ensures temporary directory is removed on script exit (normal or error)
cleanup() {
  echo "Cleaning up temporary directory: $TEMP_DIR"
  rm -rf "$TEMP_DIR"
}
# Register the cleanup function to run on script exit
trap cleanup EXIT

# --- Usage Function ---
usage() {
  echo "Usage: $0 <parent_directory>"
  echo "  <parent_directory>: The directory containing histone subdirectories (e.g., H2AFZ, H2AK5ac)."
  echo "  Processes all .csv files found in the histone subdirectories."
  echo "  Groups files by pattern extracted from filename (e.g., lncrna-exon1, protein-exon2-NC)."
  echo "  For each pattern, extracts columns 5-8 from matching files and pastes them horizontally."
  echo "  Output files are saved to:"
  echo "    ${OUTPUT_DIR}/combined_PATTERN.csv"
  echo "  Assumes CSV files have a header row which is skipped."
  exit 1
}

# --- Input Validation ---
if [ "$#" -ne 1 ]; then
  echo "Error: Incorrect number of arguments." >&2
  usage
fi

PARENT_DIR="$1"

if [ ! -d "$PARENT_DIR" ]; then
  echo "Error: Parent directory not found: $PARENT_DIR" >&2
  exit 1
fi

# --- Processing ---
echo "Starting CSV column pasting based on filename patterns..."
echo "Parent directory: $PARENT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Temporary directory: $TEMP_DIR"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# First pass: Find all unique base patterns
echo "Identifying unique patterns..."
unique_patterns=$(find "$PARENT_DIR" -mindepth 2 -maxdepth 2 -type f -name '*.csv' -print0 | while IFS= read -r -d $'\0' csv_file; do
    filename=$(basename "$csv_file")
    # Extract the base pattern from the filename
    tmp="${filename#*_}"
    base_pattern="${tmp%-histone-feature.csv}"
    # Validate pattern extraction
    if [[ "$base_pattern" != "$filename" && "$base_pattern" != "$tmp" ]]; then
        echo "$base_pattern" # Print valid patterns
    fi
done | sort -u) # Sort and get unique patterns

if [ -z "$unique_patterns" ]; then
    echo "Error: No valid filename patterns found. Ensure files match '*_PATTERN-histone-feature.csv' format." >&2
    exit 1
fi

# Second pass: Process each unique pattern
echo "Processing patterns and pasting columns..."
# Use printf to handle potential newlines in unique_patterns correctly if any pattern extraction failed unexpectedly
printf "%s\n" "$unique_patterns" | while IFS= read -r base_pattern; do
    # Skip empty lines potentially introduced if pattern extraction had issues
    [ -z "$base_pattern" ] && continue

    echo "Processing pattern: '$base_pattern'"
    output_file="${OUTPUT_DIR}/combined_${base_pattern}.csv"

    # Find all files matching this specific pattern
    # Use an array to store found files
    files_for_pattern=()
    while IFS= read -r -d $'\0' matching_file; do
        files_for_pattern+=("$matching_file")
    done < <(find "$PARENT_DIR" -mindepth 2 -maxdepth 2 -type f -name "*_${base_pattern}-histone-feature.csv" -print0)

    # Sort files alphabetically to ensure consistent column order
    IFS=$'\n' sorted_files=($(sort <<<"${files_for_pattern[*]}"))
    unset IFS

    if [ ${#sorted_files[@]} -eq 0 ]; then
        echo "Warning: No files found for pattern '$base_pattern'. Skipping." >&2
        continue
    fi

    echo "  -> Found ${#sorted_files[@]} files. Output file: '$output_file'"

    # --- Generate Header ---
    header_cols=()
    for file in "${sorted_files[@]}"; do
        # Extract histone prefix (parent directory name)
        histone_prefix=$(basename "$(dirname "$file")")
        # Add column headers for this file
        header_cols+=("${histone_prefix}_AvgSignal" "${histone_prefix}_MaxSignal" "${histone_prefix}_MaxScaledSignal")
    done
    # Join header columns with the specified separator
    header_line=$(IFS="$FIELD_SEPARATOR"; echo "${header_cols[*]}")
    # Write header to the output file (overwrite if exists)
    echo "$header_line" > "$output_file"
    echo "  -> Header generated."

    # --- Prepare temporary files containing extracted columns ---
    temp_files_for_paste=() # Array to hold names of temporary files
    echo "  -> Extracting columns to temporary files..."
    for file in "${sorted_files[@]}"; do
        # Create a unique temporary file for this input file's columns
        temp_col_file=$(mktemp "${TEMP_DIR}/${base_pattern}_$(basename "$file")_cols.XXXXXX") || { echo "Failed to create temp column file"; exit 1; }
        temp_files_for_paste+=("$temp_col_file") # Add its name to the list for paste

        # awk extracts columns 5-8, skips header (FNR>1), uses correct FS and OFS
        # Output is redirected to the temporary file
        awk -F"$FIELD_SEPARATOR" -v OFS="$FIELD_SEPARATOR" '
            FNR>1 {
                if(NF>=8) {
                    print $6,$7,$8
                } else {
                    # Print empty fields to maintain row count for paste alignment
                    print "","",""
                }
            }' "$file" > "$temp_col_file" || { echo "Awk failed for $file"; exit 1; }
    done
    echo "  -> Column extraction complete."

    # --- Execute 'paste' using temporary files and append to output file ---
    echo "  -> Pasting data from temporary files..."
    # Ensure temp_files_for_paste are expanded correctly
    if paste -d "$FIELD_SEPARATOR" "${temp_files_for_paste[@]}" >> "$output_file"; then
      echo "  -> Data pasted successfully."
    else
      # Note: Paste might still fail if row counts differ significantly between files,
      # even with the empty field printing in awk.
      echo "Warning: paste command failed for pattern '$base_pattern'. Output file might be incomplete." >&2
      # Consider removing the potentially incomplete file: rm -f "$output_file"
    fi

    # Temporary files for this pattern are cleaned up by the main trap EXIT

done

echo "Processing complete. Output files are in '$OUTPUT_DIR'."

# Exit successfully. (Cleanup trap will run automatically)
exit 0
