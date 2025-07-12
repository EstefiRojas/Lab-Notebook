#!/opt/local/bin/bash
############################################################################################################################
# Dependencies:
# - bedtools: Required for intersecting genomic intervals (specifically `bedtools intersect`).
#   Installation: https://bedtools.readthedocs.io/en/latest/content/installation.html
#
# Script Name: histone_processing.sh (Adapted Version)
#
# Author: Estefania Rojas
#
# Description: This script calculates enrichment statistics (average signal, max signal,
# max scaled signal) from histone mark narrowPeak BED files (e.g., downloaded from ENCODE)
# for genomic regions defined in an input CSV file.
# It intersects the filtered marks with the input regions, and calculates
# summary statistics for each region based on overlapping marks.
# Uses temporary files for intermediate steps and includes cleanup.
#
# Input: 1. CSV file containing the regions dataset to be analysed.
#           Format expected: Header line, then columns ID,Functional,Chromosome,Start,End,...
#           (Script specifically uses columns 3, 4, 5 for chr, start, end and 1, 2 for ID, Functional)
#           e.g., my_regions.csv
#        2. Name of the histone mark (used to find input BED files).
#           e.g., H3K4me3 (script will look for *.bed files in ../data/datasets/histone_marks/H3K4me3/)
#        3. Base name for the output CSV file and temporary file prefixes.
#           e.g., experiment1 (will create ../data/histone_feature/H3K4me3/H3K4me3_experiment1.csv)
#
# Output: A CSV file named <HISTONE_NAME>_<OUTPUT_BASE_NAME>.csv in the
#         ../data/datasets/histone_feature/<HISTONE_NAME>/ directory, containing the original
#         region information plus calculated AvgSignal, MaxSignal, and MaxScaledSignal.
#
############################################################################################################################


# --- Configuration & Input ---

# Check for correct number of arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <regions_to_extract.csv> <histone_name> <output_base_name>"
    echo "Example: $0 my_regions.csv 0.05 H3K4me3 experiment1"
    exit 1
fi

# Store input arguments with descriptive names
REGIONS_FILE=$1       # CSV file containing regions (needs header: ignored, ignored, Chromosome, Start, End, ID, Functional)
HISTONE_NAME=$2       # Name of the histone mark (used for directory structure)
OUTPUT_BASE_NAME=$3   # Base name for output files

# Define paths relative to the script's location or a base directory
# Assumes data is in ../data relative to where the script is run
BASE_DATA_PATH="../data/datasets"
MARKS_PATH="$BASE_DATA_PATH/epigenetic_data/histone_marks/$HISTONE_NAME" # Path to individual histone mark BED files
DATA_PATH="$BASE_DATA_PATH/histone_feature/$HISTONE_NAME" # Path for processed data
OUTPUT_CSV_FILE="$DATA_PATH/${HISTONE_NAME}_${OUTPUT_BASE_NAME}.csv"
OUTPUT_TEMP_PREFIX="${HISTONE_NAME}_${OUTPUT_BASE_NAME}" # Prefix for temporary files

# --- Temporary Files ---
# Using process substitution where possible, but defining temp files for clarity
SORTED_REGIONS_BED=$(mktemp "${OUTPUT_TEMP_PREFIX}_regions_sorted.XXXXXX.bed")
COMBINED_FILTERED_SORTED_MARKS_BED="$DATA_PATH/${HISTONE_NAME}_combined_filtered_sorted_marks.bed" # Combined, filtered, sorted marks
OVERLAPS_BED=$(mktemp "${OUTPUT_TEMP_PREFIX}_overlaps.XXXXXX.bed")
CALCULATED_VALUES_BED=$(mktemp "${OUTPUT_TEMP_PREFIX}_calculated.XXXXXX.bed")
AUGMENTED_REGIONS_BED=$(mktemp "${OUTPUT_TEMP_PREFIX}_augmented.XXXXXX.bed")

# --- Helper Function for Cleanup ---
cleanup() {
  echo "Cleaning up temporary files..."
  rm -f "$SORTED_REGIONS_BED"
  rm -f "$OVERLAPS_BED"
  rm -f "$CALCULATED_VALUES_BED"
  rm -f "$AUGMENTED_REGIONS_BED"
  # Optionally remove the combined marks file if it's considered intermediate
  # rm -f "$COMBINED_FILTERED_SORTED_MARKS_BED"
  echo "Cleanup complete."
}

# Trap signals to ensure cleanup runs even if the script exits unexpectedly
trap cleanup EXIT SIGINT SIGTERM

# --- Preprocessing ---
# Create the output data directory if it doesn't exist
mkdir -p "$DATA_PATH"

# --- Unzip, Combine, Filter, and Sort Histone Mark Files ---
# This block first unzips any .bed.gz files, then combines all resulting
# .bed files for the histone mark, filters columns, and sorts them.
# It checks if the final sorted file exists to avoid reprocessing.

if [ ! -f "$COMBINED_FILTERED_SORTED_MARKS_BED" ]; then
    echo "Processing histone mark files..."
    TEMP_COMBINED_UNFILTERED=$(mktemp "${OUTPUT_TEMP_PREFIX}_combined_unfiltered.XXXXXX.bed")

    # --- NEW: Unzip any .bed.gz files first ---
    echo "Unzipping any .bed.gz files found in $MARKS_PATH..."
    # Find all .bed.gz files within the specified path (up to 2 levels deep)
    # and unzip them in place. '-f' forces overwrite if the .bed file already exists.
    # Using \; executes gunzip for each file found.
    find "$MARKS_PATH" -maxdepth 2 -name '*.bed.gz' -exec gunzip -f {} \;
    echo "Unzipping complete."
    # --- END NEW ---

    # Combine all BED files for the specified histone mark
    # This will now include any files that were just unzipped.
    # Using find is more robust than *.bed if there are many files
    echo "Combining all .bed files..."
    find "$MARKS_PATH" -maxdepth 2 -name '*.bed' -exec cat {} + > "$TEMP_COMBINED_UNFILTERED"
    echo "Combining complete."

    echo "Filtering and sorting combined data..."
    awk -F '\t' '
        BEGIN { OFS="\t" }
        {
            # Basic check for expected number of fields before printing
            # This helps avoid errors if some BED files are malformed
            if (NF >= 9) {
                 print $1, $2, $3, $7, $8, $9; # Print selected columns (chrom, start, end, signalValue, pValue, qValue)
            } else {
                 # Optionally print a warning for malformed lines
                 # print "Warning: Skipping malformed line " NR " in combined input." > "/dev/stderr"
            }
        }
    ' "$TEMP_COMBINED_UNFILTERED" | \
    sort -k1,1 -k2,2n -k3,3n > "$COMBINED_FILTERED_SORTED_MARKS_BED"
    echo "Filtering and sorting complete."

    rm "$TEMP_COMBINED_UNFILTERED" # Clean up intermediate combined file
    echo "Histone mark processing complete. Saved to $COMBINED_FILTERED_SORTED_MARKS_BED"
else
    echo "Using existing combined, filtered, and sorted histone mark file: $COMBINED_FILTERED_SORTED_MARKS_BED"
fi

# --- Prepare Regions File ---

# Remove header from the input regions CSV file
# Extract columns: Chromosome (3), Start (4), End (5), ID (1), Functional (2)
# Convert to BED-like format (tab-separated) and sort
echo "Processing regions file..."
awk -F, 'BEGIN { OFS="\t" } NR > 1 {print $3, $4, $5, $1, $2}' "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$SORTED_REGIONS_BED"
echo "Regions file processed and sorted."

# --- Perform Intersection ---

# Use bedtools intersect to find overlaps between histone marks and regions
# -a: Filtered histone marks (chr, start, end, signal, pval, qval)
# -b: Sorted regions (chr, start, end, ID, Functional)
# -wo: Write the original A and B entries plus the number of overlapping bases
echo "Performing intersection with bedtools..."
bedtools intersect -a "$COMBINED_FILTERED_SORTED_MARKS_BED" -b "$SORTED_REGIONS_BED" -wo > "$OVERLAPS_BED"
echo "Intersection complete. Overlaps stored in $OVERLAPS_BED"

# --- Calculate Enrichment Statistics ---

# Process the overlaps file to calculate:
# 1. AvgSignal: Weighted average signal = sum(signal_i * overlap_i) / region_length
# 2. MaxSignal: Maximum signal value observed in any overlapping mark
# 3. MaxScaledSignal: Maximum value of (signal_i * overlap_i / region_length) across overlaps

echo "Calculating enrichment statistics..."
awk -F '\t' '
    BEGIN { OFS="\t" }
    {
        # Columns from -wo output:
        # 1-6: Fields from file A ( histone marks: chr, start, end, signal, pval, qval)
        # 7-11: Fields from file B ( regions: chr, start, end, ID, Functional)
        # 12: Length of overlap

        # Define the key for the region (Chromosome, StartRegion, EndRegion from file B)
        region_chr = $7;
        region_start = $8;
        region_end = $9;
        key = region_chr FS region_start FS region_end;

        # Get values for calculation
        signal = $4;         # Signal value from histone mark (col 4)
        overlap_len = $12;   # Length of the overlap (col 12)
        region_len = region_end - region_start; # Length of the region from file B

        # Avoid division by zero if region length is 0 (should not happen with valid BED)
        if (region_len > 0) {
            # Calculate weighted signal component for this overlap (signal * overlap_fraction)
            weighted_signal_component = signal * (overlap_len / region_len);
            sum_weighted_signal[key] += weighted_signal_component;

            # Calculate scaled signal for this overlap (signal * overlap / region_length)
            # Note: This is the same as weighted_signal_component in this context
            scaled_signal = weighted_signal_component;

            # Update total overlap length for calculating the final average later (optional, can divide sum by region_len)
            # sum_overlap[key] += overlap_len; # Not strictly needed if calculating weighted average as above

            # Update maximum signal if current signal is higher
            if (!(key in max_signal) || signal > max_signal[key]) {
                max_signal[key] = signal;
            }

            # Update maximum scaled signal if current scaled signal is higher
            if (!(key in max_scaled_signal) || scaled_signal > max_scaled_signal[key]) {
                max_scaled_signal[key] = scaled_signal;
            }
        } else {
            # Handle potential zero-length regions if necessary
            print "Warning: Zero length region encountered:", $7, $8, $9 > "/dev/stderr";
        }
    }
    END {
        # Print the results for each region
        # Output format: Chr, StartRegion, EndRegion, AvgSignal, MaxSignal, MaxScaledSignal
        for (k in sum_weighted_signal) {
            # AvgSignal is the sum of weighted components (already normalized by region length)
            avg_sig = sum_weighted_signal[k];

            # Get max values, defaulting to 0 if not found (e.g., if only zero-length regions overlapped)
            max_sig = (k in max_signal) ? max_signal[k] : 0;
            max_scaled_sig = (k in max_scaled_signal) ? max_scaled_signal[k] : 0;

            print k, avg_sig, max_sig, max_scaled_sig;
        }
    }
' "$OVERLAPS_BED" | sort -k1,1 -k2,2n -k3,3n > "$CALCULATED_VALUES_BED"
echo "Enrichment calculations complete. Results stored in $CALCULATED_VALUES_BED"

# --- Merge Calculated Values with Original Regions ---

# Use awk to join the calculated values back to the sorted regions file.
# If a region had no overlaps, assign 0 to the calculated fields.
echo "Merging calculated values with regions..."
awk -F '\t' '
    BEGIN { OFS="\t" }
    # Process the first file (calculated values)
    NR==FNR {
        # Store AvgSignal, MaxSignal, and MaxScaledSignal using the region key (chr, start, end)
        key = $1 FS $2 FS $3;
        avg_signal[key] = $4;
        max_signal[key] = $5;
        max_scaled_signal[key] = $6;
        next
    }
    # Process the second file (sorted regions: chr, start, end, ID, Functional)
    {
        # Define the key for the current region
        region_key = $1 FS $2 FS $3;

        # Check if the region key exists in our calculated values map
        if (region_key in avg_signal) {
            # Print original region info + calculated values
            print $1, $2, $3, $4, $5, avg_signal[region_key], max_signal[region_key], max_scaled_signal[region_key];
        } else {
            # If the region had no overlaps, print original info + 0 for calculated values
            print $1, $2, $3, $4, $5, 0, 0, 0;
        }
    }
' "$CALCULATED_VALUES_BED" "$SORTED_REGIONS_BED" > "$AUGMENTED_REGIONS_BED"
echo "Merging complete. Augmented regions stored in $AUGMENTED_REGIONS_BED"

# --- Final Output Formatting ---

# Create the final CSV file with a header row
echo "Creating final CSV output file: $OUTPUT_CSV_FILE"
echo "Chromosome,Start,End,ID,Functional,SumSignal,MaxSignal,MaxScaledSignal" > "$OUTPUT_CSV_FILE"

# Append the data from the augmented regions file, converting tabs to commas
awk 'BEGIN { FS="\t"; OFS="," } { $1=$1; print }' "$AUGMENTED_REGIONS_BED" >> "$OUTPUT_CSV_FILE"
# The `$1=$1` trick forces awk to rebuild the record using the new OFS

echo "Script finished successfully. Output written to $OUTPUT_CSV_FILE"

# Cleanup is handled by the trap function upon exit

exit 0
