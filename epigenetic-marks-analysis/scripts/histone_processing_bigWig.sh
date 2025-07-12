#!/opt/local/bin/bash
#
# Dependencies:
# - bigWigAverageOverBed: (Part of UCSC Kent Utilities) Required for calculating average signal from bigWig files over BED regions.
#   Installation: http://hgdownload.soe.ucsc.edu/admin/exe/
# - awk: Standard Unix text processing tool.
# - sort: Standard Unix sorting tool.
# - mktemp: Standard Unix temporary file creation tool (BSD syntax compatible).
#
# Script Name: histone_processing_bigwig.sh (Final Version)
#
# Author: Estefania Rojas (Modified for bigWig usage, custom temp dir, BSD mktemp, BED4 input)
#
# Description: This script calculates the average enrichment signal from histone
# marks (provided as bigWig files from ENCODE) for sequences defined in an input CSV file.
# It creates an intermediate BED4 file (chrom, start, end, name) for compatibility
# with bigWigAverageOverBed, avoiding potential issues with non-numeric score columns.
# Uses a configurable temporary directory with BSD mktemp compatibility.
#
# Input: 1. CSV file containing the sequences dataset to be analysed.
#           Format expected: Header line, then columns like ID,Functional,Chromosome,Start,End,...
#           e.g., data/functional-lncrna-exon1-dataset.csv
#        2. Name of the histone mark directory containing the bigWig files.
#           e.g., H3K36me3 (script will look in ../data/datasets/epigenetic_data/histone_marks/bigWig/H3K36me3/)
#        3. Base name for the output CSV file.
#           e.g., lncrna-exon1-histone-feature (will create H3K36me3_lncrna-exon1-histone-feature.csv)
#
############################################################################################################################

# Bash strict mode. -e exit on error, -u error on unset variables, -o pipefail ensures pipeline errors are caught.
# Use 'set -x' to trace execution if needed.
set -ueo pipefail
# set -x

# --- Configuration and Input Validation ---

# Get the date of execution. Format example: Tuesday April 03, 2025
DATE=$(date "+%A %B %d, %Y")
echo "Run started on: ${DATE}"
echo "=================================================="

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <regions_input.csv> <histone_name> <output_base_name>"
    echo "  Example: $0 data/regions.csv H3K36me3 region_features"
    exit 1
fi

# Store input arguments with descriptive names
REGIONS_CSV_FILE="$1"
HISTONE_NAME="$2"
OUTPUT_BASE_NAME="$3"

# --- Temporary Directory Configuration ---
# Define the directory for temporary files.
# Modify './tmp' to an absolute path if needed (e.g., /path/to/large/storage/tmp)
SCRIPT_TMP_DIR=../tmp

# Create the temporary directory if it doesn't exist
mkdir -p "$SCRIPT_TMP_DIR"
if [ ! -d "$SCRIPT_TMP_DIR" ] || [ ! -w "$SCRIPT_TMP_DIR" ]; then
    echo "Error: Temporary directory '$SCRIPT_TMP_DIR' could not be created or is not writable."
    exit 1
fi
echo "Using temporary directory: $SCRIPT_TMP_DIR"
# --- End Temporary Directory Configuration ---


# Define base paths (adjust if your directory structure is different)
EPIGENETIC_DATA_BASE_DIR="../data/datasets/epigenetic_data/histone_marks/bigWig"
OUTPUT_FEATURES_BASE_DIR="../data/datasets/histone_feature/bigWig" # Changed from broadPeak to bigWig

# Construct specific paths based on histone name
MARKS_PATH="${EPIGENETIC_DATA_BASE_DIR}/${HISTONE_NAME}"
OUTPUT_DIR="${OUTPUT_FEATURES_BASE_DIR}/${HISTONE_NAME}"
FINAL_OUTPUT_CSV="${OUTPUT_DIR}/${HISTONE_NAME}_${OUTPUT_BASE_NAME}.csv"

# Validate input file and directory existence
if [ ! -f "$REGIONS_CSV_FILE" ]; then
    echo "Error: Input regions CSV file not found: '$REGIONS_CSV_FILE'"
    exit 1
fi

if [ ! -d "$MARKS_PATH" ]; then
    echo "Error: Histone marks directory not found: '$MARKS_PATH'"
    exit 1
fi

# --- Temporary Files ---
# Use mktemp with BSD-compatible syntax: mktemp directory/prefix.XXXXXX
# Note: Suffixes are part of the prefix here, not a separate option.
sorted_regions_bed=$(mktemp "${SCRIPT_TMP_DIR}/sorted_regions.XXXXXX")
aggregated_signals_tmp=$(mktemp "${SCRIPT_TMP_DIR}/agg_signals.XXXXXX")
final_augmented_bed=$(mktemp "${SCRIPT_TMP_DIR}/aug_regions.XXXXXX")
bw_output_tmp=$(mktemp "${SCRIPT_TMP_DIR}/bw_output.XXXXXX") # Temporary output for bigWigAverageOverBed

# Setup trap to ensure temporary files are cleaned up on exit (normal or error)
# The trap works correctly as it uses the variables holding the full paths to the temp files.
# Added .avg file cleanup
trap 'echo "Cleaning up temporary files..."; rm -f "$sorted_regions_bed" "$aggregated_signals_tmp" "$final_augmented_bed" "$bw_output_tmp" "${aggregated_signals_tmp}.avg"; echo "Cleanup complete."' EXIT


# --- Find bigWig Files (Using mapfile for Bash 4.0+) ---
echo "Finding bigWig files in '$MARKS_PATH'..."
# Use mapfile and process substitution to read find results into array
mapfile -t BIGWIG_FILES < <(find "$MARKS_PATH" -type f -name "*.bigWig")

if [ ${#BIGWIG_FILES[@]} -eq 0 ]; then
    echo "Error: No .bigWig files found in '$MARKS_PATH'"
    exit 1
fi
echo "Found ${#BIGWIG_FILES[@]} bigWig file(s) for histone mark '$HISTONE_NAME':"
printf "  %s\n" "${BIGWIG_FILES[@]}" # printf is safe and works across versions
# --- End Finding bigWig Files ---


# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
echo "Output will be saved to: '$FINAL_OUTPUT_CSV'"


# --- Data Processing ---

echo "--------------------------------------------------"
echo "1. Preparing and sorting regions from '$REGIONS_CSV_FILE'..."
# Extract relevant columns (Chromosome, Start, End, ID, Functional) from CSV (skip header)
# Output as BED4: chr, start, end, name (ID)
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1}' OFS='\t' "$REGIONS_CSV_FILE" | \
    sort -k1,1 -k2,2n > "$sorted_regions_bed"
echo "  Regions prepared and saved to temporary file: $sorted_regions_bed"
REGION_COUNT=$(wc -l < "$sorted_regions_bed")
echo "  Processed $REGION_COUNT regions."

echo "--------------------------------------------------"
echo "2. Calculating average signal across ${#BIGWIG_FILES[@]} bigWig file(s)..."
# Initialize temporary file for aggregation
> "$aggregated_signals_tmp"

# Process each bigWig file
file_index=0
for bw_file in "${BIGWIG_FILES[@]}"; do
    file_index=$((file_index + 1))
    echo "  Processing bigWig file ($file_index/${#BIGWIG_FILES[@]}): $(basename "$bw_file")"

    # Run bigWigAverageOverBed
    # Output format: name(region ID), size, covered, sum, mean0(over size), mean(over covered)
    if bigWigAverageOverBed "$bw_file" "$sorted_regions_bed" "$bw_output_tmp"; then
        # Append the results (region ID and mean signal) to the aggregation file
        # We use the 'mean' value (column 6) which is sum/coveredBases
        awk -v OFS='\t' '{print $1, $6}' "$bw_output_tmp" >> "$aggregated_signals_tmp"
        echo "    Successfully processed."
    else
        echo "    Warning: bigWigAverageOverBed failed for file $(basename "$bw_file"). Skipping this file."
        # Consider whether failure should be fatal or just a warning
    fi
done

echo "--------------------------------------------------"
echo "3. Aggregating signals and joining with regions..."
# Aggregate the signals using awk: calculate the maximum signal for each region ID
# Input to awk ($aggregated_signals_tmp): region_id <tab> mean_signal (repeated for each bigWig)
# Output of awk: region_id <tab> average_mean_signal
awk '
BEGIN { OFS="\t" }
{
    region_id = $1
    signal = $2

    if (!(region_id in max_signal) || signal > max_signal[region_id]) {
        max_signal[region_id] = signal;
    }

}
END {
    for (id in max_signal) {
        # Select max signal
        print id, max_signal[id];
    }
}' "$aggregated_signals_tmp" | sort -k1,1 > "${aggregated_signals_tmp}.avg" # Sort by ID for potential joining efficiency

echo "  Joining aggregated signals back to the original regions..."
# Join the averaged signal back to the sorted regions BED file
# Use awk to perform a left join: print all regions, adding the signal if found, otherwise add 0
awk -F'\t' '
BEGIN { OFS="\t" }
NR==FNR { # Process the first file (average signals)
    signals[$1] = $2 # Store signal keyed by region ID
    next
}
{ # Process the second file (sorted regions BED)
    region_id = $4 # The ID is in the 4th column of our BED file
    signal_value = (region_id in signals) ? signals[region_id] : NA # Default to NA if no signal found
    print $1, $2, $3, $4, $5, signal_value # Print original BED fields + signal
}' "${aggregated_signals_tmp}.avg" "$sorted_regions_bed" > "$final_augmented_bed"

echo "--------------------------------------------------"
echo "4. Formatting final output CSV file..."
# Add header row to the final output CSV file
# Columns: Chromosome, Start, End, ID, Functional, AvgSignal(<HistoneName>)
echo "Chromosome,Start,End,ID,Functional,${HISTONE_NAME}_AvgSignal" > "$FINAL_OUTPUT_CSV"

# Convert the final augmented BED (tab-separated) to CSV and append
sed 's/\t/,/g' "$final_augmented_bed" >> "$FINAL_OUTPUT_CSV"

echo "=================================================="
echo "Processing complete."
echo "Final output CSV file created: $FINAL_OUTPUT_CSV"
echo "Run finished on: $(date '+%A %B %d, %Y %H:%M:%S')"
echo "=================================================="

# Trap takes care of cleanup (including the .avg file), but explicit exit 0 for success
exit 0
