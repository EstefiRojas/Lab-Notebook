#!/opt/local/bin/bash
#
# Dependencies: 
# - bedtools:
#       https://bedtools.readthedocs.io/en/latest/content/installation.html
#
# Script Name: chrm_acc_processing.sh
#
# Author: Estefania Rojas
#
# Description: This script calculates the weighted sum of enrichment for sequences of functional long and short-ncRNA 
# databases and protein-coding-RNA databases for chromatin accessibility. The script computes the intersection between the sequences
# and chromatin accessibility marks reported in ENCODE (bed narrowPeak or broadPeak file format) using bedtools intersect.
#
# Input: 1. Csv file containing the sequeces dataset to be analysed. Ex. data/functional-lncrna-exon1-dataset.csv
#        2. Name of the output csv file to save weighted sum of enrichment feature for each sequence in input file 1.
#           Ex. lncrna-exon1-histone-feature
#
#
# ENCODE narrowPeak BED Format Reference:
# 1. chromosome
# 2. start
# 3. end
# 4. name (optional)
# 5. score (indicates statistical significance)
# 6. strand
# 7. signalValue (indicates magnitude of signal)
# 8. pValue
# 9. qValue
# 10. peak (Offset of the peak summit from start coordinate, 0-based, -1 if not available)
#
#
# ENCODE broadPeak BED Format Reference:
# 1. chromosome
# 2. start
# 3. end
# 4. name (optional)
# 5. score (indicates statistical significance)
# 6. strand
# 7. signalValue (indicates magnitude of signal)
# 8. pValue
# 9. qValue
#
############################################################################################################################

# Bash strict mode. -e exit on error, -u error on unset variables, -x print every command executed.
set -uex

# Get the date of execution. Format example: Tuesday January 17, 2024
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}
##############################################################


# Check for correct usage
if [ $# -ne 2 ]; then
    echo "Usage: $0 regions_to_extract.csv outname"
    exit 1
fi

# Store input files with descriptive names
REGIONS_FILE=$1
#EXPERIMENT_TYPE=$2
#BIOSAMPLE_CLASS=$3
OUTPUT_NAME=$2


# Temporary files
#temp_chrm_acc_marks_bed=$(mktemp)  # Temporary file for BED conversion
#sorted_chrm_acc_marks_bed=$(mktemp)

# Variables
data_path=../data/datasets/chrm_acc_feature/narrowPeak
marks_path=../data/datasets/epigenetic_data/chromatin_accessibility/narrowPeak
tmp_path="${data_path}/tmp"

# Create tmp directory
mkdir -p "${tmp_path}"

# Define temporary files with process ID for uniqueness
sorted_regions_file_bed="${tmp_path}/sorted_regions_$$.bed"
augmented_regions_bed="${tmp_path}/augmented_regions_$$.bed"


# Create a directory if not present already
if [ ! -d "$data_path" ]; then
    mkdir -p "$data_path"
fi

export TMPDIR="${tmp_path}"
# Remove header from regions file and sort regions
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"


# Heavy processing ahead. Check for previously created files and skip if present.
if [ ! -f "$data_path"/chrm_acc_sorted_marks.bed ]; then
    if [ ! -f "$marks_path"/temp_chrm_acc_marks.bed ]; then
        touch "$marks_path"/temp_chrm_acc_marks.bed
        for file in $(find "$marks_path" -type f -mindepth 2 -maxdepth 3 -print | grep .bed); do
            #gunzip "$file"
            #filename=$(basename "$file" .bed)  # Extract filename without .bed extension
            awk -v OFS="\t" '{print $0}' "${file}" >> "$marks_path"/temp_chrm_acc_marks.bed
        done
    fi

    awk -F '\t' '{
             print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9
         }' OFS=, "$marks_path"/temp_chrm_acc_marks.bed | \
    sort -k1,1 -k2,2n -k3,3n > "$data_path"/chrm_acc_sorted_marks.bed
fi


# Perform interval intersection using BEDTools
bedtools intersect -a "$data_path"/chrm_acc_sorted_marks.bed -b "$sorted_regions_file_bed" -wo > overlaps_chrm_acc_"$OUTPUT_NAME".bed

## AvgSignal = sum(wi * ei)
## Calculate the weighted sum of enrichment, grouped by Chromosome, StartExon, EndExon
#awk -F '\t' '
#    BEGIN { OFS="\t" }
#    {
#        key = $1 FS $8 FS $9;
#        sum[key] += $4 * $12 / ($9 - $8);
#   }
#    END {
#        for (k in sum) {
#            print k, sum[k];
#        }
#    }
#' overlaps_chrm_acc_"$OUTPUT_NAME".bed | sort -k1,1 -k2,2n -k3,3n > "$OUTPUT_NAME".bed

# Calculate maximum enrichment using signalValue, grouped by Chromosome, StartExon, EndExon
# Formula: weighted_sum = sum(signal_i * overlap_length_i) / region_length
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
' overlaps_chrm_acc_"$OUTPUT_NAME".bed | sort -k1,1 -k2,2n -k3,3n > "$OUTPUT_NAME".bed


# Join results
# Note: Regions without intersections are assigned a value of 0,
# which assumes these regions are not accessibile.
# This may affect downstream statistical analyses.
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
    }' \
        "$OUTPUT_NAME".bed "$sorted_regions_file_bed" > "$augmented_regions_bed"

# Add header row
echo Chromosome','Start','End','ID','Functional','chrm_acc_SumSignal','chrm_acc_MaxSignal','chrm_acc_MaxScaledSignal > "$data_path"/chrm_acc_"$OUTPUT_NAME".csv

# Convert BEDTools output to desired CSV format
sed 's/\t/,/g' "$augmented_regions_bed" >> "$data_path"/chrm_acc_"$OUTPUT_NAME".csv


# Cleanup is now handled by the trap
# Add cleanup trap
cleanup() {
    rm -f "${sorted_regions_file_bed}"
    rm -f "${augmented_regions_bed}"
    rm -f "overlaps_chrm_acc_${OUTPUT_NAME}.bed"
    rm -f "${OUTPUT_NAME}.bed"
}
trap cleanup EXIT