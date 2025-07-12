#!/bin/bash
# This script retrieves the minimum p-value of SNPs present at lncRNAs from 
# GWAS summary-statistics API using corresponding ENSEMBL coordinates as filter.
# It also filters by p-value < 1e-5.
#
# Dependencies: 
# - jq:
#       https://jqlang.github.io/jq/
# - parallel:
#       https://www.gnu.org/software/parallel/
#
# Script Name: get_GWAS_associations.sh
#
# Author: Estefania Rojas
#
# Input: 1. Csv file: containing the dataset of sequences to be analysed. This 
# 			file should have the exon 1 and 2 coordinates at columns 3-9.
#           Ex. ../data/lncRNAs_ensembl_exon_coords.csv
#        2. Name of the output csv file: minimum p-value for each coordinate in 
# 			input file 1.
#           Ex. lncrna-gwas-pval-feature
#
###############################################################################
# Bash strict mode. -e report all errors, -u exit on error, -x print every command executed.
#set -x

# Get the date of execution. Format example: Tuesday January 17, 2024
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}

# Check for correct usage
if [ $# -ne 2 ]; then
    echo "Usage: $0 regions_to_extract.csv outname.csv"
    exit 1
fi
##############################################################
#Function to convert X to 23 and Y to 24 while leaving numeric values unchanged
convert_chromosome() {
    local chrm=$1
    
    case "$chrm" in
        "X"|"x")
            echo "23"
            ;;
        "Y"|"y")
            echo "24"
            ;;
        *)
            echo "$chrm"
            ;;
    esac
}

# Function to validate and process a line of data
process_line() {
    local line="$1"
    local -a fields
    
    # Read fields into array, handling empty fields
    IFS=',' read -ra fields <<< "$line"
    
    # Ensure array has at least 12 elements, pad with empty strings if needed
    while [[ ${#fields[@]} -lt 12 ]]; do
        fields+=("")
    done
    
    # Extract fields (using empty string if field is unset)
    local tl_exon1_chrm="${fields[9]}"
    local tl_exon1_start="${fields[10]}"
    local tl_exon1_end="${fields[11]}"
    local tl_exon2_chrm="${fields[12]}"
    local tl_exon2_start="${fields[13]}"
    local tl_exon2_end="${fields[14]}"
    local tl_chrm="${fields[15]}"
    local tl_start="${fields[16]}"
    local tl_end="${fields[17]}"
    local tr_chrm="${fields[18]}"
    local tr_start="${fields[19]}"
    local tr_end="${fields[20]}"
    
    # Output fields space-separated, using NULL for empty values
    echo "${tl_exon1_chrm:-NULL} ${tl_exon1_start:-0} ${tl_exon1_end:-0}" \
         "${tl_exon2_chrm:-NULL} ${tl_exon2_start:-0} ${tl_exon2_end:-0}" \
         "${tl_chrm:-NULL} ${tl_start:-0} ${tl_end:-0}" \
         "${tr_chrm:-NULL} ${tr_start:-0} ${tr_end:-0}"
}

# Function to get GWAS variants with rate limiting and better error handling
get_unique_gwas_variants() {
   local chrm=$(convert_chromosome "$1")
   local start_pos=$2
   local end_pos=$3
   local retry_count=0
   local max_retries=3
   local wait_time=2
   
   # Initialize empty array for storing variant IDs and p-values
   local variant_ids=()
   local p_values=()
   
   # Initialize current url
   local current_url="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/${chrm}/associations?bp_lower=${start_pos}&bp_upper=${end_pos}&p_upper=1e-5&size=40"
   
   while true; do
       # Rate limiting
       sleep 0.5
       
       # Create temp file for headers
       local header_file=$(mktemp)
       
       # Query current page with retry logic
       local gwas_response=""
       while [ $retry_count -lt $max_retries ]; do
           gwas_response=$(curl -s -D "${header_file}" --retry 2 --max-time 30 "${current_url}" -H "Accept:application/json")
           local http_code=$(grep "HTTP/" "${header_file}" | cut -d' ' -f2)
           
           if [ "$http_code" = "200" ]; then
               break
           elif [ "$http_code" = "429" ]; then
               # Rate limit hit - wait longer
               wait_time=$((wait_time * 2))
               sleep $wait_time
           else
               # Other error - retry with backoff
               sleep $wait_time
           fi
           ((retry_count++))
       done
       
       rm "${header_file}"
       
       # Check for errors after retries
       if [ $retry_count -eq $max_retries ]; then
           echo "0,NA"
           return
       fi
       
       # Extract variant_ids and pvalues
       local page_variants=$(echo "$gwas_response" | jq -r '._embedded.associations | to_entries[] | .value.variant_id' 2>/dev/null)
       local page_pvalues=$(echo "$gwas_response" | jq -r '._embedded.associations | to_entries[] | .value.p_value' 2>/dev/null)
       
       # Convert to arrays
       IFS=' ' read -ra new_ids <<< "$page_variants"
       IFS=' ' read -ra new_pvals <<< "$page_pvalues"
       
       # Filter out -99 values
       for i in "${!new_pvals[@]}"; do
           if [[ "${new_pvals[$i]}" != "-99" ]]; then
               variant_ids+=("${new_ids[$i]}")
               p_values+=("${new_pvals[$i]}")
           fi
       done
       
       # Check for next page
       local next_url=$(echo "$gwas_response" | jq -r '._links.next.href // empty')
       
       # Break if no next page
       if [ -z "$next_url" ]; then
           break
       fi
       
       current_url=$next_url
   done

   # Check if arrays are empty
   if [ ${#variant_ids[@]} -eq 0 ]; then
       echo "0,NA"
       return
   fi

   # Count unique variants and get minimum p-value
   local unique_count=$(printf '%s\n' "${variant_ids[@]}" | sort | uniq | wc -l | tr -d ' ')
   local min_p_value=$(printf '%s\n' "${p_values[@]}" | sort -g | head -n1)
   
   echo "${unique_count},${min_p_value}"
}

# Store input variables
REGIONS_FILE=$1
OUTPUT_NAME=$2

echo "Processing file $REGIONS_FILE"

# Define path to store output
output_path=../data/datasets/gwas_pval_feature

# Check previous runs
echo "Checking for previous runs..."
start_line=1
if [ -w $output_path/$OUTPUT_NAME.csv ]; then
	echo "Found a previous run with specified parameters."
	start_line=$(wc -l $output_path/$OUTPUT_NAME.csv | awk '{ print $1 }')
fi

# Create output directory if needed
mkdir -p "$output_path"
echo "Created/verified output directory: $output_path"

echo "Output file will be: $output_path/$OUTPUT_NAME.csv"

# Add header to output file if new run
if [ ! -e "$output_path/$OUTPUT_NAME.csv" ]; then
    echo "Creating new output file with header..."
	echo "tl_exon1_gwas_associations,tl_exon1_min_p_value,tl_exon2_gwas_associations,tl_exon2_min_p_value,tl_gwas_assoc,tl_min_p_val,tr_gwas_assoc,tr_min_p_val" > "$output_path/$OUTPUT_NAME.csv"
fi

# Export functions for parallel processing
echo "Exporting functions for parallel processing..."
export -f convert_chromosome
export -f process_line
export -f get_unique_gwas_variants

# Create a temporary directory for parallel processing
TEMP_DIR=$(mktemp -d)
echo "Created temporary directory: $TEMP_DIR"
trap 'rm -rf "$TEMP_DIR"' EXIT

# Split the input file into chunks for parallel processing
echo "Splitting input file into chunks..."
split -l 100 <(tail -n +2 "$REGIONS_FILE") "$TEMP_DIR/chunk_"
chunk_count=$(ls "$TEMP_DIR"/chunk_* | wc -l)
echo "Created $chunk_count chunks"

# Process chunks in parallel
echo "Starting parallel processing of chunks..."

# First, let's verify the chunk files exist and have content
echo "Checking chunk files:"
ls -l "$TEMP_DIR"/chunk_*
echo "Sample of first chunk:"
head -n 2 "$TEMP_DIR"/chunk_aa

parallel --jobs 4 --verbose '
    echo "Processing chunk: {}"
    echo "First few lines of chunk {}:"
    head -n 2 {}
    while IFS= read -r line || [[ -n "$line" ]]; do
        echo "Processing line: $line"
        if processed_line=$(process_line "$line"); then
            echo "Successfully processed line into: $processed_line"
            IFS=" " read -r tl_exon1_chrm tl_exon1_start tl_exon1_end \
                          tl_exon2_chrm tl_exon2_start tl_exon2_end \
                          tl_chrm tl_start tl_end \
                          tr_chrm tr_start tr_end <<< "$processed_line"
            
            if [[ "$tl_exon1_start" != "NA" ]]; then
                tl_total_associations_exon1=$(get_unique_gwas_variants $(echo "${tl_exon1_chrm}" | tr -d "\"") ${tl_exon1_start} ${tl_exon1_end})
            else
                tl_total_associations_exon1="0,NA"
            fi

            if [[ "$tl_exon2_start" != "NA" ]]; then
                tl_total_associations_exon2=$(get_unique_gwas_variants $(echo "${tl_exon2_chrm}" | tr -d "\"") ${tl_exon2_start} ${tl_exon2_end})
            else
                tl_total_associations_exon2="0,NA"
            fi

            if [[ "$tl_start" != "NA" ]]; then
                tl_total_associations=$(get_unique_gwas_variants $(echo "${tl_chrm}" | tr -d "\"") ${tl_start} ${tl_end})
            else
                tl_total_associations="0,NA"
            fi

            if [[ "$tr_start" != "NA" ]]; then
                tr_total_associations=$(get_unique_gwas_variants $(echo "${tr_chrm}" | tr -d "\"") ${tr_start} ${tr_end})
            else
                tr_total_associations="0,NA"
            fi

            echo "${tl_total_associations_exon1},${tl_total_associations_exon2},${tl_total_associations},${tr_total_associations}" >> "$TEMP_DIR/results_$(basename {})"
            echo "Wrote results for line"
        else
            echo "Failed to process line: $line"
        fi
    done < {}' ::: "$TEMP_DIR"/chunk_*

echo "Parallel processing complete. Combining results..."

# Combine results in order
if [ ! -e "$output_path/$OUTPUT_NAME.csv" ]; then
    echo "Creating final output file with header..."
    echo "tl_exon1_gwas_associations,tl_exon1_min_p_value,tl_exon2_gwas_associations,tl_exon2_min_p_value,tl_gwas_assoc,tl_min_p_val,tr_gwas_assoc,tr_min_p_val" > "$output_path/$OUTPUT_NAME.csv"
fi

echo "Concatenating chunk results..."
cat "$TEMP_DIR"/results_chunk_* >> "$output_path/$OUTPUT_NAME.csv"

# Clean up
echo "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

echo "Script completed successfully."
echo "Results written to: $output_path/$OUTPUT_NAME.csv"
