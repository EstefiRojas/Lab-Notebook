#!/bin/bash
# This script retrieves the minimum p-value of SNPs present at lncRNAs from 
# GWAS summary-statistics API using corresponding ENSEMBL coordinates as filter.
# It also filters by p-value < 1e-5.
#
# Dependencies: 
# - jq:
#       https://jqlang.github.io/jq/
#
# Script Name: get_GWAS_associations.sh
#
# Author: Estefania Rojas
#
# Input: 1. Csv file: containing the dataset of sequences to be analysed. This 
# 			file should have the exon 1 and 2 coordinates at columns 10-21.
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

    # Remove 'chr' prefix if present (case insensitive)
    chrm=$(echo "$chrm" | sed -E 's/^chr//i')
    
    # Convert to uppercase for consistency
    chrm=$(echo "$chrm" | tr '[:lower:]' '[:upper:]')
    
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

# Util function to return the smallest of two numbers
min() {
    echo $(( $1 < $2 ? $1 : $2 ))
}

get_unique_gwas_variants() {
   local chrm=$(convert_chromosome "$1" | tr -d '\r\n')
   local start_pos=$(echo "$2" | tr -d '\r\n')
   local end_pos=$(echo "$3" | tr -d '\r\n')

   # Define cache file location
   local cache_file="../data/cache/gwas_variants_cache.txt"
   
   # Generate cache key
   local cache_key=$(generate_cache_key "$chrm" "$start_pos" "$end_pos")
   
   # Try to get results from cache
   local cached_result=$(read_from_cache "$cache_file" "$cache_key")
   if [ ! -z "$cached_result" ]; then
       echo "Found in cache: ${cached_result}" >&2
       echo "$cached_result"
       return
   fi
   
   echo "Cache miss - querying API..." >&2

   local page_size=500
   local current_page_size="$page_size"
   local retry_offset=1
   
   # Create temporary files
   local variant_file=$(mktemp)
   local pvalue_file=$(mktemp)
   
   # Initialize current url
   local current_url="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/${chrm}/associations?bp_lower=${start_pos}&bp_upper=${end_pos}&p_upper=1e-5&size=${current_page_size}"
   local max_retries=1
   local page_count=0
   
   while true; do
       page_count=$((page_count + 1))
       echo "Processing page ${page_count}, URL: ${current_url}" >&2
       
       # Create temp file for headers
       local header_file=$(mktemp)
       local retry_count=0
       local success=false
       
       while [ $retry_count -lt $max_retries ] && [ "$success" = false ]; do
           # Query current page and save headers separately, following redirects
           local gwas_response=$(curl -s -L -D "${header_file}" --retry 1 "${current_url}" -H "Accept:application/json")
           
           # Get HTTP status from headers (get the last status code in case of redirects)
           local http_code=$(grep -i "^HTTP/" "${header_file}" | tail -n1 | cut -d' ' -f2)
           echo "HTTP response code: ${http_code}" >&2
           
           # Log the full response headers for debugging
           #echo "Response headers:" >&2
           #cat "${header_file}" >&2
           
           if [ "$http_code" = "200" ]; then
               success=true
               # If we were in retry mode (page size < 500), increment the page size
               if [ $current_page_size -lt $page_size ]; then
                   current_page_size=$(min $((current_page_size + 10)) $page_size)
                   echo "Successfully queried with reduced page size. Incrementing to: ${current_page_size}" >&2
               fi
           else
               retry_count=$((retry_count + 1))
               echo "Request failed with HTTP code ${http_code}, attempt ${retry_count}/${max_retries}" >&2
               [ $retry_count -lt $max_retries ] && sleep 0.1
           fi
       done
       
       rm "${header_file}"
       
       # If all retries failed, break the loop
       if [ "$success" = false ]; then
            message=$(echo "$gwas_response" | jq -r '.message // empty')
            if grep -q "Study GCST.*does not exist" <<< "$gwas_response"; then
                # Reduce page size to 10
                current_page_size=10
                current_url=$(echo "$current_url" | sed "s/size=[0-9]*/size=${current_page_size}/")

                # Extract current start value from URL
                current_start=$(echo "$current_url" | grep -o 'start=[0-9]*' | cut -d'=' -f2)
                if [ -z "$current_start" ]; then
                    current_start=0
                fi
                
                # Calculate new start value
                new_start=$((current_start + retry_offset))
                
                # Update URL with new start value
                if [[ "$current_url" =~ start= ]]; then
                    current_url=$(echo "$current_url" | sed "s/start=[0-9]*/start=${new_start}/")
                else
                    current_url="${current_url}&start=${new_start}"
                fi
                
                echo "GCST error encountered. Retrying with reduced page size=${current_page_size} and new start position: ${new_start}" >&2
                continue
            else
                echo "Warning: Failed to fetch data after $max_retries retries" >&2
                break
            fi
        fi
       
       # Check if JSON is valid before processing
       if ! echo "$gwas_response" | jq empty 2>/dev/null; then
           echo "Warning: Invalid JSON response received" >&2
           break
       fi
       
       # Count elements in current page
       local page_variants=$(echo "$gwas_response" | jq -r '._embedded.associations | length')
       echo "Variants in current page: ${page_variants}" >&2
       
       # Extract variant_ids and pvalues and append to files, counting how many were extracted
       local extracted_variants=$(echo "$gwas_response" | jq -r '._embedded.associations | to_entries[] | select(.value.p_value != -99) | .value.variant_id' | wc -l | tr -d ' ')
       local extracted_pvalues=$(echo "$gwas_response" | jq -r '._embedded.associations | to_entries[] | select(.value.p_value != -99) | .value.p_value' | wc -l | tr -d ' ')
       echo "Extracted ${extracted_variants} variants and ${extracted_pvalues} p-values after filtering" >&2
       
       if ! echo "$gwas_response" | jq -r '._embedded.associations | to_entries[] | select(.value.p_value != -99) | .value.variant_id' >> "$variant_file"; then
           echo "Warning: Error processing variant IDs" >&2
           break
       fi
       
       if ! echo "$gwas_response" | jq -r '._embedded.associations | to_entries[] | select(.value.p_value != -99) | .value.p_value' >> "$pvalue_file"; then
           echo "Warning: Error processing p-values" >&2
           break
       fi
       
       # Check for next page and handle URL properly
       local next_url=$(echo "$gwas_response" | jq -r '._links.next.href // empty')
       
       # Break if no next page
       if [ -z "$next_url" ]; then
           echo "No more pages to process" >&2
           break
       fi
       
       # Always use HTTPS
       next_url=$(echo "$next_url" | sed 's|^http:|https:|')
       
       # When getting next URL, use current_page_size
        next_url=$(echo "$next_url" | sed "s/size=[0-9]*/size=${current_page_size}/")
       
       # If the URL is relative (starts with /), construct the full URL
       if [[ "$next_url" == /* ]]; then
           next_url="https://www.ebi.ac.uk/gwas/summary-statistics/api${next_url}"
       elif [[ ! "$next_url" =~ ^https?:// ]]; then
           # If it's not a full URL and not starting with /, construct based on current URL
           base_url=$(echo "$current_url" | grep -o 'https://[^/]*/[^?]*')
           next_url="${base_url}${next_url}"
       fi
       
       # Validate the URL format
       if [[ ! "$next_url" =~ ^https://www\.ebi\.ac\.uk/gwas/summary-statistics/api/chromosomes/ ]]; then
           echo "Warning: Invalid next URL format" >&2
           break
       fi
       
       # Ensure the URL is properly encoded
       next_url=$(echo "$next_url" | sed 's/ /%20/g')
       echo "Next URL: ${next_url}" >&2
       
       current_url="$next_url"
       
       # Show current accumulated unique variants
       local current_unique=$(sort "$variant_file" | uniq | wc -l | tr -d ' ')
       echo "Accumulated unique variants so far: ${current_unique}" >&2
   done

   # Check if files are empty
   if [ ! -s "$variant_file" ]; then
       rm "$variant_file" "$pvalue_file"
       echo "0,NA"
       return
   fi

   # Count unique variants
   local unique_count=$(sort "$variant_file" | uniq | wc -l | tr -d ' ')
   
   # Find minimum p-value using awk with proper numeric comparison
   # The +$1 forces awk to treat the input as numeric, handling both scientific and decimal notation
   local min_p_value=$(awk '
       BEGIN { min=1 }
       NR==1 { min=+$1 }  # Initialize with first value
       {
           val = +$1      # Force numeric interpretation
           if (val > 0 && val < min) min = val
       }
       END { printf "%.2e\n", min }
   ' "$pvalue_file")

   # Clean up temporary files
   rm "$variant_file" "$pvalue_file"

   # Before returning results, write to cache
   local result="${unique_count},${min_p_value}"
   write_to_cache "$cache_file" "$cache_key" "$result"
   
   echo "$result"
}

# Simple progress bar function using built-in commands
show_progress() {
    local current=$1
    local total=$2
    local width=50  # Width of the progress bar
    local percentage=$((current * 100 / total))
    local completed=$((width * current / total))
    local remaining=$((width - completed))
    
    printf "\rProgress: ["
    printf "%${completed}s" | tr ' ' '#'
    printf "%${remaining}s" | tr ' ' '-'
    printf "] %d%% (%d/%d)" "$percentage" "$current" "$total"
}

# Function to generate a unique key for the cache
generate_cache_key() {
    local chrm="$1"
    local start="$2"
    local end="$3"
    echo "${chrm}_${start}_${end}"
}

# Function to read from cache
read_from_cache() {
    local cache_file="$1"
    local cache_key="$2"
    
    if [ -f "$cache_file" ]; then
        # Look for the key in cache and return the corresponding value
        grep "^${cache_key}:" "$cache_file" | cut -d':' -f2
    fi
}

# Function to write to cache
write_to_cache() {
    local cache_file="$1"
    local cache_key="$2"
    local value="$3"
    
    # Create cache directory if it doesn't exist
    local cache_dir=$(dirname "$cache_file")
    mkdir -p "$cache_dir"
    
    # Append to cache file
    echo "${cache_key}:${value}" >> "$cache_file"
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

echo "Output file: $output_path/$OUTPUT_NAME.csv"

# Add header to output file if new run
if [ ! -e $output_path/$OUTPUT_NAME.csv ]; then
	echo "tl_exon1_gwas_associations,tl_exon1_min_p_value,tl_exon2_gwas_associations,tl_exon2_min_p_value,tl_gwas_assoc,tl_min_p_val,tr_gwas_assoc,tr_min_p_val" > "$output_path"/"$OUTPUT_NAME".csv
fi

# Calculate total lines to process
total_lines=$(awk -F, 'NR > 1 { count++ } END { print count }' "$REGIONS_FILE")
echo "Total lines to process: $total_lines"
echo "Starting from line: $start_line"
processed_lines=$((start_line - 1))

# Loop through all ensembl coordinates
while IFS= read -r line || [[ -n "$line" ]]; do
    # Process the line and capture the output
    if processed_line=$(process_line "$line"); then
        # Read the processed fields into variables
        IFS=' ' read -r tl_exon1_chrm tl_exon1_start tl_exon1_end \
                     tl_exon2_chrm tl_exon2_start tl_exon2_end \
                     tl_chrm tl_start tl_end \
                     tr_chrm tr_start tr_end <<< "$processed_line"

		# Transcript 1, exon 1
        echo "Exon 1:" >&2
		if [[ "$tl_exon1_start" != "NA" ]]; then
			tl_total_associations_exon1=$(get_unique_gwas_variants $(echo "${tl_exon1_chrm}" | tr -d '"') ${tl_exon1_start} ${tl_exon1_end})
		else
			tl_total_associations_exon1='0,NA'
		fi
		#echo "Found $tl_total_associations_exon1 associations for transcript1 exon1"

		# Transcript 1, exon 2
        echo "Exon 2:" >&2
		if [[ "$tl_exon2_start" != "NA" ]]; then
			tl_total_associations_exon2=$(get_unique_gwas_variants $(echo "${tl_exon2_chrm}" | tr -d '"') ${tl_exon2_start} ${tl_exon2_end})
		else
			tl_total_associations_exon2='0,NA'
		fi
		#echo "Found $tl_total_associations_exon2 associations for transcript1 exon2"

		# Transcript 1
        echo "Transcript 1:" >&2
		if [[ "$tl_start" != "NA" ]]; then
			tl_total_associations=$(get_unique_gwas_variants $(echo "${tl_chrm}" | tr -d '"') ${tl_start} ${tl_end})
		else
			tl_total_associations='0,NA'
		fi
		#echo "Found $tl_total_associations associations for transcript2 exon2"

		# Transcript 2
        echo "Transcript 2:" >&2
		if [[ "$tr_start" != "NA" ]]; then
			tr_total_associations=$(get_unique_gwas_variants $(echo "${tr_chrm}" | tr -d '"') ${tr_start} ${tr_end})
		else
			tr_total_associations='0,NA'
		fi
		#echo "Found $tr_total_associations associations for transcript2 exon2"

	    # Send to output file
	    echo "${tl_total_associations_exon1},${tl_total_associations_exon2},${tl_total_associations},${tr_total_associations}" >> "$output_path/$OUTPUT_NAME.csv"

	    # Report progress
	    ((processed_lines++))
	    show_progress "$processed_lines" "$total_lines"
	fi

done < <(awk -F, -v start="$start_line" 'NR > start { print $0 }' "$REGIONS_FILE")
