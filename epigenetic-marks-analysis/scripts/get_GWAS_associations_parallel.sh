#!/bin/bash
# This script retrieves GWAS association data in parallel, using a file-based
# cache and a rate-limiter to avoid re-querying or overloading the API.
#
# Dependencies:
# - jq:       https://jqlang.github.io/jq/
# - parallel: https://www.gnu.org/software/parallel/
#
# Author: Estefania Rojas
#
###############################################################################

# Get the date of execution
DATE=$(date "+%A %B %d, %Y")
echo "Execution date: ${DATE}"

# Default value for verbose
VERBOSE=false

# Parse options
while getopts "v" opt; do
  case ${opt} in
    v )
      VERBOSE=true
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

# Now, the remaining arguments are the positional ones.
#if [ "$VERBOSE" = true ]; then
#    set -x
#fi

# Check for correct usage
if [ $# -ne 2 ]; then
    echo "Usage: $0 [-v] <regions_to_extract.csv> <outname.csv>"
    exit 1
fi

##############################################################
# HELPER FUNCTIONS
##############################################################

# Function to sanitize chromosome names
convert_chromosome() {
    local chrm=$1
    # Remove 'chr' prefix (case-insensitive) and convert to uppercase
    chrm=$(echo "$chrm" | sed -E 's/^chr//i' | tr '[:lower:]' '[:upper:]')
    case "$chrm" in
        "X") echo "23" ;;
        "Y") echo "24" ;;
        *) echo "$chrm" ;;
    esac
}

# Util function to return the smallest of two numbers
min() {
    echo $(( $1 < $2 ? $1 : $2 ))
}

# Function to generate a unique key for caching
generate_cache_key() {
    echo "$1_$2_$3"
}

# Function to read a result from the cache file
read_from_cache() {
    local cache_file="$1"
    local cache_key="$2"
    if [ -f "$cache_file" ]; then
        grep "^${cache_key}:" "$cache_file" | cut -d':' -f2-
    fi
}

# Function to write a result to the cache file
write_to_cache() {
    local cache_file="$1"
    local cache_key="$2"
    local value="$3"
    local cache_dir
    cache_dir=$(dirname "$cache_file")
    mkdir -p "$cache_dir"
    echo "${cache_key}:${value}" >> "$cache_file"
}

# A portable and race-condition-safe rate-limiting function.
throttle_api_call() {
    local max_calls_per_second=10
    local lock_dir="$LOCK_DIR"
    local stat_flags
    
    # Determine the correct flags for `stat` based on the operating system
    if [[ "$(uname)" == "Darwin" ]]; then
        stat_flags="-f %m" # macOS/BSD `stat`
    else
        stat_flags="-c %Y" # Linux/GNU `stat`
    fi

    while true; do
        local now
        now=$(date +%s)
        
        # Manually clean up lock files older than 1 second
        for token_file in "$lock_dir"/*; do
            if [ -f "$token_file" ]; then
                local file_epoch
                # Redirect stderr to /dev/null to silence "No such file" errors
                file_epoch=$(stat $stat_flags "$token_file" 2>/dev/null)
                
                # Only proceed if stat was successful
                if [ -n "$file_epoch" ]; then
                    if (( now - file_epoch > 1 )); then
                        rm -f "$token_file"
                    fi
                fi
            fi
        done

        # Count current lock files (tokens) in a robust way
        local current_calls
        current_calls=$(find "$lock_dir" -type f | wc -l)
        
        if [ "$current_calls" -lt "$max_calls_per_second" ]; then
            # If below the limit, create a new token and proceed
            touch "${lock_dir}/$(date +%s.%N)"
            break
        else
            # If the limit is reached, wait for a fraction of a second
            sleep 0.1
        fi
    done
}


# Main function to get GWAS data, i.e. number of unique variants, min p-value, odds-ratio, max odds-ratio, and min odds-ratio.
get_unique_gwas_variants() {
    local chrm
    local start_pos
    local end_pos
    chrm=$(convert_chromosome "$1")
    start_pos=$2
    end_pos=$3

    # Return immediately if input is invalid
    if [[ -z "$chrm" || "$chrm" == "NULL" || -z "$start_pos" || "$start_pos" == "0" ]]; then
        if [ "$VERBOSE" = true ]; then echo "Invalid input for get_unique_gwas_variants: chrm='$chrm', start_pos='$start_pos', end_pos='$end_pos'. Skipping." >&2; fi
        echo "NA,NA,NA,NA,NA"
        return
    fi

    # Define cache file and key
    local cache_file="../data/cache/gwas_variants_cache.txt"
    local cache_key
    cache_key=$(generate_cache_key "$chrm" "$start_pos" "$end_pos")

    # Try to get results from cache
    local cached_result
    cached_result=$(read_from_cache "$cache_file" "$cache_key")
    if [[ -n "$cached_result" ]]; then
        if [ "$VERBOSE" = true ]; then echo "Cache hit for key '$cache_key'. Result: $cached_result" >&2; fi
        echo "$cached_result"
        return
    fi
    if [ "$VERBOSE" = true ]; then echo "Cache miss for key '$cache_key'. Querying API." >&2; fi

    local data_file
    data_file=$(mktemp)
    trap 'rm -f "$data_file"' RETURN

    # API query logic with dynamic page size adjustment
    local max_page_size=500
    local current_page_size=$max_page_size
    local current_url="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/${chrm}/associations?bp_lower=${start_pos}&bp_upper=${end_pos}&p_upper=1e-8&size=${current_page_size}"

    while true; do
        if [ "$VERBOSE" = true ]; then echo "Calling API with URL: $current_url" >&2; fi
        # Call the rate-limiter before every API call
        throttle_api_call

        local gwas_response
        local header_file
        header_file=$(mktemp)
        gwas_response=$(curl -s -L -D "$header_file" --retry 2 --max-time 60 "${current_url}" -H "Accept:application/json")
        local http_code
        http_code=$(grep -i "^HTTP/" "${header_file}" | tail -n1 | cut -d' ' -f2)
        rm "$header_file"
        if [ "$VERBOSE" = true ]; then echo "API response code: $http_code" >&2; fi

        if [ "$http_code" = "200" ]; then
            # First, validate that the response is well-formed JSON.
            if ! echo "$gwas_response" | jq empty 2>/dev/null; then
                if [ "$VERBOSE" = true ]; then echo "Invalid JSON response received. Breaking loop." >&2; fi
                break
            fi

            # Now that we know it's valid JSON, process it.
            echo "$gwas_response" | jq -r '(._embedded.associations // [])[] | select(try (.p_value|tonumber) catch false) | select(.p_value != -99) | [.variant_id, .p_value, .odds_ratio] | @tsv' >> "$data_file"
            
            # Gradually increase page size back to max if it was reduced
            if [ "$current_page_size" -lt "$max_page_size" ]; then
                current_page_size=$(min $((current_page_size + 10)) $max_page_size)
            fi

            # Get the next page URL
            local next_url
            next_url=$(echo "$gwas_response" | jq -r '._links.next.href // empty')
            if [[ -z "$next_url" ]]; then
                if [ "$VERBOSE" = true ]; then echo "No next page. Exiting API loop." >&2; fi
                break # No more pages, exit loop
            fi
            current_url=$(echo "$next_url" | sed "s/size=[0-9]*/size=${current_page_size}/")
        
        # Handle "Too Many Requests" error from the API
        elif [ "$http_code" = "429" ]; then
            if [ "$VERBOSE" = true ]; then echo "API rate limit hit (429). Waiting 60 seconds." >&2; fi
            # Wait for 60 seconds before the next attempt.
            sleep 60
            continue # Retry the same URL

        # Check for the specific API error that requires intervention
        elif echo "$gwas_response" | grep -q "Study GCST.*does not exist"; then
            if [ "$VERBOSE" = true ]; then echo "Specific API error detected. Reducing page size and retrying." >&2; fi
            local retry_offset=1
            current_page_size=1 # Drastically reduce page size
            local current_start
            current_start=$(echo "$current_url" | grep -o 'start=[0-9]*' | cut -d'=' -f2)
            [ -z "$current_start" ] && current_start=0
            local new_start=$((current_start + retry_offset))
            
            # Update URL with new start and size
            if echo "$current_url" | grep -q 'start='; then
                 current_url=$(echo "$current_url" | sed -e "s/size=[0-9]*/size=${current_page_size}/" -e "s/start=[0-9]*/start=${new_start}/")
            else
                 current_url="${current_url}&start=${new_start}"
                 current_url=$(echo "$current_url" | sed "s/size=[0-9]*/size=${current_page_size}/")
            fi
            # Continue to retry with the new URL
            continue
        else
            if [ "$VERBOSE" = true ]; then echo "Unhandled API error ($http_code). Breaking loop." >&2; fi
            # For any other error, stop trying for this region
            break
        fi
    done

    # Process the collected data if the file is not empty
    if [ ! -s "$data_file" ]; then
        if [ "$VERBOSE" = true ]; then echo "No data collected for key '$cache_key'. Returning NA." >&2; fi
        echo "NA,NA,NA,NA,NA"
        return
    fi
    
    local unique_count
    unique_count=$(cut -f1 "$data_file" | sort -u | wc -l | tr -d ' ')

    local min_p_line
    min_p_line=$(cut -f2,3 "$data_file" | sort -t$'\t' -k1,1g | head -n1)
    
    local min_p_value="NA"
    local odds_ratio_value="NA"

    if [[ -n "$min_p_line" ]]; then
        IFS=$'\t' read -r min_p_value or_val <<< "$min_p_line"
        if [[ "$or_val" == "null" || -z "$or_val" ]]; then
            odds_ratio_value="NA"
        else
            odds_ratio_value="$or_val"
        fi
    fi

    local odds_ratio_max=$(cut -f3 "$data_file" | sort -k1,1g -r | head -n1)
    local odds_ratio_min=$(cut -f3 "$data_file" | sort -k1,1g | head -n1)

    local final_result="${unique_count},${min_p_value},${odds_ratio_value},${odds_ratio_min},${odds_ratio_max}"
    if [ "$VERBOSE" = true ]; then echo "Final result for key '$cache_key': $final_result" >&2; fi
    
    write_to_cache "$cache_file" "$cache_key" "$final_result"
    if [ "$VERBOSE" = true ]; then echo "Wrote result for key '$cache_key' to cache." >&2; fi
    
    echo "$final_result"
}

# Wrapper function to process a single line from the input CSV in parallel
process_gwas_region_line() {
    local line="$1"
    if [ "$VERBOSE" = true ]; then echo "Processing line: $line" >&2; fi
    local -a fields
    IFS=',' read -r -a fields <<< "$line"
    
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

    local tl_exon1_results
    local tl_exon2_results
    local tl_results
    local tr_results
    tl_exon1_results=$(get_unique_gwas_variants "$tl_exon1_chrm" "$tl_exon1_start" "$tl_exon1_end")
    if [ "$VERBOSE" = true ]; then echo "Result for tl_exon1: $tl_exon1_results" >&2; fi
    tl_exon2_results=$(get_unique_gwas_variants "$tl_exon2_chrm" "$tl_exon2_start" "$tl_exon2_end")
    if [ "$VERBOSE" = true ]; then echo "Result for tl_exon2: $tl_exon2_results" >&2; fi
    tl_results=$(get_unique_gwas_variants "$tl_chrm" "$tl_start" "$tl_end")
    if [ "$VERBOSE" = true ]; then echo "Result for tl: $tl_results" >&2; fi
    tr_results=$(get_unique_gwas_variants "$tr_chrm" "$tr_start" "$tr_end")
    if [ "$VERBOSE" = true ]; then echo "Result for tr: $tr_results" >&2; fi

    echo "${tl_exon1_results},${tl_exon2_results},${tl_results},${tr_results}"
}

##############################################################
# MAIN SCRIPT LOGIC
##############################################################

# Store input variables
REGIONS_FILE=$1
OUTPUT_NAME=$2

echo "Processing file: $REGIONS_FILE"

# Define paths
output_path="../data/datasets/gwas_pval_feature"
log_file="${output_path}/${OUTPUT_NAME}.log"
output_file="${output_path}/${OUTPUT_NAME}.csv"
if [ "$VERBOSE" = true ]; then
    echo "Log file: $log_file" >&2
    echo "Output file: $output_file" >&2
fi

# Set up a unique lock directory for rate limiting this specific run
export LOCK_DIR="/tmp/gwas_api_lock_$$"
mkdir -p "$LOCK_DIR"
# Ensure the lock directory is cleaned up on script exit, error, or interrupt
trap 'rm -rf "$LOCK_DIR"' EXIT

# Create output directory
mkdir -p "$output_path"
echo "Output will be written to: $output_path"

# Add header to output file if it doesn't exist
if [ ! -f "$output_file" ]; then
    echo "Creating new output file with header: $output_file"
    echo "tl_exon1_gwas_associations,tl_exon1_min_p_value,tl_exon1_odds_ratio,tl_exon1_min_odds_ratio,tl_exon1_max_odds_ratio,tl_exon2_gwas_associations,tl_exon2_min_p_value,tl_exon2_odds_ratio,tl_exon2_min_odds_ratio,tl_exon2_max_odds_ratio,tl_gwas_assoc,tl_min_p_val,tl_odds_ratio,tl_min_odds_ratio,tl_max_odds_ratio,tr_gwas_assoc,tr_min_p_val,tr_odds_ratio,tr_min_odds_ratio,tr_max_odds_ratio" > "$output_file"
fi

# Check for previous runs to resume processing
start_from_line=2 # Default: start after the header of the input file
if [ -f "$output_file" ]; then
    lines_in_output=$(wc -l < "$output_file")
    if [ "$lines_in_output" -gt 1 ]; then
        processed_rows=$((lines_in_output - 1))
        start_from_line=$((processed_rows + 2))
        echo "Resuming run. Found $processed_rows processed lines in $output_file."
        echo "Starting from line $start_from_line in $REGIONS_FILE."
    else
        echo "No previous output file found. Starting from the beginning."
    fi
fi

# Export the VERBOSE variable so it's available in subshells
export VERBOSE

# Export the new throttle_api_call function
export -f convert_chromosome get_unique_gwas_variants process_gwas_region_line \
        generate_cache_key read_from_cache write_to_cache min throttle_api_call

# Start parallel processing
echo "Starting parallel processing... A log will be available at: $log_file"
tail -n +${start_from_line} "$REGIONS_FILE" | tr -d '\r' | parallel --keep-order --jobs 4 --joblog "$log_file" --eta process_gwas_region_line >> "$output_file"

echo
echo "------------------------------------------------"
echo "Script completed."
echo "Results have been written to: $output_file"
echo "------------------------------------------------"