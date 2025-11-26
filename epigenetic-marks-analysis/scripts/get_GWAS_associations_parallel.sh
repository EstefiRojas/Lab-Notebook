#!/bin/bash
#
# GWAS Association Data Parallel Retrieval Script
# ===============================================
#
# DESCRIPTION:
#   This script retrieves GWAS (Genome-Wide Association Studies) association data
#   in parallel from the EBI GWAS Summary Statistics API. It processes genomic regions
#   from CSV input files and extracts statistical measures including p-values,
#   odds ratios, and beta coefficients.
#
# FEATURES:
#   - Parallel processing with GNU parallel for improved performance
#   - File-based JSON caching to avoid redundant API calls
#   - Rate limiting to respect API constraints (max 10 calls/second)
#   - Automatic retry logic with exponential backoff
#   - Resume capability for interrupted runs
#   - Cross-platform compatibility (macOS/Linux)
#   - Comprehensive error handling and logging
#
# USAGE:
#   ./get_GWAS_associations_parallel.sh [-v] <regions_to_extract.csv> <output_name>
#
# ARGUMENTS:
#   regions_to_extract.csv  Input CSV file containing genomic regions to query
#   output_name            Base name for output files (without .csv extension)
#
# OPTIONS:
#   -v    Enable verbose mode for detailed logging output
#
# INPUT FORMAT:
#   CSV file with columns including chromosome positions for:
#   - TL exon1 regions (columns 9-11: chrm, start, end)
#   - TL exon2 regions (columns 12-14: chrm, start, end)  
#   - TL regions (columns 15-17: chrm, start, end)
#   - TR regions (columns 18-20: chrm, start, end)
#
# OUTPUT:
#   - CSV file with GWAS association statistics for each region
#   - Log file with processing details and job information
#   - Cached JSON responses in ../data/cache/gwas_responses/
#
# DEPENDENCIES:
#   - curl:     HTTP client for API requests
#   - jq:       JSON processor (https://jqlang.github.io/jq/)
#   - parallel: GNU parallel for job parallelization
#   - sed:      Stream editor for text processing
#
# API ENDPOINT:
#   https://www.ebi.ac.uk/gwas/summary-statistics/api/
#
# AUTHOR: Estefania Rojas
# DATE: 31 Jul 2025
#
###############################################################################

# Get the date of execution
DATE=$(date "+%A %B %d, %Y")
echo "Execution date: ${DATE}"

# Default values for options
VERBOSE=false
OFFLINE_MODE=false

# Parse options for verbose mode and offline mode
while getopts "v-:" opt; do
  case ${opt} in
    v )
      VERBOSE=true
      ;;
    - )
      case "${OPTARG}" in
        offline-mode)
          OFFLINE_MODE=true
          ;;
        *)
          echo "Invalid option: --$OPTARG" 1>&2
          exit 1
          ;;
      esac
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

# Update the usage message
if [ $# -ne 2 ]; then
    echo "Usage: $0 [-v] [--offline-mode] <regions_to_extract.csv> <outname.csv>"
    echo "  -v              Enable verbose mode for detailed logging output"
    echo "  --offline-mode  Only use cached data, never query the API"
    exit 1
fi

# Export the OFFLINE_MODE variable so it's available in subshells
export OFFLINE_MODE

##############################################################
# HELPER FUNCTIONS
##############################################################

# Function: check_dependencies
# Purpose: Verify that all required command-line tools are available
# Returns: Exits with code 1 if any dependencies are missing
# Dependencies checked: curl, jq, parallel, sed
check_dependencies() {
    local missing=0
    for cmd in curl jq parallel sed; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "Error: Required command '$cmd' is not installed or not in your PATH." >&2
            missing=1
        fi
    done
    if [ "$missing" -eq 1 ]; then
        exit 1
    fi
    echo "All dependencies (curl, jq, parallel, sed) are present."
}

# Function: convert_chromosome
# Purpose: Normalize chromosome names to numerical format for API compatibility
# Arguments: $1 - chromosome name (e.g., "chr1", "X", "23")
# Returns: Numerical chromosome identifier (X->23, Y->24)
# Note: Removes 'chr' prefix and handles sex chromosomes
convert_chromosome() {
    local chrm=$1
    # Remove 'chr' prefix (case-insensitive) and convert to uppercase
    chrm=$(echo "$chrm" | sed -E 's/^chr//i' | tr '[:lower:]' '[:upper:]')
    case "$chrm" in
        "X") echo "23" ;;  # X chromosome -> 23 for API compatibility
        "Y") echo "24" ;;  # Y chromosome -> 24 for API compatibility
        *) echo "$chrm" ;; # Return as-is for autosomes
    esac
}

# Function: min
# Purpose: Mathematical utility to return the smaller of two integers
# Arguments: $1, $2 - two integers to compare
# Returns: The smaller value
min() {
    echo $(( $1 < $2 ? $1 : $2 ))
}

# Function: generate_cache_key
# Purpose: Create a unique cache key from genomic coordinates
# Arguments: $1 - chromosome, $2 - start position, $3 - end position
# Returns: Cache key string in format "chrm_start_end"
# Usage: Used to name cache files for API responses
generate_cache_key() {
    echo "$1_$2_$3"
}

# A portable and race-condition-safe rate-limiting function.
throttle_api_call() {
    local max_calls_per_second=2
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
                file_epoch=$(stat $stat_flags "$token_file" 2>/dev/null)
                
                if [ -n "$file_epoch" ]; then
                    if (( now - file_epoch > 1 )); then
                        rm -f "$token_file"
                    fi
                fi
            fi
        done

        local current_calls
        current_calls=$(find "$lock_dir" -type f | wc -l)
        
        if [ "$current_calls" -lt "$max_calls_per_second" ]; then
            touch "${lock_dir}/$(date +%s.%N)"
            break
        else
            sleep 0.1
        fi
    done
}

# Function to process GWAS JSON data and compute summary statistics
process_gwas_json() {
    local gwas_json="$1"
    local job_slot="$2"
    local data_file

    data_file=$(mktemp)
    # Ensure the temp file is removed on function exit
    trap 'rm -f "$data_file"' RETURN

    # Convert JSON array of associations to a TSV file for processing
    echo "$gwas_json" | jq -r '.[] | select(try (.p_value|tonumber) catch false) | select(.p_value != -99) | [.variant_id, .p_value, (.odds_ratio // "null"), (.beta // "null")] | @tsv' > "$data_file"

    # If no valid data, return NA for all fields
    if [ ! -s "$data_file" ]; then
        if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] No valid association data found in JSON. Returning NA." >&2; fi
        echo "NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA"
        return
    fi
    
    # Calculate number of unique variants
    local unique_count
    unique_count=$(cut -f1 "$data_file" | sort -u | wc -l | tr -d ' ')

    # Find the minimum p-value and its corresponding odds ratio and beta
    local min_p_line
    min_p_line=$(cut -f2,3,4 "$data_file" | sort -t$'\t' -k1,1g | head -n1)
    
    local min_p_value="NA"
    local odds_ratio_for_min_p="NA"
    local beta_for_min_p="NA"

    if [[ -n "$min_p_line" ]]; then
        IFS=$'\t' read -r min_p_value or_val beta_val <<< "$min_p_line"
        if [[ "$or_val" == "null" || -z "$or_val" ]]; then
            odds_ratio_for_min_p="NA"
        else
            odds_ratio_for_min_p="$or_val"
        fi
        if [[ "$beta_val" == "null" || -z "$beta_val" ]]; then
            beta_for_min_p="NA"
        else
            beta_for_min_p="$beta_val"
        fi
    fi

    # Filter for valid odds ratio pairs to find min/max
    local valid_or_pairs
    valid_or_pairs=$(cut -f2,3 "$data_file" | grep -E $'\t[0-9.eE+-]+$')

    local odds_ratio_min="NA"
    local pval_for_min_or="NA"
    local odds_ratio_max="NA"
    local pval_for_max_or="NA"

    if [ -n "$valid_or_pairs" ]; then
        local sorted_or_pairs
        sorted_or_pairs=$(echo "$valid_or_pairs" | sort -t$'\t' -k2,2g)
        
        local min_or_line
        min_or_line=$(echo "$sorted_or_pairs" | head -n1)
        IFS=$'\t' read -r pval_for_min_or odds_ratio_min <<< "$min_or_line"

        local max_or_line
        max_or_line=$(echo "$sorted_or_pairs" | tail -n1)
        IFS=$'\t' read -r pval_for_max_or odds_ratio_max <<< "$max_or_line"
    fi
    
    # Filter for valid beta pairs to find min/max
    local valid_beta_pairs
    valid_beta_pairs=$(cut -f2,4 "$data_file" | grep -E $'\t[0-9.eE+-]+$')

    local beta_min="NA"
    local pval_for_min_beta="NA"
    local beta_max="NA"
    local pval_for_max_beta="NA"

    if [ -n "$valid_beta_pairs" ]; then
        local sorted_beta_pairs
        sorted_beta_pairs=$(echo "$valid_beta_pairs" | sort -t$'\t' -k2,2g)
        
        local min_beta_line
        min_beta_line=$(echo "$sorted_beta_pairs" | head -n1)
        IFS=$'\t' read -r pval_for_min_beta beta_min <<< "$min_beta_line"

        local max_beta_line
        max_beta_line=$(echo "$sorted_beta_pairs" | tail -n1)
        IFS=$'\t' read -r pval_for_max_beta beta_max <<< "$max_beta_line"
    fi

    # Assemble the final result string
    local final_result="${unique_count},${min_p_value},${odds_ratio_for_min_p},${beta_for_min_p},${odds_ratio_min},${pval_for_min_or},${odds_ratio_max},${pval_for_max_or},${beta_min},${pval_for_min_beta},${beta_max},${pval_for_max_beta}"
    echo "$final_result"
}


# Main function to get GWAS data. It uses a file-based cache for full JSON
# responses to avoid re-querying the API. It computes number of unique
# variants, minimum p-value with its corresponding odds-ratio, max odds-ratio,
# and min odds-ratio with their corresponding p-values.
get_unique_gwas_variants() {
    local chrm
    local start_pos
    local end_pos
    local job_slot="$4" # Capture the job slot number
    chrm=$(convert_chromosome "$1")
    start_pos=$2
    end_pos=$3

    # Return immediately if input is invalid
    if [[ -z "$chrm" || "$chrm" == "NULL" || -z "$start_pos" || "$start_pos" == "0" ]]; then
        if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Invalid input for get_unique_gwas_variants: chrm='$chrm', start_pos='$start_pos', end_pos='$end_pos'. Skipping." >&2; fi
        echo "NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA"
        return
    fi

    # Define cache path and key
    local cache_dir="../data/cache/gwas_responses"
    mkdir -p "$cache_dir"
    local cache_key
    cache_key=$(generate_cache_key "$chrm" "$start_pos" "$end_pos")
    local cache_file="${cache_dir}/${cache_key}.json"

    # Try to get results from cache
    if [[ -f "$cache_file" ]]; then
        if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Cache hit for key '$cache_key'. Processing cached JSON." >&2; fi
        local all_responses
        all_responses=$(cat "$cache_file")
        # Extract associations from the cached full responses
        local all_associations
        all_associations=$(echo "$all_responses" | jq 'map(._embedded.associations // []) | add')
        process_gwas_json "$all_associations" "$job_slot"
        return
    fi
    
    # Check if offline mode is enabled
    if [ "$OFFLINE_MODE" = true ]; then
        if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Offline mode enabled. Cache miss for key '$cache_key'. Returning NA values." >&2; fi
        echo "NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA"
        return
    fi

    if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Cache miss for key '$cache_key'. Querying API." >&2; fi

    # API query logic with pagination handling
    local max_page_size=500
    local current_page_size=$max_page_size
    local current_url="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/${chrm}/associations?bp_lower=${start_pos}&bp_upper=${end_pos}&p_upper=5e-8&size=${current_page_size}"
    local all_responses="[]"
    local success=false

    while true; do
        if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Calling API with URL: $current_url" >&2; fi
        throttle_api_call

        local gwas_response
        local header_file
        header_file=$(mktemp)
        gwas_response=$(curl -s -L -D "$header_file" --retry 2 --max-time 60 "${current_url}" -H "Accept:application/json")
        local http_code
        http_code=$(grep -i "^HTTP/" "${header_file}" | tail -n1 | cut -d' ' -f2)
        rm "$header_file"
        if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] API response code: $http_code" >&2; fi

        if [ "$http_code" = "200" ]; then
            if ! echo "$gwas_response" | jq empty 2>/dev/null; then
                if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Invalid JSON response received. Breaking loop." >&2; fi
                break
            fi

            # Add the full JSON response to the list of responses
            all_responses=$(jq -s '.[0] + [.[1]]' <(echo "$all_responses") <(echo "$gwas_response"))
            
            if [ "$current_page_size" -lt "$max_page_size" ]; then
                current_page_size=$(min $((current_page_size + 10)) $max_page_size)
            fi

            local next_url
            next_url=$(echo "$gwas_response" | jq -r '._links.next.href // empty')
            if [[ -z "$next_url" ]]; then
                if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] No next page. Exiting API loop." >&2; fi
                success=true
                break
            fi
            current_url=$(echo "$next_url" | sed "s/size=[0-9]*/size=${current_page_size}/")
        
        elif [ "$http_code" = "429" ]; then
            if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] API rate limit hit (429). Waiting 60 seconds." >&2; fi
            sleep 60
            continue

        elif echo "$gwas_response" | grep -q "Study GCST.*does not exist"; then
            if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Specific API error detected. Reducing page size and retrying." >&2; fi
            local retry_offset=1
            current_page_size=1 # Drastically reduce page size
            local current_start
            current_start=$(echo "$current_url" | grep -o 'start=[0-9]*' | cut -d'=' -f2)
            [ -z "$current_start" ] && current_start=0
            local new_start=$((current_start + retry_offset))
            
            if echo "$current_url" | grep -q 'start='; then
                 current_url=$(echo "$current_url" | sed -e "s/size=[0-9]*/size=${current_page_size}/" -e "s/start=[0-9]*/start=${new_start}/")
            else
                 current_url="${current_url}&start=${new_start}"
                 current_url=$(echo "$current_url" | sed "s/size=[0-9]*/size=${current_page_size}/")
            fi
            continue
        else
            if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Unhandled API error ($http_code). Breaking loop." >&2; fi
            break
        fi
    done

    local all_associations
    all_associations=$(echo "$all_responses" | jq 'map(._embedded.associations // []) | add')

    # Write collected JSON to cache if the fetch was successful
    if [ "$success" = true ]; then
        echo "$all_responses" > "$cache_file"
        if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Wrote response for key '$cache_key' to cache." >&2; fi
    else
        if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Did not write to cache for key '$cache_key' due to incomplete fetch." >&2; fi
    fi

    # Process the collected JSON to get the final result
    local final_result
    final_result=$(process_gwas_json "$all_associations" "$job_slot")
    if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Final result for key '$cache_key': $final_result" >&2; fi
    
    echo "$final_result"
}

# Wrapper function to process a single line from the input CSV in parallel
process_gwas_region_line() {
    local input_with_num="$1"
    local job_slot="$2" # Capture the job slot number

    # Robustly split the input (e.g., "2\tcontent...") into line number and content
    # This is safer than using 'read' with IFS in some shell environments.
    local line_num="${input_with_num%%$'\t'*}"
    local line="${input_with_num#*$'\t'}"

    if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Processing line number ${line_num}: $line" >&2; fi
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
    # Pass the job slot number to the main function
    tl_exon1_results=$(get_unique_gwas_variants "$tl_exon1_chrm" "$tl_exon1_start" "$tl_exon1_end" "$job_slot")
    if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Result for tl_exon1: $tl_exon1_results" >&2; fi
    tl_exon2_results=$(get_unique_gwas_variants "$tl_exon2_chrm" "$tl_exon2_start" "$tl_exon2_end" "$job_slot")
    if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Result for tl_exon2: $tl_exon2_results" >&2; fi
    tl_results=$(get_unique_gwas_variants "$tl_chrm" "$tl_start" "$tl_end" "$job_slot")
    if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Result for tl: $tl_results" >&2; fi
    tr_results=$(get_unique_gwas_variants "$tr_chrm" "$tr_start" "$tr_end" "$job_slot")
    if [ "$VERBOSE" = true ]; then echo "[Job ${job_slot}] Result for tr: $tr_results" >&2; fi

    echo "${line_num},${tl_exon1_results},${tl_exon2_results},${tl_results},${tr_results}"
}

##############################################################
# MAIN SCRIPT LOGIC
##############################################################

# Run dependency check at the start
check_dependencies

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
    header="tl_exon1_gwas_associations,tl_exon1_min_p_value,tl_exon1_odds_ratio_for_min_p,tl_exon1_beta_for_min_p,tl_exon1_min_odds_ratio,tl_exon1_pval_for_min_or,tl_exon1_max_odds_ratio,tl_exon1_pval_for_max_or,tl_exon1_min_beta,tl_exon1_pval_for_min_beta,tl_exon1_max_beta,tl_exon1_pval_for_max_beta,"
    header+="tl_exon2_gwas_associations,tl_exon2_min_p_value,tl_exon2_odds_ratio_for_min_p,tl_exon2_beta_for_min_p,tl_exon2_min_odds_ratio,tl_exon2_pval_for_min_or,tl_exon2_max_odds_ratio,tl_exon2_pval_for_max_or,tl_exon2_min_beta,tl_exon2_pval_for_min_beta,tl_exon2_max_beta,tl_exon2_pval_for_max_beta,"
    header+="tl_gwas_assoc,tl_min_p_val,tl_odds_ratio_for_min_p,tl_beta_for_min_p,tl_min_odds_ratio,tl_pval_for_min_or,tl_max_odds_ratio,tl_pval_for_max_or,tl_min_beta,tl_pval_for_min_beta,tl_max_beta,tl_pval_for_max_beta,"
    header+="tr_gwas_assoc,tr_min_p_val,tr_odds_ratio_for_min_p,tr_beta_for_min_p,tr_min_odds_ratio,tr_pval_for_min_or,tr_max_odds_ratio,tr_pval_for_max_or,tr_min_beta,tr_pval_for_min_beta,tr_max_beta,tr_pval_for_max_beta"
    echo "$header" > "$output_file"
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

# Export functions
export -f convert_chromosome get_unique_gwas_variants process_gwas_region_line \
        generate_cache_key min throttle_api_call process_gwas_json

# Define a temporary file for unordered results
temp_output_file="${output_path}/${OUTPUT_NAME}.tmp"
if [ -f "$temp_output_file" ]; then
    echo "Removing old temporary file: $temp_output_file"
    rm -f "$temp_output_file"
fi

# Start parallel processing without keeping order
echo "Starting parallel processing... A log will be available at: $log_file"
tail -n +${start_from_line} "$REGIONS_FILE" | tr -d '\r' | tr -dc '[:print:]\n' | sed 's/"//g' | \
    nl -v "${start_from_line}" -w1 -s $'\t' | \
    parallel --jobs 1 --joblog "$log_file" --eta "/bin/bash -c 'process_gwas_region_line \"\$1\" \"\$2\"' _ {} {%}" >> "$temp_output_file"

# Sort the temporary file by line number, remove the line number, and append to the final output file
if [ -f "$temp_output_file" ] && [ -s "$temp_output_file" ]; then
    echo "Sorting results and finalizing output..."
    sort -t, -k1,1n "$temp_output_file" | cut -d, -f2- >> "$output_file"
    rm "$temp_output_file"
else
    echo "No new data was processed."
fi

echo
echo "------------------------------------------------"
echo "Script completed."
echo "Results have been written to: $output_file"
echo "------------------------------------------------"
