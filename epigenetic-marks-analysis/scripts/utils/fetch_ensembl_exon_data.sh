#!/bin/bash

# Script Name: fetch_ensembl_exon_data.sh
# Author: Estefania Rojas
# Description: Retrieves exon 1 and 2 coords for lncRNAs from ENSEMBL API using ENSEMBL transcript ids
#
# Dependencies: 
# - jq: https://jqlang.github.io/jq/
# - GNU parallel (optional, for concurrent execution)
#
# Usage: ./fetch_ensembl_exon_data.sh regions_to_extract.csv outname.csv

set -eo pipefail

# Configuration
readonly MAX_RETRIES=3
readonly RETRY_DELAY=5
readonly RATE_LIMIT_DELAY=0.5  # Delay between API calls in seconds
readonly MAX_PARALLEL_JOBS=5
readonly TEMP_DIR="/tmp/ensembl_fetch_$$"
readonly LOG_FILE="../logs/fetch_ensembl_$(date +%Y%m%d_%H%M%S).log"

# Ensure cleanup on script exit
cleanup() {
    local exit_code=$?
    rm -rf "${TEMP_DIR}"
    exit $exit_code
}

trap cleanup EXIT
trap 'trap - EXIT; cleanup' INT TERM

# Logging function
log() {
    local level=$1
    shift
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [${level}] $*" | tee -a "${LOG_FILE}"
}

# Function to make API requests with retries and rate limiting
make_api_request() {
    local url=$1
    local retry_count=0
    local response=""

    while [ $retry_count -lt $MAX_RETRIES ]; do
        response=$(curl -s --max-time 30 --retry 2 "${url}")
        if [ $? -eq 0 ] && [ -n "$response" ] && ! echo "$response" | grep -q "^<html>"; then
            echo "$response"
            return 0
        fi
        retry_count=$((retry_count + 1))
        log "WARN" "API request failed, attempt $retry_count of $MAX_RETRIES. URL: $url"
        sleep $RETRY_DELAY
    done
    log "ERROR" "API request failed after $MAX_RETRIES attempts. URL: $url"
    echo "NA"
    return 1
}

# Process response function with improved error handling
process_response() {
    local response="$1"
    local result
    
    [ -z "$response" ] && { echo "NA"; return 0; }
    echo "$response" | grep -q "^<html>" && { echo "NA"; return 0; }
    
    result=$(echo "$response" | jq -r '
        try (
            if (.results | length > 0) then
                .results[0].rnacentral_id
            else
                "NA"
            end
        ) catch "NA"
    ' 2>/dev/null) || { echo "NA"; return 0; }
    
    echo "$result"
}

# Progress tracking function
show_progress() {
    local current=$1
    local total=$2
    local width=50
    local percentage=$((current * 100 / total))
    local completed=$((width * current / total))
    local remaining=$((width - completed))
    
    printf "\rProgress: ["
    printf "%${completed}s" | tr ' ' '#'
    printf "%${remaining}s" | tr ' ' '-'
    printf "] %d%% (%d/%d)" "$percentage" "$current" "$total"
}

# Input validation
validate_input() {
    local regions_file=$1
    local output_name=$2

    [ $# -ne 2 ] && { 
        log "ERROR" "Usage: $0 regions_to_extract.csv outname.csv"
        exit 1
    }

    [ ! -f "$regions_file" ] && {
        log "ERROR" "Input file not found: $regions_file"
        exit 1
    }

    [ ! -r "$regions_file" ] && {
        log "ERROR" "Input file not readable: $regions_file"
        exit 1
    }
}

# Process single transcript
process_transcript() {
    local trans_id=$1
    local api_version=$2  # 'current' or 'archive'
    local base_url

    [ "$api_version" = "archive" ] && base_url="https://may2024.rest.ensembl.org" || base_url="https://rest.ensembl.org"
    
    local ensembl_response
    ensembl_response=$(make_api_request "${base_url}/lookup/id/$(echo "${trans_id%.*}")?expand=1;content-type=application/json")
    
    # Process response and extract fields
    local name=$(echo "$ensembl_response" | jq -r '.display_name // "NA"')
    local fields=$(echo "$ensembl_response" | jq -r '[.seq_region_name // "NA", .start // "NA", .end // "NA"] | @csv')
    local exon1_fields=$(echo "$ensembl_response" | jq -r '[ .Exon[0].start // "NA", .Exon[0].end // "NA"] | @csv')
    local exon2_fields=$(echo "$ensembl_response" | jq -r '[ .Exon[1].start // "NA", .Exon[1].end // "NA"] | @csv')
    
    # Rate limiting
    sleep $RATE_LIMIT_DELAY
    
    echo "${name}|${fields}|${exon1_fields}|${exon2_fields}"
}

# Main execution
main() {
    validate_input "$@"
    
    local regions_file=$1
    local output_name=$2
    local output_path="../data/datasets/ensembl_feature"
    local processed=0
    
    # Create required directories
    mkdir -p "$output_path" "../logs"
    mkdir -p "${TEMP_DIR}"
    
    log "INFO" "Starting processing of $regions_file"
    
    # Initialize or resume from existing output
    if [ ! -f "$output_path/$output_name.csv" ]; then
        log "INFO" "Creating new output file"
        echo "ensembl_trans_id,rnacentral_id,trans_name,tl_chromosome,tl_start,tl_end,tl_exon1_start,tl_exon1_end,tl_exon2_start,tl_exon2_end,tr_chromosome,tr_start,tr_end,tr_exon1_start,tr_exon1_end,tr_exon2_start,tr_exon2_end" > "$output_path/$output_name.csv"
        processed=1
    else
        processed=$(wc -l < "$output_path/$output_name.csv")
        log "INFO" "Resuming from line $processed"
    fi
    
    total_lines=$(wc -l < "$regions_file")
    
    # Process each line
    tail -n +$processed "$regions_file" | while IFS='|' read -r t_id1 t_id2; do
        [ -z "$t_id1" ] && continue
        
        # Process first transcript
        result1=$(process_transcript "$t_id1" "current")
        
        # Process second transcript if different
        if [ "$t_id1" = "$t_id2" ] || [ "$t_id2" = "NA" ]; then
            result2=$result1
        else
            result2=$(process_transcript "$t_id2" "current")
        fi
        
        # Combine results and append to output
        echo "${t_id1}|${t_id2},${result1},${result2}" >> "$output_path/$output_name.csv"
        
        processed=$((processed + 1))
        show_progress $processed $total_lines
    done
    
    log "INFO" "Processing completed. Output saved to $output_path/$output_name.csv"
}

main "$@"
