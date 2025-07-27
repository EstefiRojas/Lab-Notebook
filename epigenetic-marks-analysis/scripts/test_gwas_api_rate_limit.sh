#!/bin/bash
# Script to empirically test the GWAS Summary Statistics API rate limits
# This script sends requests until it receives a 409 (rate limit) error
# and measures the number of successful requests per time period.
#
# Author: Generated for rate limit testing
# Date: $(date "+%Y-%m-%d")
#
###############################################################################

# Default values
VERBOSE=false
TEST_DURATION=60  # Default test duration in seconds
OUTPUT_FILE="gwas_api_rate_limit_test.log"

# Parse command line options
while getopts "vd:o:h" opt; do
    case ${opt} in
        v )
            VERBOSE=true
            ;;
        d )
            TEST_DURATION=$OPTARG
            ;;
        o )
            OUTPUT_FILE=$OPTARG
            ;;
        h )
            echo "Usage: $0 [-v] [-d duration] [-o output_file]"
            echo "  -v: Verbose mode"
            echo "  -d: Test duration in seconds (default: 60)"
            echo "  -o: Output log file (default: gwas_api_rate_limit_test.log)"
            exit 0
            ;;
        \? )
            echo "Invalid option: -$OPTARG" 1>&2
            exit 1
            ;;
    esac
done

# GWAS API endpoint - using a simple query that should return minimal data
API_URL="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/2/associations?bp_lower=54661012&bp_upper=54661375&p_upper=1e-8&size=10"

# Initialize counters
SUCCESS_COUNT=0
ERROR_COUNT=0
RATE_LIMIT_HIT=false
START_TIME=$(date +%s)
LAST_SECOND_COUNT=0
REQUESTS_PER_SECOND=()

echo "=== GWAS API Rate Limit Test ==="
echo "Start time: $(date)"
echo "API URL: $API_URL"
echo "Test duration: ${TEST_DURATION} seconds"
echo "Output file: $OUTPUT_FILE"
echo "=================================="
echo

# Create/clear the output file
echo "# GWAS API Rate Limit Test Results" > "$OUTPUT_FILE"
echo "# Start time: $(date)" >> "$OUTPUT_FILE"
echo "# API URL: $API_URL" >> "$OUTPUT_FILE"
echo "# Format: timestamp,elapsed_seconds,http_code,success_count,requests_this_second" >> "$OUTPUT_FILE"

# Function to log verbose messages
log_verbose() {
    if [ "$VERBOSE" = true ]; then
        echo "[$(date '+%H:%M:%S')] $1"
    fi
}

# Main testing loop
while true; do
    CURRENT_TIME=$(date +%s)
    ELAPSED=$((CURRENT_TIME - START_TIME))
    
    # Check if test duration exceeded
    if [ $ELAPSED -ge $TEST_DURATION ]; then
        echo "Test duration of ${TEST_DURATION} seconds reached."
        break
    fi
    
    # Track requests per second
    if [ $ELAPSED -ne $LAST_SECOND_COUNT ]; then
        REQUESTS_THIS_SECOND=0
        LAST_SECOND_COUNT=$ELAPSED
    fi
    
    log_verbose "Sending request #$((SUCCESS_COUNT + ERROR_COUNT + 1))"
    
    # Create temporary files for response and headers
    RESPONSE_FILE=$(mktemp)
    HEADER_FILE=$(mktemp)
    
    # Send API request with timeout and capture both response and headers
    HTTP_CODE=$(curl -s -w "%{http_code}" \
                     -D "$HEADER_FILE" \
                     -o "$RESPONSE_FILE" \
                     --max-time 10 \
                     --retry 0 \
                     -H "Accept: application/json" \
                     "$API_URL")
    
    CURRENT_TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
    REQUESTS_THIS_SECOND=$((REQUESTS_THIS_SECOND + 1))
    
    # Process the response based on HTTP code
    case "$HTTP_CODE" in
        200)
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            log_verbose "✓ Success #$SUCCESS_COUNT (HTTP 200)"
            echo "$CURRENT_TIMESTAMP,$ELAPSED,$HTTP_CODE,$SUCCESS_COUNT,$REQUESTS_THIS_SECOND" >> "$OUTPUT_FILE"
            ;;
        429)
            echo "🚫 Rate limit reached! (HTTP 409)"
            echo "Total successful requests before rate limit: $SUCCESS_COUNT"
            echo "Time elapsed: $ELAPSED seconds"
            echo "Average requests per second: $(echo "scale=2; $SUCCESS_COUNT / $ELAPSED" | bc -l 2>/dev/null || echo "N/A")"
            echo "$CURRENT_TIMESTAMP,$ELAPSED,$HTTP_CODE,$SUCCESS_COUNT,$REQUESTS_THIS_SECOND" >> "$OUTPUT_FILE"
            RATE_LIMIT_HIT=true
            break
            ;;
        000)
            echo "🔌 Connection timeout or network error"
            ERROR_COUNT=$((ERROR_COUNT + 1))
            echo "$CURRENT_TIMESTAMP,$ELAPSED,$HTTP_CODE,$SUCCESS_COUNT,$REQUESTS_THIS_SECOND" >> "$OUTPUT_FILE"
            sleep 0.5
            ;;
        *)
            echo "❌ Unexpected HTTP code: $HTTP_CODE"
            ERROR_COUNT=$((ERROR_COUNT + 1))
            echo "$CURRENT_TIMESTAMP,$ELAPSED,$HTTP_CODE,$SUCCESS_COUNT,$REQUESTS_THIS_SECOND" >> "$OUTPUT_FILE"
            
            if [ "$VERBOSE" = true ] && [ -f "$RESPONSE_FILE" ]; then
                echo "Response body:"
                cat "$RESPONSE_FILE"
                echo
            fi
            ;;
    esac
    
    # Clean up temporary files
    rm -f "$RESPONSE_FILE" "$HEADER_FILE"
    
    # Small delay to prevent overwhelming the server immediately
    # This can be adjusted based on testing needs
    sleep 0.05
done

# Final summary
END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - START_TIME))

echo
echo "=== Test Summary ==="
echo "End time: $(date)"
echo "Total duration: $TOTAL_ELAPSED seconds"
echo "Successful requests: $SUCCESS_COUNT"
echo "Error responses: $ERROR_COUNT"
echo "Rate limit hit: $RATE_LIMIT_HIT"

if [ $TOTAL_ELAPSED -gt 0 ]; then
    AVG_REQUESTS_PER_SEC=$(echo "scale=2; $SUCCESS_COUNT / $TOTAL_ELAPSED" | bc -l 2>/dev/null)
    echo "Average requests per second: ${AVG_REQUESTS_PER_SEC:-N/A}"
fi

echo "Detailed log saved to: $OUTPUT_FILE"
echo "===================="

# Add summary to log file
echo "" >> "$OUTPUT_FILE"
echo "# === SUMMARY ===" >> "$OUTPUT_FILE"
echo "# Total duration: $TOTAL_ELAPSED seconds" >> "$OUTPUT_FILE"
echo "# Successful requests: $SUCCESS_COUNT" >> "$OUTPUT_FILE"
echo "# Error responses: $ERROR_COUNT" >> "$OUTPUT_FILE"
echo "# Rate limit hit: $RATE_LIMIT_HIT" >> "$OUTPUT_FILE"
echo "# Average requests per second: ${AVG_REQUESTS_PER_SEC:-N/A}" >> "$OUTPUT_FILE"
