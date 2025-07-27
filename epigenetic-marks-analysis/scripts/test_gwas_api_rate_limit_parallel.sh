#!/bin/bash
# Script to empirically test the GWAS Summary Statistics API rate limits
# using multiple concurrent requests to determine maximum parallel throughput
#
# This script spawns multiple background processes to send simultaneous requests
# and measures the collective rate limit threshold.
#
# Author: Generated for parallel rate limit testing
# Date: $(date "+%Y-%m-%d")
#
###############################################################################

# Default values
VERBOSE=false
TEST_DURATION=60  # Default test duration in seconds
OUTPUT_FILE="gwas_api_parallel_rate_limit_test.log"
CONCURRENT_WORKERS=10  # Number of parallel workers
MAX_REQUESTS_PER_WORKER=1000  # Max requests per worker to prevent runaway

# Parse command line options
while getopts "vd:o:w:m:h" opt; do
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
        w )
            CONCURRENT_WORKERS=$OPTARG
            ;;
        m )
            MAX_REQUESTS_PER_WORKER=$OPTARG
            ;;
        h )
            echo "Usage: $0 [-v] [-d duration] [-o output_file] [-w workers] [-m max_requests]"
            echo "  -v: Verbose mode"
            echo "  -d: Test duration in seconds (default: 60)"
            echo "  -o: Output log file (default: gwas_api_parallel_rate_limit_test.log)"
            echo "  -w: Number of concurrent workers (default: 10)"
            echo "  -m: Max requests per worker (default: 1000)"
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

# Global variables for coordination
WORKER_PIDS=()
STOP_FLAG_FILE="/tmp/gwas_test_stop_$$"
RESULTS_DIR="/tmp/gwas_parallel_results_$$"
START_TIME=$(date +%s)

echo "=== GWAS API Parallel Rate Limit Test ==="
echo "Start time: $(date)"
echo "API URL: $API_URL"
echo "Test duration: ${TEST_DURATION} seconds"
echo "Concurrent workers: $CONCURRENT_WORKERS"
echo "Max requests per worker: $MAX_REQUESTS_PER_WORKER"
echo "Output file: $OUTPUT_FILE"
echo "=========================================="
echo

# Create results directory and initialize files
mkdir -p "$RESULTS_DIR"
rm -f "$STOP_FLAG_FILE"

# Create/clear the main output file
echo "# GWAS API Parallel Rate Limit Test Results" > "$OUTPUT_FILE"
echo "# Start time: $(date)" >> "$OUTPUT_FILE"
echo "# API URL: $API_URL" >> "$OUTPUT_FILE"
echo "# Concurrent workers: $CONCURRENT_WORKERS" >> "$OUTPUT_FILE"
echo "# Format: timestamp,worker_id,elapsed_seconds,http_code,cumulative_success,cumulative_errors" >> "$OUTPUT_FILE"

# Function to log verbose messages
log_verbose() {
    if [ "$VERBOSE" = true ]; then
        echo "[$(date '+%H:%M:%S')] $1"
    fi
}

# Worker function that runs in background
worker_function() {
    local worker_id=$1
    local worker_success=0
    local worker_errors=0
    local worker_results_file="$RESULTS_DIR/worker_${worker_id}.log"
    
    log_verbose "Worker $worker_id started"
    
    # Initialize worker results file
    echo "# Worker $worker_id results" > "$worker_results_file"
    
    local request_count=0
    while [ $request_count -lt $MAX_REQUESTS_PER_WORKER ]; do
        # Check if stop flag exists
        if [ -f "$STOP_FLAG_FILE" ]; then
            log_verbose "Worker $worker_id received stop signal"
            break
        fi
        
        # Check if test duration exceeded
        local current_time=$(date +%s)
        local elapsed=$((current_time - START_TIME))
        if [ $elapsed -ge $TEST_DURATION ]; then
            log_verbose "Worker $worker_id: test duration exceeded"
            break
        fi
        
        request_count=$((request_count + 1))
        
        # Create temporary files for this request
        local response_file=$(mktemp)
        local header_file=$(mktemp)
        
        # Send API request
        local http_code=$(curl -s -w "%{http_code}" \
                             -D "$header_file" \
                             -o "$response_file" \
                             --max-time 10 \
                             --retry 0 \
                             -H "Accept: application/json" \
                             "$API_URL" 2>/dev/null)
        
        local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
        
        # Process response
        case "$http_code" in
            200)
                worker_success=$((worker_success + 1))
                log_verbose "Worker $worker_id: Success #$worker_success"
                ;;
            409|429)
                worker_errors=$((worker_errors + 1))
                log_verbose "Worker $worker_id: Rate limit hit (HTTP $http_code)"
                echo "$timestamp,$worker_id,$elapsed,$http_code,$worker_success,$worker_errors" >> "$worker_results_file"
                # Signal other workers to stop
                touch "$STOP_FLAG_FILE"
                break
                ;;
            *)
                worker_errors=$((worker_errors + 1))
                log_verbose "Worker $worker_id: Error HTTP $http_code"
                ;;
        esac
        
        # Log result
        echo "$timestamp,$worker_id,$elapsed,$http_code,$worker_success,$worker_errors" >> "$worker_results_file"
        
        # Clean up temporary files
        rm -f "$response_file" "$header_file"
        
        # Small delay to prevent completely overwhelming the server
        sleep 0.01
    done
    
    log_verbose "Worker $worker_id finished: $worker_success successes, $worker_errors errors"
    echo "# Worker $worker_id final: $worker_success successes, $worker_errors errors" >> "$worker_results_file"
}

# Start all worker processes
echo "Starting $CONCURRENT_WORKERS parallel workers..."
for ((i=1; i<=CONCURRENT_WORKERS; i++)); do
    worker_function $i &
    WORKER_PIDS+=($!)
    log_verbose "Started worker $i (PID: ${WORKER_PIDS[$((i-1))]})"
done

# Monitor workers and collect results
echo "Workers started. Monitoring progress..."
echo "Press Ctrl+C to stop early"

# Set up signal handler for clean shutdown
cleanup() {
    echo "Stopping workers..."
    touch "$STOP_FLAG_FILE"
    for pid in "${WORKER_PIDS[@]}"; do
        kill $pid 2>/dev/null
    done
    wait
}
trap cleanup INT TERM

# Wait for workers to complete or timeout
wait

# Collect and analyze results
echo
echo "=== Collecting Results ==="

total_success=0
total_errors=0
first_rate_limit_time=""
workers_hit_limit=0

# Process each worker's results
for ((i=1; i<=CONCURRENT_WORKERS; i++)); do
    worker_file="$RESULTS_DIR/worker_${i}.log"
    if [ -f "$worker_file" ]; then
        # Extract worker summary
        worker_success=$(grep "^# Worker $i final:" "$worker_file" | cut -d' ' -f5)
        worker_errors=$(grep "^# Worker $i final:" "$worker_file" | cut -d' ' -f7)
        
        # Add to totals
        total_success=$((total_success + ${worker_success:-0}))
        total_errors=$((total_errors + ${worker_errors:-0}))
        
        # Check if this worker hit rate limit
        if grep -q ",409," "$worker_file" || grep -q ",429," "$worker_file"; then
            workers_hit_limit=$((workers_hit_limit + 1))
            if [ -z "$first_rate_limit_time" ]; then
                first_rate_limit_time=$(grep -E ",409,|,429," "$worker_file" | head -1 | cut -d',' -f1)
            fi
        fi
        
        # Append worker results to main file
        echo "# === Worker $i Results ===" >> "$OUTPUT_FILE"
        grep -v "^#" "$worker_file" >> "$OUTPUT_FILE"
    fi
done

# Calculate final statistics
end_time=$(date +%s)
total_elapsed=$((end_time - START_TIME))
total_requests=$((total_success + total_errors))

echo
echo "=== Final Results ==="
echo "End time: $(date)"
echo "Total duration: $total_elapsed seconds"
echo "Total successful requests: $total_success"
echo "Total error responses: $total_errors"
echo "Total requests sent: $total_requests"
echo "Workers that hit rate limits: $workers_hit_limit"

if [ $total_elapsed -gt 0 ]; then
    avg_requests_per_sec=$(echo "scale=2; $total_success / $total_elapsed" | bc -l 2>/dev/null)
    avg_total_per_sec=$(echo "scale=2; $total_requests / $total_elapsed" | bc -l 2>/dev/null)
    echo "Average successful requests per second: ${avg_requests_per_sec:-N/A}"
    echo "Average total requests per second: ${avg_total_per_sec:-N/A}"
fi

if [ -n "$first_rate_limit_time" ]; then
    echo "First rate limit hit at: $first_rate_limit_time"
fi

echo "Detailed log saved to: $OUTPUT_FILE"
echo "======================="

# Add summary to log file
echo "" >> "$OUTPUT_FILE"
echo "# === FINAL SUMMARY ===" >> "$OUTPUT_FILE"
echo "# Total duration: $total_elapsed seconds" >> "$OUTPUT_FILE"
echo "# Concurrent workers: $CONCURRENT_WORKERS" >> "$OUTPUT_FILE"
echo "# Total successful requests: $total_success" >> "$OUTPUT_FILE"
echo "# Total error responses: $total_errors" >> "$OUTPUT_FILE"
echo "# Total requests sent: $total_requests" >> "$OUTPUT_FILE"
echo "# Workers that hit rate limits: $workers_hit_limit" >> "$OUTPUT_FILE"
echo "# Average successful requests per second: ${avg_requests_per_sec:-N/A}" >> "$OUTPUT_FILE"
echo "# Average total requests per second: ${avg_total_per_sec:-N/A}" >> "$OUTPUT_FILE"
if [ -n "$first_rate_limit_time" ]; then
    echo "# First rate limit hit at: $first_rate_limit_time" >> "$OUTPUT_FILE"
fi

# Clean up temporary files
rm -rf "$RESULTS_DIR"
rm -f "$STOP_FLAG_FILE"
