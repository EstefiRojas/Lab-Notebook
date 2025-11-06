#!/opt/local/bin/bash

# Optimized script using associative arrays for better performance
# Usage: ./add_match_column_optimized.sh input.csv input.psl output.csv

# Check bash version (associative arrays require bash 4+)
if [ "${BASH_VERSION%%.*}" -lt 4 ]; then
    echo "Error: This script requires bash version 4 or higher"
    echo "Current version: $BASH_VERSION"
    exit 1
fi

# Check if correct number of arguments provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_csv> <input_psl> <output_csv>"
    echo "Example: $0 file1.csv blat_output.psl file1_with_match.csv"
    exit 1
fi

INPUT_CSV="$1"
INPUT_PSL="$2"
OUTPUT_CSV="$3"

# Check if input files exist
if [ ! -f "$INPUT_CSV" ]; then
    echo "Error: Input CSV file '$INPUT_CSV' not found!"
    exit 1
fi

if [ ! -f "$INPUT_PSL" ]; then
    echo "Error: Input PSL file '$INPUT_PSL' not found!"
    exit 1
fi

# Declare associative array to store Target_Gene_IDs from PSL
declare -A psl_targets

echo "Processing PSL file to extract Target_Gene_IDs..."

# Function to extract Target_Gene_ID from query name
extract_target_id() {
    local query="$1"
    local target_id=""
    
    # Handle different possible formats
    # Format 1: gRNA_ID_Target_Gene_ID_lncRNA (e.g., gL_000047_Hum_XLOC_000060_lncRNA)
    # Format 2: May have additional underscores in Target_Gene_ID
    
    # Remove leading gRNA_ID (starts with gL_)
    local temp="${query#gL_*_}"
    
    # Remove trailing _lncRNA or other suffixes
    target_id="${temp%_lncRNA*}"
    target_id="${temp%_lncRNA}"
    
    # If we still have a valid ID, use it
    if [[ "$target_id" =~ ^Hum_XLOC_[0-9]+ ]]; then
        echo "$target_id"
    fi
}

# Process PSL file
{
    # Skip headers (lines starting with 'psLayout', 'match', dashes, or empty lines)
    while IFS=$'\t' read -r line; do
        # Skip header lines
        if [[ "$line" =~ ^psLayout ]] || [[ "$line" =~ ^match ]] || [[ "$line" =~ ^-+ ]] || [[ -z "$line" ]]; then
            continue
        fi
        
        # Parse the line - column 10 (index 9) contains the query name
        IFS=$'\t' read -ra fields <<< "$line"
        if [ ${#fields[@]} -ge 10 ]; then
            query_name="${fields[9]}"
            
            # Extract Target_Gene_ID using more robust parsing
            if [[ "$query_name" =~ gL_[0-9]+_([^_]+_[^_]+_[0-9]+)_lncRNA ]]; then
                target_id="${BASH_REMATCH[1]}"
                psl_targets["$target_id"]=1
            elif [[ "$query_name" =~ _((Hum_XLOC_[0-9]+))_ ]]; then
                target_id="${BASH_REMATCH[1]}"
                psl_targets["$target_id"]=1
            fi
        fi
    done
} < "$INPUT_PSL"

NUM_IDS=${#psl_targets[@]}
echo "Found $NUM_IDS unique Target_Gene_IDs in PSL file"

# Debug: Show first few IDs found (optional)
if [ $NUM_IDS -gt 0 ]; then
    echo "Sample Target_Gene_IDs found:"
    count=0
    for id in "${!psl_targets[@]}"; do
        echo "  - $id"
        ((count++))
        [ $count -ge 3 ] && break
    done
fi

echo "Processing CSV file and adding match_pc column..."

# Process CSV file
{
    # Read and modify header
    IFS=$'\t' read -r header
    echo -e "${header}\tmatch_pc"
    
    # Process data rows
    row_count=0
    matched_count=0
    unmatched_count=0
    
    while IFS=$'\t' read -r line; do
        # Extract Target_Gene_ID (first field)
        target_gene_id="${line%%$'\t'*}"
        
        # Check if ID exists in PSL targets
        if [[ -n "${psl_targets[$target_gene_id]}" ]]; then
            match_value="YES"
            ((matched_count++))
        else
            match_value="NO"
            ((unmatched_count++))
        fi
        
        # Output line with match value
        echo -e "${line}\t${match_value}"
        ((row_count++))
        
        # Progress indicator for large files
        if [ $((row_count % 1000)) -eq 0 ]; then
            echo "  Processed $row_count rows..." >&2
        fi
    done
    
    echo "Processing complete!" >&2
    echo "Total data rows processed: $row_count" >&2
    echo "Rows with matches: $matched_count" >&2
    echo "Rows without matches: $unmatched_count" >&2
    
} < "$INPUT_CSV" > "$OUTPUT_CSV"

echo "Output written to: $OUTPUT_CSV"

# Validation check
output_lines=$(wc -l < "$OUTPUT_CSV")
input_lines=$(wc -l < "$INPUT_CSV")

if [ "$output_lines" -eq "$input_lines" ]; then
    echo "✓ Validation passed: Output has same number of lines as input"
else
    echo "⚠ Warning: Line count mismatch (input: $input_lines, output: $output_lines)"
fi