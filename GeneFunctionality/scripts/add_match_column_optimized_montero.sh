#!/opt/local/bin/bash

# CORRECTED SCRIPT (v4)
# Fixes:
# 1. Removed the 'gL_' regex check that was causing all PSL lines to be ignored.
# 2. Uses PSL Column 10 directly as the matching key.
# 3. Ensures alignment between CSV Col 3 and PSL Col 10.

# Check bash version
if [ "${BASH_VERSION%%.*}" -lt 4 ]; then
    echo "Error: This script requires bash version 4 or higher"
    exit 1
fi

# Check arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_csv> <input_psl> <output_csv>"
    exit 1
fi

INPUT_CSV="$1"
INPUT_PSL="$2"
OUTPUT_CSV="$3"

# Check files
if [ ! -f "$INPUT_CSV" ]; then echo "Error: '$INPUT_CSV' not found!"; exit 1; fi
if [ ! -f "$INPUT_PSL" ]; then echo "Error: '$INPUT_PSL' not found!"; exit 1; fi

declare -A psl_strands
declare -A psl_col14

echo "Processing PSL file..."

# Process PSL file
{
    while IFS=$'\t' read -r line; do
        # Skip headers
        if [[ "$line" =~ ^psLayout ]] || [[ "$line" =~ ^match ]] || [[ "$line" =~ ^-+ ]] || [[ -z "$line" ]]; then
            continue
        fi
        
        IFS=$'\t' read -ra fields <<< "$line"
        
        # Check for at least 14 columns
        if [ ${#fields[@]} -ge 14 ]; then
            # Column 10 is the ID (e.g., human_lncrna_fused_10138_3_gRNA1)
            query_name="${fields[9]}"
            
            # Column 14 contains the target strand info
            target_col14="${fields[13]}" 
            
            # Extract strand from column 14
            psl_strand=""
            if [[ "$target_col14" == *"(+)" ]]; then
                psl_strand="+"
            elif [[ "$target_col14" == *"(-)" ]]; then
                psl_strand="-"
            fi
            
            # FIX: Use the query_name directly. 
            # Removed the 'if [[ =~ ^(gL_...) ]]' check which was filtering out your data.
            if [[ -n "$query_name" ]]; then
                psl_strands["$query_name"]="$psl_strand"
                psl_col14["$query_name"]="$target_col14"
            fi
        fi
    done
} < "$INPUT_PSL"

NUM_IDS=${#psl_strands[@]}
echo "Found $NUM_IDS unique gRNA_IDs in PSL file"

# Debug: Print first 3 keys to verify they look like your CSV data
if [ $NUM_IDS -gt 0 ]; then
    echo "Debug - First 3 IDs stored:"
    count=0
    for id in "${!psl_strands[@]}"; do
        echo "  Key: '$id' | Strand: '${psl_strands[$id]}'"
        ((count++))
        [ $count -ge 3 ] && break
    done
fi

echo "Processing CSV file..."

# Process CSV file
{
    # Read header
    IFS=$'\t' read -r header
    IFS=$'\t' read -ra header_fields <<< "$header"
    
    # Print Header + new columns
    # Note: We reconstruct the header to ensure clean tab separation
    if [ ${#header_fields[@]} -ge 7 ]; then
        echo -e "${header_fields[0]}\t${header_fields[1]}\t${header_fields[2]}\t${header_fields[3]}\t${header_fields[4]}\t${header_fields[5]}\t${header_fields[6]}\tmatch_pc\tpsl_target_col14"
    else
        echo -e "${header}\tmatch_pc\tpsl_target_col14"
    fi

    row_count=0
    matched_count=0
    unmatched_count=0
    
    while IFS=$'\t' read -r line; do
        IFS=$'\t' read -ra csv_fields <<< "$line"
        
        if [ ${#csv_fields[@]} -lt 6 ]; then
            continue
        fi
        
        # KEY MATCHING:
        # In your file1 sample, "Target_Gene_ID" (the one matching PSL) is Column 3.
        # Bash array index 2 = Column 3.
        gRNA_id="${csv_fields[2]}" 
        
        # Strand is Column 6 (index 5)
        csv_strand_raw="${csv_fields[5]}" 
        
        # Normalize CSV strand
        csv_strand_norm=""
        if [[ "$csv_strand_raw" == *"(+)"* ]] || [[ "$csv_strand_raw" == "+" ]]; then
            csv_strand_norm="+"
        elif [[ "$csv_strand_raw" == *"(-)"* ]] || [[ "$csv_strand_raw" == "-" ]]; then
            csv_strand_norm="-"
        fi
        
        match_value="NO"
        col14_value="NA"
        
        # Check lookup map
        if [[ -n "${psl_strands[$gRNA_id]}" ]]; then
            psl_strand_val="${psl_strands[$gRNA_id]}"
            col14_value="${psl_col14[$gRNA_id]}"
            
            # Compare normalized strands
            if [[ -n "$csv_strand_norm" && "$csv_strand_norm" == "$psl_strand_val" ]]; then
                match_value="YES"
                ((matched_count++))
            else
                ((unmatched_count++))
            fi
        else
            ((unmatched_count++))
        fi
        
        # Output reconstruction (handling 7 existing columns)
        if [ ${#csv_fields[@]} -ge 7 ]; then
            echo -e "${csv_fields[0]}\t${csv_fields[1]}\t${csv_fields[2]}\t${csv_fields[3]}\t${csv_fields[4]}\t${csv_fields[5]}\t${csv_fields[6]}\t${match_value}\t${col14_value}"
        else
            echo -e "${line}\t${match_value}\t${col14_value}"
        fi

        ((row_count++))
        if [ $((row_count % 5000)) -eq 0 ]; then echo "  Processed $row_count rows..." >&2; fi

    done
    
    echo "Finished. Rows matched: $matched_count | Rows unmatched: $unmatched_count" >&2
    
} < "$INPUT_CSV" > "$OUTPUT_CSV".tmp

# Final Awk replacement (Col 2 replacement)
awk 'BEGIN{
    FS="\t"; OFS="\t";
    map["Shared"] = "Core";
    map["Partially shared"] = "Common";
    map["Cell-type specific"] = "Specific";
}
{
    if ($2 in map) { $2 = map[$2]; }
    print
}' "$OUTPUT_CSV".tmp > "$OUTPUT_CSV"

rm -rf "$OUTPUT_CSV".tmp
echo "Output saved to $OUTPUT_CSV"