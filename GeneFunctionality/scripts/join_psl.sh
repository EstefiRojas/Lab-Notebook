#!/bin/bash

# Input Files
FILE1=$1
FILE2=$2
OUTPUT=$3

# Header for CSV
echo "Key,ENSG_ID,Col14_File1,Col14_File2,Result_YesNo" > "$OUTPUT"

# AWK command to process files
awk '
    # ---------------------------------------------------------
    # Processing File 1
    # ---------------------------------------------------------
    NR==FNR {
        # Store Column 14 (value) keyed by Column 10 (id)
        # $10 is the key, $14 is the value containing strand info
        file1_data[$10] = $14
        next
    }

    # ---------------------------------------------------------
    # Processing File 2
    # ---------------------------------------------------------
    ($10 in file1_data) {
        # Retrieve stored Col 14 from File 1
        col14_f1 = file1_data[$10]
        col14_f2 = $14

        # 1. Extract Strand from File 1 (Last 3 chars: (+) or (-))
        # We look for (+) or (-) at the end of the string
        match(col14_f1, /\([+-]\)$/)
        strand1 = substr(col14_f1, RSTART+1, 1) # Extract just + or -

        # 2. Extract Strand from File 2
        match(col14_f2, /\([+-]\)$/)
        strand2 = substr(col14_f2, RSTART+1, 1) # Extract just + or -

        # 3. Extract ENSG ID from File 2
        # Looks for pattern ENSG followed by digits
        match(col14_f2, /ENSG[0-9]+/)
        ensg_id = substr(col14_f2, RSTART, RLENGTH)

        # 4. Compare Strands
        # "depending on whether both strands are different or not... YES or NO respectively"
        # Different -> YES, Same -> NO
        if (strand1 != strand2) {
            result = "YES"
        } else {
            result = "NO"
        }

        # 5. Output to CSV format
        # OFS is comma
        print $10, ensg_id, col14_f1, col14_f2, result
    }
' OFS="," "$FILE1" "$FILE2" >> "$OUTPUT"

echo "Processing complete. Results saved to $OUTPUT"