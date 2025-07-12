#!/bin/bash
#
###############################################################################
# Get the date of execution. Format example: Tuesday January 17, 2024
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}

# Check for correct usage
if [ $# -ne 1 ]; then
    echo "Usage: $0 regions_to_extract.csv"
    exit 1
fi
##############################################################

INPUT1=$1

awk -F, 'NR > 1 {
    # Thoroughly clean the strings - especially handling newlines
    gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", $3)
    gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", $4)
    
    # Debug version - uncomment to see cleaned values
    # print "CLEANED: Col3=[" $3 "], Col4=[" $4 "], Equal=" ($3 == $4)
    
    if (length($3) > 0 && $3 == $4) {
        print $3 
    } else if (length($3) == 0) {
        print $4
    } else if (length($4) == 0) {
        print $3
    } else {
        print $3"|"$4
    }
}' "$INPUT1" > "unique_ids.csv"

