#!/bin/bash

set -uex

# Check for correct usage
if [ $# -ne 3 ]; then
    echo "Usage: $0 feature_to_add_file.csv features_file.csv postfix"
    exit 1
fi

# Store input file
FEATURE_TO_ADD=$1
FEATURES_FILE=$2
POSTFIX=$3

sort -t, -k4,4 -V "$FEATURE_TO_ADD" | \
cut -d',' -f6  | \
paste -d, "$FEATURES_FILE" - > "${FEATURES_FILE%.csv}"-"$POSTFIX".csv

