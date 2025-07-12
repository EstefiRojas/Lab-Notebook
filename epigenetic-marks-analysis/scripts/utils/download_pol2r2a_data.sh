#!/bin/bash
set -uex

echo "Starting POL2R2A marks download..."

for file in $(find ../data/POL2R2A/ -type f -print | grep .txt); do
    mkdir -p "${file%.txt}"  # Create directory if not present

    for link in $(cat "$file" | grep bed.gz); do
        code=$(basename "$link" .bed.gz)
        if [ ! -f "${file%.txt}"/"$code".bed.gz ]; then
            wget -q "$link" -O "${file%.txt}"/"$code".bed.gz
            echo "Download successful: " "$code"
        else
            echo "File present: $code"
        fi
    done
    
done
echo "Finished downloading POL2R2A marks."

