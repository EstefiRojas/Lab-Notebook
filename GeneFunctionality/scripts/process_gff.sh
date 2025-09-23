#!/bin/bash

# A bash script to find coordinates in a GFF3 file based on RNAcentral IDs,
# fetch the corresponding genomic sequences, and output the results to CSV files.

# --- Argument Check ---
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <gff_file> <ids_csv_file> <output_basename>"
    echo "Example: $0 sample.gff3 rnacentral_ids.csv my_analysis"
    exit 1
fi

GFF_FILE="$1"
IDS_FILE="$2"
OUTPUT_BASENAME="$3"

# --- Configuration ---
# Output file names are now constructed from the command-line argument.
EXON1_OUTPUT="${OUTPUT_BASENAME}_exon1_data.csv"
EXON2_OUTPUT="${OUTPUT_BASENAME}_exon2_data.csv"


# Check if the provided files exist
if [ ! -f "$GFF_FILE" ]; then
    echo "Error: GFF file not found at '$GFF_FILE'"
    exit 1
fi

if [ ! -f "$IDS_FILE" ]; then
    echo "Error: IDs CSV file not found at '$IDS_FILE'"
    exit 1
fi


# --- Temporary files ---
# Using mktemp for safer temporary file creation
TARGET_IDS_FILE=$(mktemp)
EXON1_COORDS=$(mktemp)
EXON2_COORDS=$(mktemp)

# --- Cleanup ---
# Ensure temporary files are removed when the script exits, even on error.
trap 'rm -f "$TARGET_IDS_FILE" "$EXON1_COORDS" "$EXON2_COORDS"' EXIT

# --- Script Start ---
echo "--- Starting Analysis ---"
echo "Using GFF file: $GFF_FILE"
echo "Using IDs file: $IDS_FILE"
echo "Output files will be named based on: $OUTPUT_BASENAME"

# --- Step 1: Extract target RNAcentral IDs from the CSV file ---
# - `cut -d, -f4 "$IDS_FILE"`: Extracts the 4th column (ENST IDs) using comma as a delimiter.
# - `tail -n +2`: Skips the header row.
# - The result is stored in a temporary file for fast lookup.
cut -d, -f4 "$IDS_FILE" | tail -n +2 > "$TARGET_IDS_FILE"
echo "Successfully extracted $(wc -l < "$TARGET_IDS_FILE") target IDs from '$IDS_FILE'."

# --- Step 2: Parse the GFF file to find matching exons ---
# We use 'awk' to do the heavy lifting in a single pass.
echo "Parsing GFF file: '$GFF_FILE'..."
awk -v exon1_coords="$EXON1_COORDS" -v exon2_coords="$EXON2_COORDS" '
    BEGIN {
        # Tell awk to use a tab as the field separator for the GFF file.
        FS="\t";
        # Load all target IDs into an associative array for O(1) lookups.
        while ( (getline id < "'$TARGET_IDS_FILE'") > 0 ) {
            targets[id] = 1;
        }
        close("'$TARGET_IDS_FILE'");
    }
    # For every line in the GFF file...
    !/^#/ && $3 == "exon" {
        # Split the attributes column (9th field) by semicolon.
        split($9, attrs, ";");
        name = "";
        id = "";
        # Loop through attributes to find "Name" and "ID".
        for (i in attrs) {
            if (attrs[i] ~ /^Name=/) {
                sub(/^Name=/, "", attrs[i]);
                name = attrs[i];
            }
            if (attrs[i] ~ /^ID=/) {
                sub(/^ID=/, "", attrs[i]);
                id = attrs[i];
            }
        }
        # If the Name is one of our targets...
        if (name in targets) {
            # Check if it is exon 1 or exon 2 and write coordinates to the respective temp file.
            if (id ~ /:ncRNA_exon1$/) {
                print $1, $4, $5, $7, name >> exon1_coords;
            } else if (id ~ /:ncRNA_exon2$/) {
                print $1, $4, $5, $7, name >> exon2_coords;
            }
        }
    }
' "$GFF_FILE"

echo "Found $(wc -l < "$EXON1_COORDS") matching records for exon 1."
echo "Found $(wc -l < "$EXON2_COORDS") matching records for exon 2."

# --- Step 3: Fetch sequences and create output CSV files ---

# Define a function to process a coordinates file.
# Arguments: $1=CoordinatesFile, $2=OutputFile, $3=ExonTypeName
process_coords_file() {
    local coords_file="$1"
    local output_file="$2"
    local exon_type_name="$3"
    local total_lines=$(wc -l < "$coords_file")
    local counter=1

    if [ "$total_lines" -eq 0 ]; then
        echo -e "\nNo coordinates found for $exon_type_name, skipping."
        return
    fi
    
    echo -e "\n--- Fetching sequences for ${exon_type_name} ---"

    # Write CSV header
    echo "ID,Functional,Chromosome,Start,End,Sequence,RNACentralID" > "$output_file"

    # Extract the exon number from the type name (e.g., "Exon 1" -> "1")
    local exon_number=$(echo "$exon_type_name" | awk '{print $2}')

    # Read the coordinates file line by line
    while read -r chromosome start end strand rnacentral_id; do
        echo "Fetching $exon_type_name ${counter}/${total_lines}: $rnacentral_id (chr${chromosome}:${start}-${end})"
        
        # Convert strand from GFF format (+/-) to Ensembl format (1/-1)
        if [ "$strand" == "+" ]; then
            ensembl_strand=1
        else
            ensembl_strand=-1
        fi

        # Construct the API URL
        url="https://rest.ensembl.org/sequence/region/human/${chromosome}:${start}..${end}:${ensembl_strand}?content-type=text/plain"

        # Fetch the sequence using curl, with silent mode and fail-fast
        sequence=$(curl -s -f "$url")
        
        # If curl fails, set a placeholder message
        if [ $? -ne 0 ]; then
            sequence="SEQUENCE_NOT_FOUND"
        fi

        # Construct the new ID in the format: URS..._exon_1
        local new_id="${rnacentral_id}_exon_${exon_number}"

        # Write the formatted data to the CSV file
        echo "${new_id},NA,chr${chromosome},${start},${end},${sequence},${rnacentral_id}" >> "$output_file"

        # Be polite to the API server
        sleep 0.2
        ((counter++))
    done < "$coords_file"
    echo "Successfully wrote $(wc -l < "$coords_file") records to '$output_file'."
}

# Process both exon files
process_coords_file "$EXON1_COORDS" "$EXON1_OUTPUT" "Exon 1"
process_coords_file "$EXON2_COORDS" "$EXON2_OUTPUT" "Exon 2"

echo -e "\n--- Analysis Complete ---"


