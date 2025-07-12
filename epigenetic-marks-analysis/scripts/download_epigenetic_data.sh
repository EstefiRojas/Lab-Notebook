#!/bin/bash
#
# Dependencies:
# - wget: Required for downloading files
# - awk: Required for extracting links from metadata
#
# Script Name: download_epigenetic_data.sh
#
# Author: Estefania Rojas
#
# Description: This script downloads epigenetic data in narrowPeak bed format
# from ENCODE database by finding all links provided in the folder
# data/epigenetic_data/.
# The script implements robust retry logic and basic file verification
# (exists, readable, not empty) to ensure downloaded files are present. If
# files are already present and seem valid (not empty), the script skips them.
#
###############################################################################

# Bash strict mode. -e exit on error, -u error on unset variables, -o pipefail
# ensures pipeline errors are caught
set -uo pipefail

# Maximum number of download attempts per file
MAX_ATTEMPTS=2

# --- Function Definitions ---

# Function to verify basic integrity of a downloaded file (exists, readable, not empty)
# Returns:
#   0: File exists, is readable, and is not empty
#   1: File exists but is empty or not readable
#   2: File does not exist
verify_bed_file() {
    local file="$1"

    # 1. Check if file exists
    if [ ! -f "$file" ]; then
        # No need to print here, the main loop handles messages for non-existent files
        return 2 # Specific code for file not found
    fi

    # 2. Check if file is readable
    if [ ! -r "$file" ]; then
        echo "i Verification failed: File not readable '$file'" >&2
        return 1 # Generic failure code
    fi

    # 3. Check if file is empty
    if [ ! -s "$file" ]; then
        echo "i Verification failed: File is empty '$file'" >&2
        return 1 # Generic failure code
    fi

    # If all checks passed:
    return 0 # File exists, is readable, and not empty
}

# Function to download a file with retry logic
download_with_retry() {
    local url="$1"
    local output_file="$2"
    local attempt=1
    local success=false

    while [ $attempt -le $MAX_ATTEMPTS ] && [ "$success" = false ]; do
        echo "Download attempt $attempt for $(basename "$output_file")..."

        # Use wget to download the file (-q for quiet, -O for output file)
        # Added --tries=1 and --timeout=60 for wget specific control per attempt
        if wget -q --tries=1 --timeout=60 "$url" -O "$output_file"; then
            # Verify if the downloaded file exists and is not empty
            if verify_bed_file "$output_file"; then
                echo "✓ Successfully downloaded and verified (exists, not empty): $(basename "$output_file")"
                success=true
            else
                # verify_bed_file already printed an error message
                echo "✗ Downloaded file '$output_file' is invalid (empty/unreadable, check details above), retrying..."
                rm -f "$output_file" # Remove invalid file before retrying
            fi
        else
            # Wget failed (e.g., network error, 404, timeout)
            local wget_exit_code=$?
            echo "✗ Download failed (wget exit code $wget_exit_code), retrying..."
            # Ensure partially downloaded file is removed if wget failed midway
            rm -f "$output_file"
        fi

        # Increase attempt counter
        attempt=$((attempt + 1))

        # If not successful and not the last attempt, wait before retrying
        if [ "$success" = false ] && [ $attempt -le $MAX_ATTEMPTS ]; then
            sleep 2 # Wait 2 seconds before retrying
        fi
    done

    # Return success status
    if [ "$success" = true ]; then
        return 0
    else
        echo "! Failed to download and verify after $MAX_ATTEMPTS attempts: $(basename "$output_file")"
        # Ensure any potentially corrupted file from the last attempt is removed
        rm -f "$output_file"
        return 1
    fi
}

# --- Main Script Logic ---

echo "Starting Epigenetic marks download..."
echo "=================================================="
# --- DEBUG: Show current working directory ---
# echo "[DEBUG] Script running from: $(pwd)"
# --- END DEBUG ---

# Counters for reporting
total_files=0
skipped_files=0
successful_downloads=0
failed_downloads=0

# Check if epigenetic type argument is provided
if [ -z "${1:-}" ]; then
    echo "Error: Please provide the epigenetic data type (e.g., histone_marks, chromatin_accessibility, methylome) as the first argument." >&2
    exit 1
fi
epigenetic_type=$1
# echo "[DEBUG] Epigenetic type provided: $epigenetic_type" # DEBUG

# Define the base directory for input files
input_base_dir="../data/datasets/epigenetic_data/${epigenetic_type}/"
# --- DEBUG: Show the calculated input base directory ---
# echo "[DEBUG] Input base directory set to: $input_base_dir"
# --- END DEBUG ---


# Check if the input directory exists
if [ ! -d "$input_base_dir" ]; then
    echo "Error: Input directory not found: $input_base_dir" >&2
    # --- DEBUG: Confirm directory check failure ---
    # echo "[DEBUG] Directory check failed for: $input_base_dir"
    # --- END DEBUG ---
    exit 1
fi
# --- DEBUG: Confirm directory check success ---
# echo "[DEBUG] Input directory exists: $input_base_dir"
# --- END DEBUG ---

# Define biomarker types and histone marks (only used if epigenetic_type is histone_marks)
biomarker_types=("cell+line" "in+vitro+differentiated+cells" "primary+cell" "tissue")
histone_marks=("H2AFZ" "H2AK5ac" "H2AK9ac" "H2BK120ac" "H2BK12ac" "H2BK15ac" \
               "H2BK20ac" "H2BK5ac" "H3F3A" "H3K14ac" "H3K18ac" "H3K23ac" \
               "H3K23me2" "H3K27ac" "H3K27me3" "H3K36me3" "H3K4ac" "H3K4me1" \
               "H3K4me2" "H3K4me3" "H3K56ac" "H3K79me1" "H3K79me2" "H3K9ac" \
               "H3K9me1" "H3K9me2" "H3K9me3" "H3T11ph" "H4K12ac" "H4K20me1" \
               "H4K5ac" "H4K8ac" "H4K91ac")


# Optional: Generate link files if the type is histone_marks
# Corrected string comparison and loop syntax
if [ "${epigenetic_type}" == "histone_marks" ]; then
    echo "Checking/Generating download link files for histone marks..."
    for histone in "${histone_marks[@]}"; do
        for biomarker in "${biomarker_types[@]}"; do
            links_dir="../data/datasets/epigenetic_data/histone_marks/${histone}"
            links_file="${links_dir}/${biomarker}.txt"
            metadata_file="${links_dir}/${biomarker}.tsv" # Define metadata file path

            mkdir -p "$links_dir" # Ensure histone directory exists

            # Check if links file already exists
            if [ -f "$links_file" ]; then
                echo "✓ Links file already present: ${histone}/${biomarker}.txt"
                continue
            fi

            echo "Generating links file for: ${histone}/${biomarker}.txt"
            # Metadata url - Ensure URL encoding is correct for '+' in biomarker
            # '+' should typically be %2B in URLs, but ENCODE might handle '+' directly. Test if needed.
            ENCODE_METADATA_URL="https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=Histone+ChIP-seq&assay_title=Mint-ChIP-seq&status=released&biosample_ontology.classification=${biomarker}&assembly=GRCh38&files.file_type=bed+narrowPeak&target.label=${histone}&type=Experiment&files.analyses.status=released&files.preferred_default=true"

            # Download metadata, exit script if download fails
            echo "Downloading metadata for ${histone}/${biomarker}..."
            if ! wget -q "${ENCODE_METADATA_URL}" -O "$metadata_file"; then
                 echo "Error: Failed to download metadata from ${ENCODE_METADATA_URL}" >&2
                 # Clean up potentially empty/partial metadata file
                 rm -f "$metadata_file"
                 # Decide whether to continue with other files or exit
                 # For now, let's skip this specific file generation and continue
                 echo "Skipping link generation for ${histone}/${biomarker} due to metadata download failure."
                 continue # Skip to the next biomarker/histone combination
                 # Alternatively, exit 1 if metadata is critical
            fi

            # Check if metadata file is empty or contains errors (e.g., HTML error page)
            if [ ! -s "$metadata_file" ] || grep -q -i "<error>" "$metadata_file"; then
                 echo "Error: Metadata file '$metadata_file' is empty or contains errors." >&2
                 rm -f "$metadata_file" # Clean up bad metadata file
                 continue # Skip this file
            fi

            # Save url to first line of link files
            echo "# Metadata URL: ${ENCODE_METADATA_URL}" > "$links_file"
            echo "# Generated on: $(date)" >> "$links_file"

            # Filter links from metadata (assuming column 48 is 'File download URL')
            # Check if the awk command produces any output
            # Updated to handle potential .bed.gz files as well
            if ! awk -F'\t' 'NR > 1 && $6 == "GRCh38" && $56 == "released" && $48 ~ /\.bed(\.gz)?$/ {print $48}' "$metadata_file" >> "$links_file"; then
                 echo "Warning: awk command failed or produced no output for $metadata_file. Links file might be incomplete." >&2
                 # Keep the file with header comments even if no links found/awk fails
            fi
            # Optional: Remove the metadata file after extracting links if not needed
            # rm -f "$metadata_file"
            echo "✓ Generated links file: ${histone}/${biomarker}.txt"
        done
    done
    echo "Finished checking/generating histone mark link files."
fi

echo "Searching for .txt files containing download links in: $input_base_dir"
# --- DEBUG: Manually run find to see output before the loop ---
# echo "[DEBUG] Running find command: find \"$input_base_dir\" -type f -name \"*.txt\" -print"
# find "$input_base_dir" -type f -name "*.txt" -print || echo "[DEBUG] find command failed or found nothing."
# echo "[DEBUG] Starting main processing loop..."
# --- END DEBUG ---


# Find all .txt files containing download links and loop through each
# Use process substitution <(...) to avoid running the while loop in a subshell
while IFS= read -r -d $'\0' file; do
    # --- DEBUG: Show which file is being processed by the loop ---
    # echo "[DEBUG] Processing file found by find: $file"
    # --- END DEBUG ---

    # Create output directory based on the txt filename (without extension)
    output_dir="${file%.txt}"
    mkdir -p "$output_dir"

    echo "Processing links from: $(basename "$file")"

    # Get all lines that look like URLs ending with .bed or .bed.gz from the txt file
    # Skip commented lines
    # Using grep with -E for extended regex. Handle case where grep finds nothing.
    links=$(grep -E -v '^#' "$file" | grep -Eo 'https://[^[:space:]"]+\.bed(\.gz)?' || true)
    # --- DEBUG: Show extracted links ---
    # if [ -z "$links" ]; then
    #     echo "[DEBUG] No links extracted from $file"
    # else
    #     echo "[DEBUG] Extracted links from $file:"
    #     echo "$links" # Print the links found
    # fi
    # --- END DEBUG ---

    # If no links found, continue to next file
    if [ -z "$links" ]; then
        echo "No valid download links found in $(basename "$file")"
        echo "--------------------------------------------------"
        continue
    fi

    # Count the total number of links found in this file
    link_count=$(echo "$links" | wc -l)
    echo "Found $link_count links to process"

    # Process each link using process substitution
    # This inner loop is fine with process substitution as it doesn't modify the main counters directly
    while IFS= read -r link; do
        # --- DEBUG: Show which link is being processed ---
        # echo "[DEBUG] Processing link: $link"
        # --- END DEBUG ---

        # Skip empty lines that might result from parsing
        if [ -z "$link" ]; then
            # echo "[DEBUG] Skipping empty line." # DEBUG
            continue
        fi

        total_files=$((total_files + 1))
        # --- DEBUG: Show counter increment ---
        # echo "[DEBUG] total_files incremented to: $total_files"
        # --- END DEBUG ---

        # Handle potential .gz extension in output filename
        code=$(basename "$link")
        output_file="${output_dir}/${code}"

        # Check if file already exists and is valid (exists, readable, not empty)
        # Initialize verify_status assuming file doesn't exist
        verify_status=2
        if [ -f "$output_file" ]; then
             # File exists, now verify it using the simplified function
             verify_bed_file "$output_file"
             verify_status=$? # Capture return code (0=valid, 1=invalid, 2=not found - though already checked)
        fi
        # If verify_bed_file returned 2, it means the file check inside the function failed (e.g. permissions)
        # Treat this as an error/invalid state for the purpose of re-downloading.
        if [ $verify_status -eq 2 ]; then
             # Message now indicates file not found or inaccessible
             echo "i File '$output_file' not found or inaccessible, proceeding to download."
             # Set status to trigger download logic below
             verify_status=1 # Treat as if invalid to trigger download logic
        fi


        if [ $verify_status -eq 0 ]; then
            # Status 0 means file exists, is readable and not empty
            echo "✓ File already present and valid (exists, not empty): $code"
            skipped_files=$((skipped_files + 1))
            # --- DEBUG: Show counter increment ---
            # echo "[DEBUG] skipped_files incremented to: $skipped_files"
            # --- END DEBUG ---
            continue # Skip to the next link
        elif [ $verify_status -eq 1 ]; then
             # Status 1 means file exists but is invalid (empty/unreadable)
             echo "i File exists but is invalid (empty/unreadable, check details above), attempting re-download: $code"
             rm -f "$output_file" # Remove the invalid file before downloading
        fi
        # If verify_status was 1 (file was invalid and removed), proceed to download

        # Download the file with retry logic
        if download_with_retry "$link" "$output_file"; then
            successful_downloads=$((successful_downloads + 1))
            # --- DEBUG: Show counter increment ---
            # echo "[DEBUG] successful_downloads incremented to: $successful_downloads"
            # --- END DEBUG ---
            # Optional: If the downloaded file is gzipped, unzip it
            # Note: If you uncomment this, ensure `gunzip` is listed as a dependency
            # if [[ "$output_file" == *.gz ]]; then
            #    echo "Unzipping $output_file..."
            #    # Use gunzip -f to force overwrite if an unzipped version somehow exists
            #    if gunzip -f "$output_file"; then
            #        echo "✓ Successfully unzipped."
            #        # Optional: You might want to verify the unzipped .bed file now
            #        # unzipped_file="${output_file%.gz}"
            #        # verify_bed_file "$unzipped_file" # Using the simplified check
            #    else
            #        echo "! Failed to unzip $output_file" >&2
            #        # Decide how to handle unzip failure - maybe increment failed_downloads?
            #        # Consider removing the failed .gz file: rm -f "$output_file"
            #    fi
            # fi
        else
            failed_downloads=$((failed_downloads + 1))
            # --- DEBUG: Show counter increment ---
            # echo "[DEBUG] failed_downloads incremented to: $failed_downloads"
            # --- END DEBUG ---
        fi

    done < <(echo "$links") # Use process substitution to feed links

    echo "--------------------------------------------------"
# Redirect the output of find (using process substitution) into the while loop
done < <(find "$input_base_dir" -type f -name "*.txt" -print0)

# --- DEBUG: Indicate that the main processing loop has finished ---
# echo "[DEBUG] Finished main processing loop."
# --- END DEBUG ---


echo "=================================================="
echo "Finished downloading Epigenetic marks."
echo "Summary:"
echo "  Total links processed: $total_files"
echo "  Already present and valid (exists, not empty): $skipped_files"
echo "  Successfully downloaded/verified: $successful_downloads"
echo "  Failed to download/verify: $failed_downloads"
echo "=================================================="

# Exit with success code if no failures, otherwise exit with an error code
if [ $failed_downloads -eq 0 ]; then
    exit 0
else
    exit 1
fi
