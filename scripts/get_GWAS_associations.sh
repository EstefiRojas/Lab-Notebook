#!/bin/bash
# This script retrieves the minimum p-value of SNPs present at lncRNAs from 
# GWAS summary-statistics API using corresponding ENSEMBL coordinates as filter.
# It also filters by p-value < 1e-5.
#
# Dependencies: 
# - jq:
#       https://jqlang.github.io/jq/
#
# Script Name: get_GWAS_associations.sh
#
# Author: Estefania Rojas
#
# Input: 1. Csv file: containing the dataset of sequences to be analysed. This 
# 			file should have the exon 1 and 2 coordinates at columns 3-9.
#           Ex. ../data/lncRNAs_ensembl_exon_coords.csv
#        2. Name of the output csv file: minimum p-value for each coordinate in 
# 			input file 1.
#           Ex. lncrna-gwas-pval-feature
#
###############################################################################
# Bash strict mode. -e report all errors, -u exit on error, -x print every command executed.
set -x

# Get the date of execution. Format example: Tuesday January 17, 2024
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}

# Check for correct usage
if [ $# -ne 2 ]; then
    echo "Usage: $0 regions_to_extract.csv outname.csv"
    exit 1
fi
##############################################################
#Function to convert X to 23 and Y to 24 while leaving numeric values unchanged
convert_chromosome() {
    local chrm=$1
    
    case "$chrm" in
        "X"|"x")
            echo "23"
            ;;
        "Y"|"y")
            echo "24"
            ;;
        *)
            echo "$chrm"
            ;;
    esac
}

# Function to get unique GWAS variants given chromosome and position range
get_unique_gwas_variants() {
   local chrm=$(convert_chromosome "$1")
   local start_pos=$2
   local end_pos=$3
   
   # Initialize empty array for storing variant IDs and p-values
   local variant_ids=()
   local p_values=()

   # If X or Y chromosome, convert to 23 or 24, respectibly.

   
   # Initialize current url
   local current_url="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/${chrm}/associations?bp_lower=${start_pos}&bp_upper=${end_pos}&p_upper=1e-5&size=40"

   while true; do
       # Create temp file for headers
       local header_file=$(mktemp)
       
       # Query current page and save headers separately
       local gwas_response=$(curl -s -D "${header_file}" --retry 2 "${current_url}" -H "Accept:application/json")
       
       # Get HTTP status from headers
       local http_code=$(grep "HTTP/" "${header_file}" | cut -d' ' -f2)
       rm "${header_file}"
       
       # Check for 404 or other errors
       if [ "$http_code" != "200" ]; then
           break
       fi
       
       # Extract variant_ids from current page and add to array
       local page_variants=$(echo "$gwas_response" | jq -r '._embedded.associations | to_entries[] | .value.variant_id' 2>/dev/null)
       if [ ! -z "$page_variants" ]; then
           while IFS= read -r variant; do
               variant_ids+=("$variant")
           done <<< "$page_variants"
       fi

       # Extract pvalues from current page and add to array
       local page_pvalues=$(echo "$gwas_response" | jq -r '._embedded.associations | to_entries[] | .value.p_value' 2>/dev/null)
       if [ ! -z "$page_pvalues" ]; then
           while IFS= read -r p_val; do
               p_values+=("$p_val")
           done <<< "$page_pvalues"
       fi
       
       # Check for next page, return empty if not found
       local next_url=$(echo "$gwas_response" | jq -r '._links.next.href // empty')
       
       # Break if no next page
       if [ -z "$next_url" ]; then
           break
       fi
       
       # Update URL for next iteration
       current_url=$next_url
   done

   # Count unique variants, properly handling empty array
   local unique_count=0
   if [ ${#variant_ids[@]} -gt 0 ]; then
       unique_count=$(echo "${variant_ids[@]}" | tr ' ' '\n' | sort | uniq | wc -l | tr -d ' ')
   fi
   
   # Initialize min_p_value with first valid p-value
   local min_p_value=""
   for p_value in "${p_values[@]}"; do
       if [ "$p_value" != "-99" ]; then
           min_p_value=$p_value
           break
       fi
   done

   # If we found a valid initial p-value, continue searching for minimum
   if [ ! -z "$min_p_value" ]; then
       for p_value in "${p_values[@]}"; do
           if [ "$p_value" != "-99" ] && (( $(echo "$p_value < $min_p_value" | bc -l) )); then
               min_p_value=$p_value
           fi
       done
   else
   	   min_p_value='"NA"'
   fi

   # Return both values separated by a comma
   echo -e "${unique_count},${min_p_value}"
}

# Store input variables
REGIONS_FILE=$1
OUTPUT_NAME=$2

echo "Processing file $REGIONS_FILE"

# Define path to store output
output_path=../data/datasets/gwas_pval_feature

# Check previous runs
echo "Checking for previous runs..."
start_line=1
if [ -w $output_path/$OUTPUT_NAME.csv ]; then
	echo "Found a previous run with specified parameters."
	start_line=$(wc -l $output_path/$OUTPUT_NAME.csv | awk '{ print $1 }')
fi

# Create output directory if needed
mkdir -p "$output_path"

echo "Output file: $output_path/$OUTPUT_NAME.csv"

# Add header to output file
# Add header to output file if new run
if [ ! -e $output_path/$OUTPUT_NAME.csv ]; then
	echo "tl_exon1_gwas_associations,tl_exon1_min_p_value,tl_exon2_gwas_associations,tl_exon2_min_p_value,tr_exon1_gwas_associations,tr_exon1_min_p_value,tr_exon2_gwas_associations,tr_exon2_min_p_value" > "$output_path"/"$OUTPUT_NAME".csv
fi

# Loop through all ensembl coordinates
while read -r line; do
	echo "Processing $line"
	# Parse coordinates
	IFS=' ' read -r tl_chrm tl_start tl_end tl_exon1_start tl_exon1_end tl_exon2_start tl_exon2_end \
					tr_chrm tr_start tr_end tr_exon1_start tr_exon1_end tr_exon2_start tr_exon2_end <<< "$line"

	# Transcript 1, exon 1
	if [[ "$tl_exon1_start" != "NA" ]]; then
		tl_total_associations_exon1=$(get_unique_gwas_variants $(echo "${tl_chrm}" | tr -d '"') ${tl_exon1_start} ${tl_exon1_end})
	else
		tl_total_associations_exon1='"NA","NA"'
	fi
	echo "Found $tl_total_associations_exon1 associations for transcript1 exon1"

	# Transcript 1, exon 2
	if [[ "$tl_exon2_start" != "NA" ]]; then
		tl_total_associations_exon2=$(get_unique_gwas_variants $(echo "${tl_chrm}" | tr -d '"') ${tl_exon2_start} ${tl_exon2_end})
	else
		tl_total_associations_exon2='"NA","NA"'
	fi
	echo "Found $tl_total_associations_exon2 associations for transcript1 exon2"

	if [[ "$tr_exon1_start" != "NA" ]]; then
		tr_total_associations_exon1=$(get_unique_gwas_variants $(echo "${tr_chrm}" | tr -d '"') ${tr_exon1_start} ${tr_exon1_end})
	else
		tr_total_associations_exon1='"NA","NA"'
	fi
	echo "Found $tr_total_associations_exon1 associations for transcript2 exon1"

	if [[ "$tr_exon2_start" != "NA" ]]; then
		tr_total_associations_exon2=$(get_unique_gwas_variants $(echo "${tr_chrm}" | tr -d '"') ${tr_exon2_start} ${tr_exon2_end})
	else
		tr_total_associations_exon2='"NA","NA"'
	fi
	echo "Found $tr_total_associations_exon2 associations for transcript2 exon2"


    # Send to output file
    echo "${tl_total_associations_exon1},${tl_total_associations_exon2},${tr_total_associations_exon1},${tr_total_associations_exon2}" >> "$output_path/$OUTPUT_NAME.csv"


done < <(awk -F, -v start="$start_line" 'NR > start { print $3,$4,$5,$6,$7,$8,$9,$11,$12,$13,$14,$15,$16,$17 }' "$REGIONS_FILE") # Process only remaining lines.
