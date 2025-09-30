#!/bin/bash

# Bash script to perform a two-stage join and a final annotation step.
# 1. Processes a PSL file to filter, clean, and reformat columns. It specifically
#    wraps the comma-containing list fields (blockSizes, etc.) in quotes.
# 2. Joins the result with a CSV data file based on a composite key (GeneID, TranscriptID),
#    appending the 'Probability_Functional' column.
# 3. Annotates the result from step 2 by parsing IDs from a third file and the data stream
#    to find a common identifier and determine if a row is "ESSENTIAL".
#
# Usage: ./join_blat_matches_Montero.sh <original.psl> <data.csv> <essentials.csv> > final_output.csv

# --- Input Files ---
ORIGINAL_PSL_FILE="$1"  # The first input file (PSL format)
DATA_FILE="$2"         # The second input file (CSV format with functional probability)
ESSENTIALS_FILE="$3"   # The third file for annotation (CSV format with essential genes)

# --- Validate Input ---
# Check if all three file arguments were provided.
if [[ -z "$ORIGINAL_PSL_FILE" || -z "$DATA_FILE" || -z "$ESSENTIALS_FILE" ]]; then
  echo "Usage: $0 <original.psl> <data.csv> <essentials.csv>"
  echo "Please provide paths to all three input files."
  exit 1
fi

# Check if the provided files actually exist.
if [[ ! -f "$ORIGINAL_PSL_FILE" ]]; then echo "Error: File not found: $ORIGINAL_PSL_FILE"; exit 1; fi
if [[ ! -f "$DATA_FILE" ]]; then echo "Error: File not found: $DATA_FILE"; exit 1; fi
if [[ ! -f "$ESSENTIALS_FILE" ]]; then echo "Error: File not found: $ESSENTIALS_FILE"; exit 1; fi


# --- Header Row ---
# Define the header for the final output CSV by combining headers from all stages.
HEADER_JOIN_1_2="match,mis-match,rep.match,N_s,Q_gap_count,Q_gap_bases,T_gap_count,T_gap_bases,strand,Q_name_prefix,Q_name_suffix,Q_size,Q_start,Q_end,GeneID,TranscriptID,ExonNum,T_size,T_start,T_end,block_count,blockSizes,qStarts,tStarts,Probability_Functional"
# Header for the new annotation column from the third file.
HEADER_ESSENTIAL="Essential_Status"
# Combine headers into the final version.
FINAL_HEADER="${HEADER_JOIN_1_2},${HEADER_ESSENTIAL}"

# Print the final header row as the first line of the output.
echo "$FINAL_HEADER"

# --- Main Process ---
# This script uses a pipeline of three 'awk' commands.
# Each stage's output is piped directly to the next.

# First stage: Process the original PSL file.
awk '
  BEGIN { OFS="," }
  FNR > 5 {
      # The fields in a PSL file are space-delimited. Default awk behavior handles this.
      if ($2 == 0 && $3 == 0 && $4 == 0 && $5 == 0 && $6 == 0 && $7 == 0 && $8 == 0) {
          q_name_original = $10
          num_q_parts = split(q_name_original, qname_parts, "_")
          q_name_prefix = ""
          q_name_suffix = ""
          if (num_q_parts >= 2) {
              q_name_prefix = qname_parts[1] "_" qname_parts[2]
          } else {
              q_name_prefix = q_name_original
          }
          for (k=3; k < num_q_parts; k++) {
              q_name_suffix = q_name_suffix (k==3 ? "" : "_") qname_parts[k]
          }
          if (q_name_suffix == "") { q_name_suffix = "NA" }

          t_name_original = $14
          gene_id = transcript_id = exon_num = "NA"
          if (split(t_name_original, parts_slash, "/") >= 1) {
              gene_id = parts_slash[1]
              if (split(parts_slash[2], parts_underscore, "_") >= 1) {
                  transcript_id = parts_underscore[1]
                  exon_num = parts_underscore[2]
              }
          }

          # Print fields one by one to construct the CSV line.
          printf "%s", $1
          for (i=2; i<=9; i++) printf "%s%s", OFS, $i
          printf "%s%s%s%s", OFS, q_name_prefix, OFS, q_name_suffix
          for (i=11; i<=13; i++) printf "%s%s", OFS, $i
          printf "%s%s%s%s%s%s", OFS, gene_id, OFS, transcript_id, OFS, exon_num
          # Print fields from T_size up to block_count.
          for (i=15; i<=18; i++) printf "%s%s", OFS, $i
          
          # **FIX**: Print the last three fields (blockSizes, qStarts, tStarts), which contain
          # commas, by wrapping them in double quotes to create valid CSV output.
          # Original PSL fields are 19, 20, 21.
          printf "%s\"%s\"", OFS, $19 # blockSizes
          printf "%s\"%s\"", OFS, $20 # qStarts
          printf "%s\"%s\"", OFS, $21 # tStarts
          
          printf "\n" # End the line.
      }
  }' "$ORIGINAL_PSL_FILE" | \
# Second stage: Join processed PSL data with the functional data file.
awk -F, -v OFS=, '
  BEGIN { sep = SUBSEP }
  NR == FNR {
      if (FNR > 1) {
          gene_id_d = $1; gsub(/\r$/, "", gene_id_d)
          split($2, trans_parts, "|"); trans_id_part1 = trans_parts[1]; gsub(/\r$/, "", trans_id_part1)
          split($4, prob_parts, "|"); prob_part1 = prob_parts[1]; gsub(/\r$/, "", prob_part1)
          key = gene_id_d sep trans_id_part1
          data_map[key] = prob_part1
      }
      next
  }
  {
      gene_id_pipe = $15; transcript_id_pipe = $16
      key = gene_id_pipe sep transcript_id_pipe
      if (key in data_map) {
          gsub(/\r$/, "", $0)
          print $0, data_map[key]
      }
  }' "$DATA_FILE" - | \
# Third stage: Annotate the result by parsing IDs from the essentials file and the data stream.
awk -F, -v OFS=, '
  # First, read the essentials file and build a map of parsed IDs.
  # NR==FNR is true only while reading the first file (ESSENTIALS_FILE).
  NR == FNR {
      if (FNR > 1) { # Skip the header row.
          # The full ID is in the second column ("unified_ID").
          full_id = $2
          gsub(/\r$/, "", full_id) # Clean the ID.

          # Extract the "fused_NUMBER" part to use as a key.
          # e.g., from "human_lncrna_fused_10388" extract "fused_10388"
          if (match(full_id, /fused_[0-9]+/)) {
              parsed_key = substr(full_id, RSTART, RLENGTH)
              essentials_map[parsed_key] = 1 # Store the parsed key in the map.
          }
      }
      next # Move to the next line of ESSENTIALS_FILE.
  }

  # This block runs for the piped input (output from the second join).
  # It performs a "left join" by processing every line from the pipe.
  {
      # The field to check is the 11th column ("Q_name_suffix").
      field_to_check = $11
      gsub(/\r$/, "", field_to_check) # Clean the field.

      parsed_key_pipe = "" # Initialize parsed key for this line
      # Extract the "fused_NUMBER" part from the suffix.
      # e.g., from "fused_10006_3" extract "fused_10006"
      if (match(field_to_check, /fused_[0-9]+/)) {
          parsed_key_pipe = substr(field_to_check, RSTART, RLENGTH)
      }

      gsub(/\r$/, "", $0) # Clean the full line from the pipe.

      # Check if the parsed key from our data exists in the map of essential IDs.
      if (parsed_key_pipe != "" && parsed_key_pipe in essentials_map) {
          # If it exists, print the line and append "ESSENTIAL".
          print $0, "ESSENTIAL"
      } else {
          # If it does not exist, print the line and append "NOT ESSENTIAL".
          print $0, "NOT ESSENTIAL"
      }
  }
' "$ESSENTIALS_FILE" -
