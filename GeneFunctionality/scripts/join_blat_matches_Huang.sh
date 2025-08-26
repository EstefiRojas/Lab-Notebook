#!/bin/bash

# Bash script to perform a two-stage join and a final annotation step.
# MODIFIED to handle new input file formats and a two-column essentiality output.
# 1. Processes a PSL file to filter and reformat columns. It now parses
#    the Q_name to extract an ENSEMBL gene ID (e.g., ENSG...).
# 2. Joins the result with a CSV data file based on a composite key (GeneID, TranscriptID),
#    appending the 'Probability_Functional' column.
# 3. Annotates the result by looking up the parsed ENSG ID in a third file.
#    - If found, it's marked "Essential" and the coPARSE value is appended.
#    - Otherwise, it's marked "Non-essential" and "NA" is appended.
#
# Usage: ./modified_join_script.sh <original.psl> <data.csv> <essentials.csv> > final_output.csv

# --- Input Files ---
ORIGINAL_PSL_FILE="$1"  # The first input file (PSL format)
DATA_FILE="$2"          # The second input file (CSV format with functional probability)
ESSENTIALS_FILE="$3"    # The third file for annotation (CSV format with essentiality groups)

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
# **MODIFIED LOGIC**: Header for the two new annotation columns.
HEADER_ESSENTIAL="Essentiality_Status,coPARSE-lncRNA"
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
          
          # This logic can be kept for context, but is not used in joins.
          q_name_prefix = "NA"
          if (split(q_name_original, qname_parts, "_") >= 2) {
              q_name_prefix = qname_parts[1]
          }

          # Parse the Q_name (e.g., lncRNA_ENSG..._crRNA_A)
          # to extract the ENSG identifier for the final annotation step.
          q_name_suffix = "NA" # Default value
          if (split(q_name_original, qname_parts_ensg, "_") >= 2) {
              q_name_suffix = qname_parts_ensg[2]
          }

          t_name_original = $14
          gene_id = transcript_id = exon_num = "NA"
          # This parsing logic correctly handles the "GeneID/TranscriptID_ExonNum" format.
          if (split(t_name_original, parts_slash, "/") >= 2) {
              gene_id = parts_slash[1]
              if (split(parts_slash[2], parts_underscore, "_") >= 2) {
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
          
          # Print the last three fields (blockSizes, qStarts, tStarts), which contain
          # commas, by wrapping them in double quotes to create valid CSV output.
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
# Third stage: Annotate the result using the new essentials file format.
awk -F, -v OFS=, '
  # First, read the essentials file and build a map of Gene -> coPARSE value.
  # NR==FNR is true only while reading the first file (ESSENTIALS_FILE).
  NR == FNR {
      if (FNR > 1) { # Skip the header row.
          # The ENSG ID is in the 5th column.
          gene_id = $5
          gsub(/\r$/, "", gene_id)

          # The coPARSE value is in the 6th column.
          coparse_val = $6
          gsub(/\r$/, "", coparse_val)

          # Check if the gene is present in any of the first four cell line columns.
          status = ""
          if ($1 != "") { status = "Essential" }
          if ($2 != "") { status = "Essential" }
          if ($3 != "") { status = "Essential" }
          if ($4 != "") { status = "Essential" }

          # Only add to map if it was found in at least one cell line.
          if (status == "Essential" && gene_id != "") {
              essentials_map[gene_id] = coparse_val
          }
      }
      next # Move to the next line of ESSENTIALS_FILE.
  }

  # This block runs for the piped input (output from the second join).
  {
      # The join key is in the 11th column ("Q_name_suffix"), which now contains the parsed ENSG ID.
      key = $11
      gsub(/\r$/, "", key) # Clean the key.
      gsub(/\r$/, "", $0) # Clean the full line from the pipe.

      # **MODIFIED LOGIC**: Check if the key from our data exists in the map of essential IDs.
      if (key in essentials_map) {
          # If it exists, print the line, "Essential", and the stored coPARSE value.
          print $0, "Essential", essentials_map[key]
      } else {
          # If it does not exist, print the line, "Non-essential", and "NA".
          print $0, "Non-essential", "NA"
      }
  }
' "$ESSENTIALS_FILE" -
