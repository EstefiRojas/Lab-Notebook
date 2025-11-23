#!/bin/bash

# Usage: ./join_unified_matches.sh <psl_file> <model_predictions_csv>

PSL_FILE="$1"
MODEL_FILE="$2"

if [[ -z "$PSL_FILE" || -z "$MODEL_FILE" ]]; then
    echo "Usage: $0 <psl_file> <model_predictions_csv>"
    exit 1
fi

# Header
echo "Study,gRNA_ID,Target_Gene_ID,gRNA_Type,Essentiality,GeneID,TranscriptID,ExonNum,Probability_Functional"

# Process PSL and Join
awk '
  BEGIN { OFS="," }
  FNR > 5 {
      # Parse PSL
      q_name = $10
      t_name = $14
      
      # Parse Query Name (Study|gRNA_ID|Target_Gene_ID|gRNA_Type|Essentiality)
      split(q_name, q_parts, "|")
      study = q_parts[1]
      grna_id = q_parts[2]
      target_gene = q_parts[3]
      type = q_parts[4]
      essentiality = q_parts[5]
      
      # Parse Target Name (GeneID/TranscriptID_ExonNum)
      # Format: ENSG.../ENST..._exon1
      split(t_name, t_parts, "/")
      gene_id = t_parts[1]
      
      split(t_parts[2], t_subparts, "_")
      transcript_id = t_subparts[1]
      exon_num = t_subparts[2]
      
      # Print intermediate format for joining: key, rest
      # Key = GeneID SUBSEP TranscriptID
      print gene_id, transcript_id, study, grna_id, target_gene, type, essentiality, exon_num
  }' "$PSL_FILE" | \
awk -F, -v OFS=, '
  BEGIN { sep = SUBSEP }
  # Read Model Predictions (File 2 in pipeline, but passed as argument to script so read first by awk logic if we structure it right)
  NR==FNR {
      if (FNR > 1) {
          # Model File: Col 1 = GeneID, Col 2 = TransID1|TransID2...
          # We need to map GeneID + TransID -> Prob
          gene_id = $1
          split($2, trans_list, "|")
          split($4, prob_list, "|")
          
          for (i in trans_list) {
              key = gene_id sep trans_list[i]
              val = prob_list[i]
              data_map[key] = val
          }
      }
      next
  }
  {
      # Process piped input (PSL parsed)
      # Input: gene_id transcript_id study grna_id target_gene type essentiality exon_num
      gene_id = $1
      transcript_id = $2
      study = $3
      grna_id = $4
      target_gene = $5
      type = $6
      essentiality = $7
      exon_num = $8
      
      key = gene_id sep transcript_id
      
      prob = "NA"
      if (key in data_map) {
          prob = data_map[key]
      }
      
      print study, grna_id, target_gene, type, essentiality, gene_id, transcript_id, exon_num, prob
  }
' "$MODEL_FILE" -
