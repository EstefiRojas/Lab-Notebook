# Montero et al. 2024 Analysis README

This document outlines the steps to process and analyze the gRNA data from the Montero et al. study.


## Step 1: Convert Supplementary Table 9 to FASTA Format

First, we convert the relevant columns from the source CSV file into two separate FASTA files, one for each gRNA.

```bash

# Convert gRNA1 data

./csv_to_fasta.sh ../data/Montero/Supplementary_table9_useThis.csv 2 4 gRNA1 > ../data/Montero/processed/Supplementary_table9_grna1.fasta

# Convert gRNA2 data

./csv_to_fasta.sh ../data/Montero/Supplementary_table9_useThis.csv 2 5 gRNA2 > ../data/Montero/processed/Supplementary_table9_grna2.fasta

```


## Step 2: Filter for lncRNA Records

Next, filter the newly created FASTA files to retain only the records that are associated with lncRNAs. This uses the `fasta_filter_Montero.sh` script.

```bash

# Filter gRNA1 file

./fasta_filter_Montero.sh ../data/Montero/processed/Supplementary_table9_grna1.fasta ../data/Montero/processed/Supplementary_table9_grna1_filtered.fasta lncrna

# Filter gRNA2 file

./fasta_filter_Montero.sh ../data/Montero/processed/Supplementary_table9_grna2.fasta ../data/Montero/processed/Supplementary_table9_grna2_filtered.fasta lncrna

```


## Step 3: Run BLAT Alignments

Run `blat` to align the filtered gRNA sequences against the sequences given a prediction by the model for `exon1` and `exon2`.

### Alignments for gRNA1

```bash

# Align gRNA1 against exon 1 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon1.fasta ../data/Montero/processed/Supplementary_table9_grna1_filtered.fasta -minScore=15 -minIdentity=100 ../data/Montero/processed/table9_grna1_vs_exon1.psl

# Align gRNA1 against exon 2 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon2.fasta ../data/Montero/processed/Supplementary_table9_grna1_filtered.fasta -minScore=15 -minIdentity=100 ../data/Montero/processed/table9_grna1_vs_exon2.psl

```

### Alignments for gRNA2

```bash

# Align gRNA2 against exon 1 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon1.fasta ../data/Montero/processed/Supplementary_table9_grna2_filtered.fasta -minScore=15 -minIdentity=100 ../data/Montero/processed/table9_grna2_vs_exon1.psl

# Align gRNA2 against exon 2 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon2.fasta ../data/Montero/processed/Supplementary_table9_grna2_filtered.fasta -minScore=15 -minIdentity=100 ../data/Montero/processed/table9_grna2_vs_exon2.psl

```


## Step 4: Finaly, join blat matches with model predictions

```bash

# Join gRNA1 with exon 1 model prediction

./join_blat_matches_Montero.sh ../data/Montero/processed/table9_grna1_vs_exon1.psl ../data/model_predictions/gencode-lncrna-ranking.csv ../data/Montero/Supplementary_table13_useThis.csv > ../data/Montero/processed/annotated_table9_grna1_vs_exon1_prob.csv

# Join gRNA2 with exon 1 model prediction

./join_blat_matches_Montero.sh ../data/Montero/processed/table9_grna2_vs_exon1.psl ../data/model_predictions/gencode-lncrna-ranking.csv ../data/Montero/Supplementary_table13_useThis.csv > ../data/Montero/processed/annotated_table9_grna2_vs_exon1_prob.csv

# Join gRNA1 with exon 2 model prediction

./join_blat_matches_Montero.sh ../data/Montero/processed/table9_grna1_vs_exon2.psl ../data/model_predictions/gencode-lncrna-ranking.csv > ../data/Montero/processed/table9_grna1_vs_exon2_prob.csv

# Join gRNA2 with exon 2 model prediction

./join_blat_matches_Montero.sh ../data/Montero/processed/table9_grna2_vs_exon2.psl ../data/model_predictions/gencode-lncrna-ranking.csv > ../data/Montero/processed/table9_grna2_vs_exon2_prob.csv

```


## Step 5: Visualize the distribution differences

Using the R script `essential_prob_distribution_Montero.R`, a violin plot is generated. This script joins exon1 and exon2 data and keeps 
the maximum probability for each gene id. Then computes a K-S stat between the two groups Essential and Non-essential.


## Step 6: Find lncRNA that align to gRNAs and compare against mRNA from UniProt



```bash

# Change names of essential labels to unify them with other studies.
./update_essential_labels.sh ../data/Montero/Supplementary_table13_useThis.csv ../data/Montero/processed/Supplementary_table13_new_labels.csv

# Align gRNA1 and gRNA2 to mRNA from uniprot
blat -t=dna -q=dna ../data/references/human_proteins_long.fasta ../data/Montero/processed/Supplementary_table9_grna1_filtered.fasta -minScore=15 -minIdentity=100 ../data/Montero/processed/table9_grna1_vs_mrna.psl -noHead
blat -t=dna -q=dna ../data/references/human_proteins_long.fasta ../data/Montero/processed/Supplementary_table9_grna2_filtered.fasta -minScore=15 -minIdentity=100 ../data/Montero/processed/table9_grna2_vs_mrna.psl -noHead

# Keep just antisense results
awk '$9 ~ /^-/' ../data/Montero/processed/table9_grna1_vs_mrna.psl > ../data/Montero/processed/table9_gRNA1_vs_mrna_antisense_only.psl
awk '$9 ~ /^-/' ../data/Montero/processed/table9_grna2_vs_mrna.psl > ../data/Montero/processed/table9_gRNA2_vs_mrna_antisense_only.psl


# Get a list of Ensembl ids. This script aligns gRNAs against the 
# full human genome (GRCh38) and then intersects the results with 
# a dedicated GENCODE lncRNA annotation file. The final report includes 
# antisense BLAT hits only; if a hit does not overlap a lncRNA, "NA" is reported.
#./find_lncRNA_guides.sh ../data/Montero/Table9_Filtered-gRNA1.csv ../data/Montero/Table3-Gene-RRA-Ranking.csv
blat -t=dna -q=dna ../data/references/gencode.v49.lncRNA_transcripts.fa ../data/Montero/processed/Supplementary_table9_grna1_filtered.fasta -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 ../tmp/montero_blat_grna1.psl.tmp -noHead
blat -t=dna -q=dna ../data/references/gencode.v49.lncRNA_transcripts.fa ../data/Montero/processed/Supplementary_table9_grna2_filtered.fasta -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 ../tmp/montero_blat_grna2.psl.tmp -noHead
# Filter to keep only lines where column 9 (strand) starts with '-' (reverse complement)
awk '$9 ~ /^-/' ../tmp/montero_blat_grna1.psl.tmp > ../tmp/montero_blat_grna1_antisense.psl
awk '$9 ~ /^-/' ../tmp/montero_blat_grna2.psl.tmp > ../tmp/montero_blat_grna2_antisense.psl

# Transform to bed
awk 'BEGIN{OFS="\t"} {print $14, $10, $9}' ../tmp/montero_blat_grna1_antisense.psl > ../tmp/montero_blat_grna1_antisense.bed
awk 'BEGIN{OFS="\t"} {print $14, $10, $9}' ../tmp/montero_blat_grna2_antisense.psl > ../tmp/montero_blat_grna2_antisense.bed

# Combine tables
touch ../tmp/montero_gRNA1_metadata_map.tsv

awk 'BEGIN {FS=","; OFS="\t"}
FNR==NR {
    sub(/\r$/, "", $8);
    if (FNR > 1) {
        essential_genes[$1] = $8
    }
    next
}
{
    sub(/\r$/, "", $4);
    # Check for the header row *first*
    if (FNR > 1) { 
        if($3 == "lncRNA_gRNA") {
            if ($1 in essential_genes) {
                # Add the sequence ($4) as the last column to the map
                print $1, essential_genes[$1], $2 "_gRNA1", $4 >> "../tmp/montero_gRNA1_metadata_map.tsv"
            } else {
                # Print non-essential gRNAs
                print $1, "non-essential", $2 "_gRNA1", $4 >> "../tmp/montero_gRNA1_metadata_map.tsv"
            }
        }
    }
}' ../data/Montero/processed/Supplementary_table13_new_labels.csv ../data/Montero/Supplementary_table9_useThis.csv


touch ../tmp/montero_gRNA2_metadata_map.tsv

awk 'BEGIN {FS=","; OFS="\t"}
FNR==NR {
    sub(/\r$/, "", $8);
    if (FNR > 1) {
        essential_genes[$1] = $8
    }
    next
}
{
    sub(/\r$/, "", $4);
    # Check for the header row *first*
    if (FNR > 1) { 
        if($3 == "lncRNA_gRNA") {
            if ($1 in essential_genes) {
                # Add the sequence ($4) as the last column to the map
                print $1, essential_genes[$1], $2 "_gRNA2", $4 >> "../tmp/montero_gRNA2_metadata_map.tsv"
            } else {
                # Print non-essential gRNAs
                print $1, "non-essential", $2 "_gRNA2", $4 >> "../tmp/montero_gRNA2_metadata_map.tsv"
            }
        }
    }
}' ../data/Montero/processed/Supplementary_table13_new_labels.csv ../data/Montero/Supplementary_table9_useThis.csv

# Obtain just unique Target Gene Id and ENSG Id combination hits.
echo -e "gRNA_ID\tEssentiality\tTarget_Gene_ID\tgRNA_Sequence\tlncRNA_ENSG_ID\tStrand" > ../results/gRNA1_lncRNA_antisense_matches_Montero.tsv

awk 'BEGIN{FS="\t"; OFS="\t"}
# 1. Read BLAT_OUTPUT_BED into memory (File 1)
FNR==NR {
    gRNA_id = $2;

    # Split the complex ID from $1 by the "|"
    n = split($1, parts, "|");
    ensg_id = parts[2];
    strand = parts[n];
    
    # Create a single string for this hit
    hit_string = ensg_id OFS strand
    
    # Use a temporary array to track unique gRNA_id/hit combinations
    # SUBSEP (Control-FS) is a guaranteed-safe separator
    unique_key = gRNA_id SUBSEP hit_string
    
    # Only add this hit if we havent seen this exact combo before
    if ( !(unique_key in _seen) ) {
        _seen[unique_key] = 1 # Mark as seen
        
        # Now, append this unique hit to the main blat_hits array
        # The value is a string of all hits, separated by newlines
        if (gRNA_id in blat_hits) {
            # Append new hit with a newline separator
            blat_hits[gRNA_id] = blat_hits[gRNA_id] "\n" hit_string
        } else {
            # This is the first hit for this gRNA_id
            blat_hits[gRNA_id] = hit_string
        }
    }
    next;
}

# 2. Process GRNA_METADATA_MAP (File 2 - the "left" file)
{
    gRNA_id = $3;
    essential = $2;
    target_gene = $1;
    sequence = $4;
    
    # Create the base output string (all the metadata)
    base_output = gRNA_id OFS essential OFS target_gene OFS sequence
    
    # Check if this gRNA_id was found in the BLAT file (fast O(1) lookup)
    if (gRNA_id in blat_hits) {
        # YES: It has hits. Get the string of all hits.
        all_hits_string = blat_hits[gRNA_id]
        
        # Split that string (by newline) into an array of individual hits
        n_hits = split(all_hits_string, hits_array, "\n")
        
        # Loop *only* through the k hits for this gRNA_id (very fast)
        for (i = 1; i <= n_hits; i++) {
            # hits_array[i] contains "ensg_id\tstrand"
            print base_output, hits_array[i]
        }
    } else {
        # NO: This gRNA_id had no BLAT hits. Print it once with placeholders.
        print base_output, "NA", "NA"
    }
}' ../tmp/montero_blat_grna1_antisense.bed ../tmp/montero_gRNA1_metadata_map.tsv >> ../results/gRNA1_lncRNA_antisense_matches_Montero.tsv


echo -e "gRNA_ID\tEssentiality\tTarget_Gene_ID\tgRNA_Sequence\tlncRNA_ENSG_ID\tStrand" > ../results/gRNA2_lncRNA_antisense_matches_Montero.tsv

awk 'BEGIN{FS="\t"; OFS="\t"}
# 1. Read BLAT_OUTPUT_BED into memory (File 1)
FNR==NR {
    gRNA_id = $2;

    # Split the complex ID from $1 by the "|"
    n = split($1, parts, "|");
    ensg_id = parts[2];
    strand = parts[n];
    
    # Create a single string for this hit
    hit_string = ensg_id OFS strand
    
    # Use a temporary array to track unique gRNA_id/hit combinations
    # SUBSEP (Control-FS) is a guaranteed-safe separator
    unique_key = gRNA_id SUBSEP hit_string
    
    # Only add this hit if we havent seen this exact combo before
    if ( !(unique_key in _seen) ) {
        _seen[unique_key] = 1 # Mark as seen
        
        # Now, append this unique hit to the main blat_hits array
        # The value is a string of all hits, separated by newlines
        if (gRNA_id in blat_hits) {
            # Append new hit with a newline separator
            blat_hits[gRNA_id] = blat_hits[gRNA_id] "\n" hit_string
        } else {
            # This is the first hit for this gRNA_id
            blat_hits[gRNA_id] = hit_string
        }
    }
    next;
}

# 2. Process GRNA_METADATA_MAP (File 2 - the "left" file)
{
    gRNA_id = $3;
    essential = $2;
    target_gene = $1;
    sequence = $4;
    
    # Create the base output string (all the metadata)
    base_output = gRNA_id OFS essential OFS target_gene OFS sequence
    
    # Check if this gRNA_id was found in the BLAT file (fast O(1) lookup)
    if (gRNA_id in blat_hits) {
        # YES: It has hits. Get the string of all hits.
        all_hits_string = blat_hits[gRNA_id]
        
        # Split that string (by newline) into an array of individual hits
        n_hits = split(all_hits_string, hits_array, "\n")
        
        # Loop *only* through the k hits for this gRNA_id (very fast)
        for (i = 1; i <= n_hits; i++) {
            # hits_array[i] contains "ensg_id\tstrand"
            print base_output, hits_array[i]
        }
    } else {
        # NO: This gRNA_id had no BLAT hits. Print it once with placeholders.
        print base_output, "NA", "NA"
    }
}' ../tmp/montero_blat_grna2_antisense.bed ../tmp/montero_gRNA2_metadata_map.tsv >> ../results/gRNA2_lncRNA_antisense_matches_Montero.tsv


./add_probability_field.sh ../results/gRNA1_lncRNA_antisense_matches_Montero.tsv ../data/model_predictions/gencode-lncrna-ranking.csv
./add_match_column_optimized.sh ../results/gRNA_lncRNA_matches_with_prob.tsv ../data/Montero/processed/table9_gRNA1_vs_mrna_antisense_only.psl ../results/gRNA1_essential_matches_optimized.tsv

./add_probability_field.sh ../results/gRNA2_lncRNA_antisense_matches_Montero.tsv ../data/model_predictions/gencode-lncrna-ranking.csv
./add_match_column_optimized.sh ../results/gRNA_lncRNA_matches_with_prob.tsv ../data/Montero/processed/table9_gRNA2_vs_mrna_antisense_only.psl ../results/gRNA2_essential_matches_optimized.tsv

# Compare both lists (lncRNA and mRNA), keep those without matches in mRNA sequences.
./add_match_column_optimized.sh ../results/gRNA_lncRNA_matches_best_hit_multi_Montero.tsv ../data/Montero/processed/table9_gRNA1_vs_mrna_antisense_only.psl ../results/gRNA1_essential_matches_optimized_Montero.tsv
./add_match_column_optimized.sh ../results/gRNA_lncRNA_matches_best_hit_multi_Montero.tsv ../data/Montero/processed/table9_gRNA2_vs_mrna_antisense_only.psl ../results/gRNA2_essential_matches_optimized_Montero.tsv

# Extract core essential matches
awk '$2 ~ /^core/' ../results/gRNA1_essential_matches_optimized.tsv > ../results/gRNA1_core_essential_matches.tsv
awk '$2 ~ /^core/' ../results/gRNA2_essential_matches_optimized.tsv > ../results/gRNA2_core_essential_matches.tsv
```
