# Liang et al. 2024 Analysis README

This document outlines the steps to process and analyze the gRNA data from the Liang et al. study.

## Step 1: Convert Supplementary Table S1 gRNAs to FASTA Format

First, we convert the relevant columns from the source CSV file into a separate FASTA file.

```bash

# Convert gRNA1 data

./csv_to_fasta.sh ../data/Liang/LiangMuller_May2025_Table1_Filtered-gRNAs.csv 1-2 4 lncrna > ../data/Liang/processed/Supplementary_tableS1_gRNAs.fasta


```


## Step 2: Run BLAT Alignments

Run `blat` to align the filtered gRNA sequences against the sequences given a prediction by our model for lncRNA `exon1` and `exon2`.

### Alignments for gRNA

```bash

# Align gRNA1 against exon 1 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon1.fasta ../data/Liang/processed/Supplementary_tableS1_gRNAs.fasta -minScore=15 -minIdentity=100 ../data/Liang/processed/tableS1_gRNAs_vs_exon1.psl

# Align gRNA1 against exon 2 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon2.fasta ../data/Liang/processed/Supplementary_tableS1_gRNAs.fasta -minScore=15 -minIdentity=100 ../data/Liang/processed/tableS1_gRNAs_vs_exon2.psl

```


## Step 3: Join blat matches with model predictions

```bash

# Join gRNAs with exon 1 model prediction

./join_blat_matches_Liang.sh ../data/Liang/processed/tableS1_gRNAs_vs_exon1.psl ../data/model_predictions/gencode-lncrna-ranking.csv ../data/Liang/LiangMuller_May2025-Table3-Gene-RRA-Ranking.csv > ../data/Liang/processed/annotated_tableS3_gRNAs_vs_exon1_prob.csv

# Join gRNAs with exon 2 model prediction

./join_blat_matches_Liang.sh ../data/Liang/processed/tableS1_gRNAs_vs_exon2.psl ../data/model_predictions/gencode-lncrna-ranking.csv ../data/Liang/LiangMuller_May2025-Table3-Gene-RRA-Ranking.csv  > ../data/Liang/processed/annotated_tableS3_gRNAs_vs_exon2_prob.csv

```


## Step 4: Visualize the distribution differences

Using the R script `essential_prob_distribution_Liang.R`, a violin plot is generated. This script joins exon1 and exon2 data and keeps 
the maximum probability for each gene id. Then computes a K-S stat between the four groups: Shared, Partially shared, Cell-type specific, 
and Non-essential.


## Step 5: Obtain a list of lncRNA sequences to test the model

In this step, we obtain a list of essential and non-essential lncRNA sequences. The corresponding exon1 and exon2 sequences are extracted. 
This list of sequences is later used in conjuction with the functionality prediction model to obtain a prediction of functionality. With these 
predictions, an effective measure of the model accuracy can be tested.

```bash

# Get a list of Ensembl ids

./find_lncRNA_guides.sh ../data/Liang/LiangMuller_May2025_Table1_Filtered-gRNAs.csv ../data/Liang/LiangMuller_May2025-Table3-Gene-RRA-Ranking.csv

# Obtain just unique Target Gene Id and ENSG Id combination hits

./filter_unique_hits.sh ../results/gRNA_lncRNA_matches.tsv

# Get the number of unique essential target gene ids

tail -n +2 ../results/gRNA_lncRNA_matches_unique_sorted.tsv | cut -f3 | sort -u | wc -l

# Add the highest probability assigned by the model by joining by ENSG Id

./add_probability_field.sh ../results/gRNA_lncRNA_matches_unique_sorted.tsv ../data/model_predictions/gencode-lncrna-ranking.csv

# Select the best hit (max probability) while keeping all gRNA IDs that match the same Target Gene ID and Ensembl Gene ID

./select_best_hit.sh ../results/gRNA_lncRNA_matches_with_prob.tsv

# TODO: Plot the resulting data using an R script.
# Lets use a similar script to essential_prob_distribution_Liang.R

# Filter just ENSG IDs without model probability assigned

./filter_na_prob.sh ../results/gRNA_lncRNA_matches_with_prob.tsv

# Obtain exon1 and exon2 sequences from ENSG IDs without probability assigned

./get_exon_sequences.sh ../results/gRNA_lncRNA_matches_unique_ensg_na_prob.tsv

# Get the number of unique essential target gene ids that matched at least one ENSG Id

tail -n +2 results/gRNA_lncRNA_matches_with_prob.tsv | awk -F'\t' '$9 != "NA"' | cut -f1 | sort -u | wc -l

# Get basic stats from the lncRNA annotations file

./analyze_gtf_stats.sh ../data/references/gencode.v49.long_noncoding_RNAs.gtf

```

