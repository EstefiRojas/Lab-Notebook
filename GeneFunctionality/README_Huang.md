# Huang et al. 2024 Analysis README

This document outlines the steps to process and analyze the crRNA data from the Huang et al. study.

## Step 1: Convert Supplementary Table 6b crRNAs to FASTA Format

First, we convert the relevant columns from the source CSV file into a separate FASTA file.

```bash

# Convert crRNA data

./csv_to_fasta.sh ../data/Huang/SuppT6b.csv 5,7 3 crRNA_A > ../data/Huang/processed/crRNA_A_seq.fasta


```


## Step 2: Run BLAT Alignments

Run `blat` to align the filtered crRNA sequences against the sequences given a prediction by our model for lncRNA `exon1` and `exon2`.

### Alignments for crRNA

```bash

# Align crRNA against exon 1 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon1.fasta ../data/Huang/processed/crRNA_A_seq.fasta -minScore=15 -minIdentity=100 ../data/Huang/processed/crRNAs_vs_exon1.psl

# Align crRNA against exon 2 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon2.fasta ../data/Huang/processed/crRNA_A_seq.fasta -minScore=15 -minIdentity=100 ../data/Huang/processed/crRNAs_vs_exon2.psl

```


## Step 3: Finaly, join blat matches with model predictions

```bash

# Join ENSEMBLE Gene IDs for essential genes with corresponding coPARSE data

awk -F, -v OFS=, 'NR==FNR{if(FNR>1){sub(/\r$/,"",$2); map[$1]=$2}; next} {sub(/\r$/,"",$0); if(FNR==1){print $0,"coPARSE-lncRNA"; next} val="NA"; for(i=1;i<=4;i++) if($i in map){val=map[$i]; break}; print $0, val}' ../data/Huang/essential_genes_coparse.csv ../data/Huang/column_matches_ensgids.csv > ../data/Huang/processed/essential_genes_ensgids.csv

# Join gRNAs with exon 1 model prediction

./join_blat_matches_Huang.sh ../data/Huang/processed/crRNAs_vs_exon1.psl ../data/model_predictions/gencode-lncrna-ranking.csv ../data/Huang/processed/essential_genes_ensgids.csv > ../data/Huang/processed/annotated_crRNAs_vs_exon1_prob.csv

# Join gRNAs with exon 2 model prediction

./join_blat_matches_Huang.sh ../data/Huang/processed/crRNAs_vs_exon2.psl ../data/model_predictions/gencode-lncrna-ranking.csv ../data/Huang/processed/essential_genes_ensgids.csv > ../data/Huang/processed/annotated_crRNAs_vs_exon2_prob.csv

```

## Step 4: Visualize the distribution differences

Using the R script `essential_prob_distribution_Huang.R`, a box plot is generated. This script joins exon1 and exon2 data and keeps 
the maximum probability for each gene id. Then computes a K-S stat between the three groups: coParse, Non coParse, and Non-essential.

