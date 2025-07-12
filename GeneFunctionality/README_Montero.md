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

