# Zhu et al. 2016 Analysis README

This document outlines the steps to process and analyze the gRNA data from the Zhu et al. study.

## Step 1: Convert Supplementary Table 9 to FASTA Format

First, we convert the relevant columns from the source CSV file into two separate FASTA files, one for each gRNA.

```bash

# Convert gRNA1 data

./csv_to_fasta.sh ../data/Zhu/Supp_Table10.csv 2,15 8 sgRNA1 > ../data/Zhu/processed/Suppl_Table10_sgRNA1.fasta

# Convert gRNA2 data

./csv_to_fasta.sh ../data/Zhu/Supp_Table10.csv 2,15 13 sgRNA2 > ../data/Zhu/processed/Suppl_Table10_sgRNA2.fasta

```


## Step 2: Run BLAT Alignments

Run `blat` to align the filtered gRNA sequences against the sequences given a prediction by the model for `exon1` and `exon2`.

### Alignments for sgRNA1

```bash

# Align gRNA1 against exon 1 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon1.fasta ../data/Zhu/processed/Suppl_Table10_sgRNA1.fasta -minScore=15 -minIdentity=100 ../data/Zhu/processed/table10_grna1_vs_exon1.psl

# Align gRNA1 against exon 2 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon2.fasta ../data/Zhu/processed/Suppl_Table10_sgRNA1.fasta -minScore=15 -minIdentity=100 ../data/Zhu/processed/table10_grna1_vs_exon2.psl

```

### Alignments for sgRNA2

```bash

# Align gRNA2 against exon 1 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon1.fasta ../data/Zhu/processed/Suppl_Table10_sgRNA2.fasta -minScore=15 -minIdentity=100 ../data/Zhu/processed/table10_grna2_vs_exon1.psl

# Align gRNA2 against exon 2 model prediction sequences

blat -t=dna -q=dna ../data/model_predictions/lncrna_exon2.fasta ../data/Zhu/processed/Suppl_Table10_sgRNA2.fasta -minScore=15 -minIdentity=100 ../data/Zhu/processed/table10_grna2_vs_exon2.psl

```

This run will reveal very little matches against the lncRNA sequences. Just HOTAIR matches. As consequence, no further analysis is possible.


