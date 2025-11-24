# Liang et al. 2024 Analysis README

This document outlines the steps to process and analyze the gRNA data from the Liang et al. study.
The objective is to use the essentiality classification to see if there is a good correlation with
the model's predictions.

## Step 1: Convert Supplementary Table S1 gRNAs to FASTA Format

First, we convert the relevant columns from the source CSV file into a separate FASTA file.

```bash

# Convert gRNA1 data

./csv_to_fasta.sh ../data/Liang/LiangMuller_May2025_Table1_Filtered-gRNAs.csv 1-2 4 3 > ../data/Liang/processed/Supplementary_tableS1_gRNAs.fasta

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


# Better approach
The previous steps just considered lncRNAs that already have a probability assigned by our models at the date of writing for the blat against the gRNAs. In the following approach, we will run a blat match against all lncRNAs reported to date. This will generate an exhaustive list of all lncRNA to which the gRNAs attach to.

## Step 5: Obtain a list of lncRNA sequences to test the model

In this step, we obtain a list of essential and non-essential lncRNA sequences. The corresponding exon1 and exon2 sequences are extracted. 
This list of sequences is later used in conjuction with the functionality prediction model to obtain a prediction of functionality. With these 
predictions, an effective measure of the model accuracy can be tested.

```bash

# Get a list of Ensembl ids. This script aligns gRNAs against the 
# full human genome (GRCh38) and then intersects the results with 
# a dedicated GENCODE lncRNA annotation file. The final report includes 
# antisense BLAT hits only; if a hit does not overlap a lncRNA, "NA" is reported in 
# lncRNA fields.
./find_lncRNA_guides.sh ../data/Liang/LiangMuller_May2025_Table1_Filtered-gRNAs.csv ../data/Liang/LiangMuller_May2025-Table3-Gene-RRA-Ranking.csv


# Obtain just unique Target Gene Id and ENSG Id combination hits.
# This script processes the output from the main find_lncRNA_guides.sh
# script to create a filtered and sorted report. For each unique
# combination of a Target_Gene_ID, lncRNA_ENSG_ID, and gRNA_ID, it
# keeps a single entry. The final output is then sorted by the Target_Gene_ID.
./filter_unique_hits.sh ../results/gRNA_lncRNA_antisense_matches.tsv


# Get the number of unique essential target gene ids.
tail -n +2 ../results/gRNA_lncRNA_matches_unique_sorted.tsv | cut -f3 | sort -u | wc -l


# Add the highest probability assigned by the model by joining by ENSG Id.
# This script joins the main gRNA analysis report with a second
# file containing functional probability data. It matches records
# by the lncRNA Ensembl Gene ID, appends the "highest_prob"
# value, and reorders columns to prioritize the Target_Gene_ID.
./add_probability_field.sh ../results/gRNA_lncRNA_matches_unique_sorted.tsv ../data/model_predictions/gencode-lncrna-ranking.csv


# Select the best hit (max probability) while keeping all gRNA IDs that 
# match the same Target Gene ID and Ensembl Gene ID.
# This script filters an enriched gRNA report to find the best
# lncRNA hit for each unique Target_Gene_ID. If a gene has
# multiple lncRNA matches, it determines the best one by finding
# the highest probability. It then keeps ALL gRNA records that
# target that single best lncRNA. If a gene only has non-matches
# (NA), it keeps one of those records. The final output is sorted.
./select_best_hit.sh ../results/gRNA_lncRNA_matches_with_prob.tsv # This script creates '../results/gRNA_lncRNA_matches_best_hit_multi.tsv' file used downstream in mRNA comparisson.


# [OPTIONAL]
# Filter ENSG IDs without model probability assigned
./filter_na_prob.sh ../results/gRNA_lncRNA_matches_with_prob.tsv


# Obtain exon1 and exon2 sequences from ENSG IDs without probability assigned so that it can be run with our model
./get_exon_sequences.sh ../results/gRNA_lncRNA_matches_unique_ensg_na_prob.tsv


# Get the number of unique essential target gene ids that matched at least one ENSG Id
tail -n +2 ../results/gRNA_lncRNA_matches_with_prob.tsv | awk -F'\t' '$5 != "NA"' | cut -f1 | sort -u | wc -l


# Get basic stats from the lncRNA annotations file
./analyze_gtf_stats.sh ../data/references/gencode.v49.long_noncoding_RNAs.gtf

```

## Step 6: Get mRNA readings to cure the lncRNA list of essentials

This step involves finding protein coding genes that will attach to the list of gRNAs from the previous step.

```bash

# Obtain mRNA readings
# Download uniprot data for human proteins
curl -X GET "https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606)&format=fasta&compressed=true" \
  -H "Accept: application/gzip" \
  -o uniprotkb_AND_model_organism_9606_AND_r_2025_10_29.fasta.gz

# Index the genome
miniprot -t8 -d GRCh38.genome.mpi GRCh38.primary_assembly.genome.fa

# Align the protein sequences from uniprot with the genome
miniprot --aln -N 0 GRCh38.genome.mpi uniprotkb_AND_model_organism_9606_AND_r_2025_10_29.fasta.gz > proteome-human-proteins.miniprot.aln

# Convert the alignment to fasta format
python3 parse_aln_to_fasta.py proteome-human-proteins.miniprot.aln GRCh38.primary_assembly.genome.fa human_proteins.fasta

# Run the blat command
blat -t=dna -q=dna ../data/references/human_proteins.fasta ../data/Liang/processed/Supplementary_tableS1_gRNAs.fasta -minScore=15 -minIdentity=100 ../data/Liang/processed/tableS1_gRNAs_vs_mRNA_no_introns.psl -noHead

# Keep just antisense results
awk '$9 ~ /^-/' ../data/Liang/processed/tableS1_gRNAs_vs_mRNA_no_introns.psl > ../data/Liang/processed/tableS1_gRNAs_vs_mrna_no_introns_antisense_only.psl

# Compare both lists (lncRNA and mRNA), mark those that matched both lncRNA and mRNA with the same gRNA in the same strand.
./add_match_column_optimized.sh ../results/gRNA_lncRNA_matches_best_hit_multi.tsv ../data/Liang/processed/tableS1_gRNAs_vs_mrna_no_introns_antisense_only.psl ../results/gRNA_essential_matches_no_introns_optimized.tsv

# Add the lncRNA gene name
./map_lncRNA_genes.sh ../results/gRNA_essential_matches_no_introns_optimized.tsv ../data/references/gencode.v49.long_noncoding_RNAs.gtf
```
## Step 7: Plot the resulting data using an R script.

This step uses the script `Liang_essentials_boxplot.R`.


