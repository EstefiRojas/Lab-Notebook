# Liu et al. 2017 Analysis README

This document outlines the steps to process and analyze the gRNA data from the Liu et al. study.

The same analyses as with the other studies were conducted with this study data. Differently to the other studies, no blat match was 
found with the functionality predictions as of 7 Oct 2025. Thus we directly describe the steps taken to extract Ensbl Ids from the 
reported gRNAs to compute a functional probability with our model.

## Step 1: Obtain a list of lncRNA sequences to test the model

In this step, we obtain a list of essential and non-essential lncRNA sequences. The corresponding exon1 and exon2 sequences are extracted. 
This list of sequences is later used in conjuction with the functionality prediction model to obtain a prediction of functionality. With these 
predictions, an effective measure of the model accuracy can be tested.

```bash

# Get a list of Ensembl Ids 

./find_lncRNA_guides_Liu.sh TableS2_gRNA_sequences.csv

# Count the number of unique essential LH ids

awk -F',' 'NR > 1 && tolower($1) != "negative_control" {ids[$1]} END {print length(ids)}' ../data/Liu/TableS2_gRNA_sequences.csv
#16401

# Count the number of unique LH ids that matched the genome

awk -F'\t' 'NR > 1 && tolower($1) != "negative_control" {ids[$2]} END {print length(ids)}' ../results/gRNA_lncRNA_matches_Liu.tsv
#16401

# Get unique gene ids and ensembl id


```
