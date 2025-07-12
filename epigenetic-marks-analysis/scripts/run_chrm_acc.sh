#!/opt/local/bin/bash

set -uex

#EXP_TYPE=$1
#BIO_CLASS=$2

# LONG NCRNA
## Process chrm_acc marks for exon1
sh chrm_acc_processing.sh ../data/datasets/functional-lncrna-exon1-dataset.csv lncrna-exon1-chrm_acc-feature
sh chrm_acc_processing.sh ../data/datasets/lncrna-exon1-negative-control-dataset.csv lncrna-exon1-NC-chrm_acc-feature
## Process chrm_acc marks for exon2
sh chrm_acc_processing.sh ../data/datasets/functional-lncrna-exon2-dataset.csv lncrna-exon2-chrm_acc-feature
sh chrm_acc_processing.sh ../data/datasets/lncrna-exon2-negative-control-dataset.csv  lncrna-exon2-NC-chrm_acc-feature
## Join datasets into one for R script
#sh join_datasets.sh ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_lncrna-exon1-chrm_acc-feature.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_lncrna-exon2-chrm_acc-feature.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_lncrna_positives_matrix
#sh join_datasets.sh ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_lncrna-exon1-NC-chrm_acc-feature.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_lncrna-exon2-NC-chrm_acc-feature.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_lncrna_negatives_matrix
#sh join_datasets.sh ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_lncrna_positives_matrix.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_lncrna_negatives_matrix.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_lncrna_matrix

#SHORT NCRNA
sh chrm_acc_processing.sh \
    ../data/datasets/functional-short-ncrna-dataset.csv \
     \
    short-ncrna-chrm_acc-feature
sh chrm_acc_processing.sh \
    ../data/datasets/short-ncrna-negative-control-dataset.csv \
     \
    short-ncrna-NC-chrm_acc-feature
#sh join_datasets.sh \
#    ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_short-ncrna-chrm_acc-feature.csv \
#    ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_short-ncrna-NC-chrm_acc-feature.csv \
#    ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_short_ncrna_matrix


# PROTEIN CODING
sh chrm_acc_processing.sh ../data/datasets/functional-protein-exon2-dataset.csv  protein-exon2-chrm_acc-feature
sh chrm_acc_processing.sh ../data/datasets/protein-exon2-negative-control-dataset.csv  protein-exon2-NC-chrm_acc-feature

sh chrm_acc_processing.sh ../data/datasets/functional-protein-exon3-dataset.csv  protein-exon3-chrm_acc-feature
sh chrm_acc_processing.sh ../data/datasets/protein-exon3-negative-control-dataset.csv  protein-exon3-NC-chrm_acc-feature

#sh join_datasets.sh ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_protein-exon2-chrm_acc-feature.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_protein-exon3-chrm_acc-feature.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_protein_positives_matrix
#sh join_datasets.sh ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_protein-exon2-NC-chrm_acc-feature.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_protein-exon3-NC-chrm_acc-feature.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_protein_negatives_matrix
#sh join_datasets.sh ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_protein_positives_matrix.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_protein_negatives_matrix.csv ../data/datasets/chrm_acc_feature/broadPeak/chrm_acc_protein_matrix

