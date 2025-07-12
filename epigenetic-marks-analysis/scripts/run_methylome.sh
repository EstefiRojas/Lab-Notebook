#!/bin/bash

set -uex


# LONG NCRNA
## Process methylome marks for exon1
sh methylome_processing.sh ../data/datasets/functional-lncrna-exon1-dataset.csv lncrna-exon1-methylome-feature
sh methylome_processing.sh ../data/datasets/lncrna-exon1-negative-control-dataset.csv lncrna-exon1-NC-methylome-feature
## Process methylome marks for exon2
sh methylome_processing.sh ../data/datasets/functional-lncrna-exon2-dataset.csv lncrna-exon2-methylome-feature
sh methylome_processing.sh ../data/datasets/lncrna-exon2-negative-control-dataset.csv lncrna-exon2-NC-methylome-feature
## Join datasets into one for R script
#sh join_datasets.sh ../data/datasets/methylome_feature/lncrna-exon1-methylome-feature.csv ../data/datasets/methylome_feature/lncrna-exon2-methylome-feature.csv ../data/datasets/methylome_feature/lncrna_positives_matrix
#sh join_datasets.sh ../data/datasets/methylome_feature/lncrna-exon1-NC-methylome-feature.csv ../data/datasets/methylome_feature/lncrna-exon2-NC-methylome-feature.csv ../data/datasets/methylome_feature/lncrna_negatives_matrix
#sh join_datasets.sh ../data/datasets/methylome_feature/lncrna_positives_matrix.csv ../data/datasets/methylome_feature/lncrna_negatives_matrix.csv ../data/datasets/methylome_feature/lncrna_matrix

#SHORT NCRNA
sh methylome_processing.sh \
    ../data/datasets/functional-short-ncrna-dataset.csv \
    \
    short-ncrna-methylome-feature
sh methylome_processing.sh \
    ../data/datasets/short-ncrna-negative-control-dataset.csv \
    \
    short-ncrna-NC-methylome-feature
#sh join_datasets.sh \
#    ../data/datasets/methylome_feature/short-ncrna-methylome-feature.csv \
#    ../data/datasets/methylome_feature/short-ncrna-NC-methylome-feature.csv \
#    ../data/datasets/methylome_feature/short_ncrna_matrix


# PROTEIN CODING
sh methylome_processing.sh ../data/datasets/functional-protein-exon2-dataset.csv protein-exon2-methylome-feature
sh methylome_processing.sh ../data/datasets/protein-exon2-negative-control-dataset.csv protein-exon2-NC-methylome-feature

sh methylome_processing.sh ../data/datasets/functional-protein-exon3-dataset.csv protein-exon3-methylome-feature
sh methylome_processing.sh ../data/datasets/protein-exon3-negative-control-dataset.csv protein-exon3-NC-methylome-feature

#sh join_datasets.sh ../data/datasets/methylome_feature/protein-exon2-methylome-feature.csv ../data/datasets/methylome_feature/protein-exon3-methylome-feature.csv ../data/datasets/methylome_feature/protein_positives_matrix
#sh join_datasets.sh ../data/datasets/methylome_feature/protein-exon2-NC-methylome-feature.csv ../data/datasets/methylome_feature/protein-exon3-NC-methylome-feature.csv ../data/datasets/methylome_feature/protein_negatives_matrix
#sh join_datasets.sh ../data/datasets/methylome_feature/protein_positives_matrix.csv ../data/datasets/methylome_feature/protein_negatives_matrix.csv ../data/datasets/methylome_feature/protein_matrix

