# LONG NCRNA
## Process histone marks for exon1
sh pol2r2a_processing.sh ../data/functional-lncrna-exon1-dataset.csv lncrna-exon1-pol2r2a-features
sh pol2r2a_processing.sh ../data/lncrna-exon1-negative-control-dataset.csv lncrna-exon1-NC-pol2r2a-features
## Process histone marks for exon2
sh pol2r2a_processing.sh ../data/functional-lncrna-exon2-dataset.csv lncrna-exon2-pol2r2a-features
sh pol2r2a_processing.sh ../data/lncrna-exon2-negative-control-dataset.csv lncrna-exon2-NC-pol2r2a-features
## Join datasets into one for R script
sh join_datasets.sh ../data/pol2r2a_feature/lncrna-exon1-pol2r2a-features.csv ../data/pol2r2a_feature/lncrna-exon2-pol2r2a-features.csv ../data/pol2r2a_feature/lncrna_positives_matrix
sh join_datasets.sh ../data/pol2r2a_feature/lncrna-exon1-NC-pol2r2a-features.csv ../data/pol2r2a_feature/lncrna-exon2-NC-pol2r2a-features.csv ../data/pol2r2a_feature/lncrna_negatives_matrix
sh join_datasets.sh ../data/pol2r2a_feature/lncrna_positives_matrix.csv ../data/pol2r2a_feature/lncrna_negatives_matrix.csv ../data/pol2r2a_feature/lncrna_matrix

#SHORT NCRNA
sh pol2r2a_processing.sh \
    ../data/functional-short-ncrna-dataset.csv \
    \
    short-ncrna-pol2r2a-features
sh pol2r2a_processing.sh \
    ../data/short-ncrna-negative-control-dataset.csv \
    \
    short-ncrna-NC-pol2r2a-features
sh join_datasets.sh \
    ../data/pol2r2a_feature/short-ncrna-pol2r2a-features.csv \
    ../data/pol2r2a_feature/short-ncrna-NC-pol2r2a-features.csv \
    ../data/pol2r2a_feature/short_ncrna_matrix


# PROTEIN CODING
sh pol2r2a_processing.sh ../data/functional-protein-exon2-dataset.csv protein-exon2-pol2r2a-features
sh pol2r2a_processing.sh ../data/protein-exon2-negative-control-dataset.csv protein-exon2-NC-pol2r2a-features

sh pol2r2a_processing.sh ../data/functional-protein-exon3-dataset.csv protein-exon3-pol2r2a-features
sh pol2r2a_processing.sh ../data/protein-exon3-negative-control-dataset.csv protein-exon3-NC-pol2r2a-features

sh join_datasets.sh ../data/pol2r2a_feature/protein-exon2-pol2r2a-features.csv ../data/pol2r2a_feature/protein-exon3-pol2r2a-features.csv ../data/pol2r2a_feature/protein_positives_matrix
sh join_datasets.sh ../data/pol2r2a_feature/protein-exon2-NC-pol2r2a-features.csv ../data/pol2r2a_feature/protein-exon3-NC-pol2r2a-features.csv ../data/pol2r2a_feature/protein_negatives_matrix
sh join_datasets.sh ../data/pol2r2a_feature/protein_positives_matrix.csv ../data/pol2r2a_feature/protein_negatives_matrix.csv ../data/pol2r2a_feature/protein_matrix
