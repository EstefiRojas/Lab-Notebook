# LONG NCRNA
## Process histone marks for exon1
sh dinucleotide_frequencies.sh ../data/datasets/functional-lncrna-exon1-dataset.csv lncrna-exon1-dinucleotide-feature
sh dinucleotide_frequencies.sh ../data/datasets/lncrna-exon1-negative-control-dataset.csv lncrna-exon1-NC-dinucleotide-feature
## Process histone marks for exon2
sh dinucleotide_frequencies.sh ../data/datasets/functional-lncrna-exon2-dataset.csv lncrna-exon2-dinucleotide-feature
sh dinucleotide_frequencies.sh ../data/datasets/lncrna-exon2-negative-control-dataset.csv lncrna-exon2-NC-dinucleotide-feature
## Join datasets into one for R script
#sh join_datasets.sh ../data/dinucleotide_feature/lncrna-exon1-dinucleotide-feature.csv ../data/dinucleotide_feature/lncrna-exon2-dinucleotide-feature.csv ../data/dinucleotide_feature/lncrna_positives_matrix
#sh join_datasets.sh ../data/dinucleotide_feature/lncrna-exon1-NC-dinucleotide-feature.csv ../data/dinucleotide_feature/lncrna-exon2-NC-dinucleotide-feature.csv ../data/dinucleotide_feature/lncrna_negatives_matrix
#sh join_datasets.sh ../data/dinucleotide_feature/lncrna_positives_matrix.csv ../data/dinucleotide_feature/lncrna_negatives_matrix.csv ../data/dinucleotide_feature/lncrna_matrix

#SHORT NCRNA
sh dinucleotide_frequencies.sh \
    ../data/datasets/functional-short-ncrna-dataset.csv \
    \
    short-ncrna-dinucleotide-feature
sh dinucleotide_frequencies.sh \
    ../data/datasets/short-ncrna-negative-control-dataset.csv \
    \
    short-ncrna-NC-dinucleotide-feature
#sh join_datasets.sh \
    ../data/dinucleotide_feature/short-ncrna-dinucleotide-feature.csv \
    ../data/dinucleotide_feature/short-ncrna-NC-dinucleotide-feature.csv \
    ../data/dinucleotide_feature/short_ncrna_matrix


# PROTEIN CODING
sh dinucleotide_frequencies.sh ../data/datasets/functional-protein-exon2-dataset.csv protein-exon2-dinucleotide-feature
sh dinucleotide_frequencies.sh ../data/datasets/protein-exon2-negative-control-dataset.csv protein-exon2-NC-dinucleotide-feature

sh dinucleotide_frequencies.sh ../data/datasets/functional-protein-exon3-dataset.csv protein-exon3-dinucleotide-feature
sh dinucleotide_frequencies.sh ../data/datasets/protein-exon3-negative-control-dataset.csv protein-exon3-NC-dinucleotide-feature

#sh join_datasets.sh ../data/dinucleotide_feature/protein-exon2-dinucleotide-feature.csv ../data/dinucleotide_feature/protein-exon3-dinucleotide-feature.csv ../data/dinucleotide_feature/protein_positives_matrix
#sh join_datasets.sh ../data/dinucleotide_feature/protein-exon2-NC-dinucleotide-feature.csv ../data/dinucleotide_feature/protein-exon3-NC-dinucleotide-feature.csv ../data/dinucleotide_feature/protein_negatives_matrix
#sh join_datasets.sh ../data/dinucleotide_feature/protein_positives_matrix.csv ../data/dinucleotide_feature/protein_negatives_matrix.csv ../data/dinucleotide_feature/protein_matrix
