#!/opt/local/bin/bash
# Helper script to run all necessary scripts to get the max signal peak for 
# rna types for a specific histone.

set -uex

# Verify bash version
echo `bash --version`
# Check for correct usage
if [ $# -ne 1 ]; then
    echo "Usage: $0 histone_name" # biosample_class"
    exit 1
fi

# Store input file
HISTONE_NAME=$1
#BIOSAMPLE_CLASS=$2


# LONG NCRNA
## Process histone marks for exon1
bash histone_processing.sh ../data/datasets/functional-lncrna-exon1-dataset.csv "$HISTONE_NAME" lncrna-exon1-histone-feature
bash histone_processing.sh ../data/datasets/lncrna-exon1-negative-control-dataset.csv "$HISTONE_NAME" lncrna-exon1-NC-histone-feature
## Process histone marks for exon2
bash histone_processing.sh ../data/datasets/functional-lncrna-exon2-dataset.csv "$HISTONE_NAME" lncrna-exon2-histone-feature
bash histone_processing.sh ../data/datasets/lncrna-exon2-negative-control-dataset.csv "$HISTONE_NAME" lncrna-exon2-NC-histone-feature


#SHORT NCRNA
bash histone_processing.sh \
    ../data/datasets/functional-short-ncrna-dataset.csv \
    "$HISTONE_NAME" \
    short-ncrna-histone-feature
bash histone_processing.sh \
    ../data/datasets/short-ncrna-negative-control-dataset.csv \
    "$HISTONE_NAME" \
    short-ncrna-NC-histone-feature


# PROTEIN CODING
bash histone_processing.sh ../data/datasets/functional-protein-exon2-dataset.csv "$HISTONE_NAME" protein-exon2-histone-feature
bash histone_processing.sh ../data/datasets/protein-exon2-negative-control-dataset.csv "$HISTONE_NAME" protein-exon2-NC-histone-feature

bash histone_processing.sh ../data/datasets/functional-protein-exon3-dataset.csv "$HISTONE_NAME" protein-exon3-histone-feature
bash histone_processing.sh ../data/datasets/protein-exon3-negative-control-dataset.csv "$HISTONE_NAME" protein-exon3-NC-histone-feature
