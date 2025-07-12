#!/opt/local/bin/bash
# Helper script to run all necessary scripts to get the max signal peak for 
# for all histones of interest.
bash histone_processing.sh ../results/short-ncrna-negative-control-dataset-features.csv H3K27ac lncrna-exon1-NC-histone-feature
bash histone_processing.sh ../results/short-ncrna-negative-control-dataset-features.csv H3K36me3 lncrna-exon1-NC-histone-feature
bash histone_processing.sh ../results/short-ncrna-negative-control-dataset-features.csv H3K79me2 lncrna-exon1-NC-histone-feature
bash chrm_acc_processing.sh ../data/datasets/lncrna-exon2-negative-control-dataset.csv  lncrna-exon2-NC-chrm_acc-feature



bash histone_processing.sh \
    ../results/short-ncrna-negative-control-dataset-features.csv \
    H3K27ac \
    short-ncrna-NC-histone-feature

bash histone_processing.sh \
    ../results/short-ncrna-negative-control-dataset-features.csv \
    H3K36me3 \
    short-ncrna-NC-histone-feature

bash histone_processing.sh \
    ../results/short-ncrna-negative-control-dataset-features.csv \
    H3K79me2 \
    short-ncrna-NC-histone-feature

bash chrm_acc_processing.sh \
    ../results/short-ncrna-negative-control-dataset-features.csv \
     \
    short-ncrna-NC-chrm_acc-feature