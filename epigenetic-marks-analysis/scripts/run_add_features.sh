# ADD Chromatin Accessibility
./add_feature_from_to.sh ../data/datasets/chrm_acc_feature/chrm_acc_lncrna-exon1-chrm_acc-feature.csv ../data/features/functional-lncrna-exon1-dataset-features.csv 1 && \
./add_feature_from_to.sh ../data/datasets/chrm_acc_feature/chrm_acc_lncrna-exon2-chrm_acc-feature.csv ../data/features/functional-lncrna-exon2-dataset-features.csv 1 && \
./add_feature_from_to.sh ../data/datasets/chrm_acc_feature/chrm_acc_protein-exon2-chrm_acc-feature.csv ../data/features/functional-protein-exon2-dataset-features.csv 1 && \
./add_feature_from_to.sh ../data/datasets/chrm_acc_feature/chrm_acc_protein-exon3-chrm_acc-feature.csv ../data/features/functional-protein-exon3-dataset-features.csv 1 && \
./add_feature_from_to.sh ../data/datasets/chrm_acc_feature/chrm_acc_short-ncrna-chrm_acc-feature.csv ../data/features/functional-short-ncrna-dataset-features.csv 1 && \

./add_feature_from_to.sh ../data/datasets/chrm_acc_feature/chrm_acc_lncrna-exon1-NC-chrm_acc-feature.csv ../data/features/lncrna-exon1-negative-control-dataset-features.csv 1 && \
./add_feature_from_to.sh ../data/datasets/chrm_acc_feature/chrm_acc_lncrna-exon2-NC-chrm_acc-feature.csv ../data/features/lncrna-exon2-negative-control-dataset-features.csv 1 && \
./add_feature_from_to.sh ../data/datasets/chrm_acc_feature/chrm_acc_protein-exon2-NC-chrm_acc-feature.csv ../data/features/protein-exon2-negative-control-dataset-features.csv 1 && \
./add_feature_from_to.sh ../data/datasets/chrm_acc_feature/chrm_acc_protein-exon3-NC-chrm_acc-feature.csv ../data/features/protein-exon3-negative-control-dataset-features.csv 1 && \
./add_feature_from_to.sh ../data/datasets/chrm_acc_feature/chrm_acc_short-ncrna-NC-chrm_acc-feature.csv ../data/features/short-ncrna-negative-control-dataset-features.csv 1

# ADD Histone H3K27ac
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K27ac/H3K27ac_lncrna-exon1-histone-feature.csv ../data/features/functional-lncrna-exon1-dataset-features-1.csv 2 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K27ac/H3K27ac_lncrna-exon2-histone-feature.csv ../data/features/functional-lncrna-exon2-dataset-features-1.csv 2 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K27ac/H3K27ac_protein-exon2-histone-feature.csv ../data/features/functional-protein-exon2-dataset-features-1.csv 2 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K27ac/H3K27ac_protein-exon3-histone-feature.csv ../data/features/functional-protein-exon3-dataset-features-1.csv 2 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K27ac/H3K27ac_short-ncrna-histone-feature.csv ../data/features/functional-short-ncrna-dataset-features-1.csv 2 && \

./add_feature_from_to.sh ../data/datasets/histone_feature/H3K27ac/H3K27ac_lncrna-exon1-NC-histone-feature.csv ../data/features/lncrna-exon1-negative-control-dataset-features-1.csv 2 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K27ac/H3K27ac_lncrna-exon2-NC-histone-feature.csv ../data/features/lncrna-exon2-negative-control-dataset-features-1.csv 2 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K27ac/H3K27ac_protein-exon2-NC-histone-feature.csv ../data/features/protein-exon2-negative-control-dataset-features-1.csv 2 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K27ac/H3K27ac_protein-exon3-NC-histone-feature.csv ../data/features/protein-exon3-negative-control-dataset-features-1.csv 2 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K27ac/H3K27ac_short-ncrna-NC-histone-feature.csv ../data/features/short-ncrna-negative-control-dataset-features-1.csv 2

# ADD Histone H3K36me3
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K36me3/H3K36me3_lncrna-exon1-histone-feature.csv ../data/features/functional-lncrna-exon1-dataset-features-1-2.csv 3 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K36me3/H3K36me3_lncrna-exon2-histone-feature.csv ../data/features/functional-lncrna-exon2-dataset-features-1-2.csv 3 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K36me3/H3K36me3_protein-exon2-histone-feature.csv ../data/features/functional-protein-exon2-dataset-features-1-2.csv 3 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K36me3/H3K36me3_protein-exon3-histone-feature.csv ../data/features/functional-protein-exon3-dataset-features-1-2.csv 3 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K36me3/H3K36me3_short-ncrna-histone-feature.csv ../data/features/functional-short-ncrna-dataset-features-1-2.csv 3 && \

./add_feature_from_to.sh ../data/datasets/histone_feature/H3K36me3/H3K36me3_lncrna-exon1-NC-histone-feature.csv ../data/features/lncrna-exon1-negative-control-dataset-features-1-2.csv 3 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K36me3/H3K36me3_lncrna-exon2-NC-histone-feature.csv ../data/features/lncrna-exon2-negative-control-dataset-features-1-2.csv 3 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K36me3/H3K36me3_protein-exon2-NC-histone-feature.csv ../data/features/protein-exon2-negative-control-dataset-features-1-2.csv 3 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K36me3/H3K36me3_protein-exon3-NC-histone-feature.csv ../data/features/protein-exon3-negative-control-dataset-features-1-2.csv 3 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K36me3/H3K36me3_short-ncrna-NC-histone-feature.csv ../data/features/short-ncrna-negative-control-dataset-features-1-2.csv 3

# ADD Histone H3K4me3
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K4me3/H3K4me3_lncrna-exon1-histone-feature.csv ../data/features/functional-lncrna-exon1-dataset-features-1-2-3.csv 4 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K4me3/H3K4me3_lncrna-exon2-histone-feature.csv ../data/features/functional-lncrna-exon2-dataset-features-1-2-3.csv 4 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K4me3/H3K4me3_protein-exon2-histone-feature.csv ../data/features/functional-protein-exon2-dataset-features-1-2-3.csv 4 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K4me3/H3K4me3_protein-exon3-histone-feature.csv ../data/features/functional-protein-exon3-dataset-features-1-2-3.csv 4 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K4me3/H3K4me3_short-ncrna-histone-feature.csv ../data/features/functional-short-ncrna-dataset-features-1-2-3.csv 4 && \

./add_feature_from_to.sh ../data/datasets/histone_feature/H3K4me3/H3K4me3_lncrna-exon1-NC-histone-feature.csv ../data/features/lncrna-exon1-negative-control-dataset-features-1-2-3.csv 4 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K4me3/H3K4me3_lncrna-exon2-NC-histone-feature.csv ../data/features/lncrna-exon2-negative-control-dataset-features-1-2-3.csv 4 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K4me3/H3K4me3_protein-exon2-NC-histone-feature.csv ../data/features/protein-exon2-negative-control-dataset-features-1-2-3.csv 4 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K4me3/H3K4me3_protein-exon3-NC-histone-feature.csv ../data/features/protein-exon3-negative-control-dataset-features-1-2-3.csv 4 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K4me3/H3K4me3_short-ncrna-NC-histone-feature.csv ../data/features/short-ncrna-negative-control-dataset-features-1-2-3.csv 4

# ADD Histone H3K9me3
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K9me3/H3K9me3_lncrna-exon1-histone-feature.csv ../data/features/functional-lncrna-exon1-dataset-features-1-2-3-4.csv 5 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K9me3/H3K9me3_lncrna-exon2-histone-feature.csv ../data/features/functional-lncrna-exon2-dataset-features-1-2-3-4.csv 5 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K9me3/H3K9me3_protein-exon2-histone-feature.csv ../data/features/functional-protein-exon2-dataset-features-1-2-3-4.csv 5 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K9me3/H3K9me3_protein-exon3-histone-feature.csv ../data/features/functional-protein-exon3-dataset-features-1-2-3-4.csv 5 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K9me3/H3K9me3_short-ncrna-histone-feature.csv ../data/features/functional-short-ncrna-dataset-features-1-2-3-4.csv 5 && \

./add_feature_from_to.sh ../data/datasets/histone_feature/H3K9me3/H3K9me3_lncrna-exon1-NC-histone-feature.csv ../data/features/lncrna-exon1-negative-control-dataset-features-1-2-3-4.csv 5 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K9me3/H3K9me3_lncrna-exon2-NC-histone-feature.csv ../data/features/lncrna-exon2-negative-control-dataset-features-1-2-3-4.csv 5 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K9me3/H3K9me3_protein-exon2-NC-histone-feature.csv ../data/features/protein-exon2-negative-control-dataset-features-1-2-3-4.csv 5 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K9me3/H3K9me3_protein-exon3-NC-histone-feature.csv ../data/features/protein-exon3-negative-control-dataset-features-1-2-3-4.csv 5 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K9me3/H3K9me3_short-ncrna-NC-histone-feature.csv ../data/features/short-ncrna-negative-control-dataset-features-1-2-3-4.csv 5

# ADD Methylome
./add_feature_from_to.sh ../data/datasets/methylome_feature/lncrna-exon1-methylome-feature.csv ../data/features/functional-lncrna-exon1-dataset-features-1-2-3-4-5.csv 6 && \
./add_feature_from_to.sh ../data/datasets/methylome_feature/lncrna-exon2-methylome-feature.csv ../data/features/functional-lncrna-exon2-dataset-features-1-2-3-4-5.csv 6 && \
./add_feature_from_to.sh ../data/datasets/methylome_feature/protein-exon2-methylome-feature.csv ../data/features/functional-protein-exon2-dataset-features-1-2-3-4-5.csv 6 && \
./add_feature_from_to.sh ../data/datasets/methylome_feature/protein-exon3-methylome-feature.csv ../data/features/functional-protein-exon3-dataset-features-1-2-3-4-5.csv 6 && \
./add_feature_from_to.sh ../data/datasets/methylome_feature/short-ncrna-methylome-feature.csv ../data/features/functional-short-ncrna-dataset-features-1-2-3-4-5.csv 6 && \

./add_feature_from_to.sh ../data/datasets/methylome_feature/lncrna-exon1-NC-methylome-feature.csv ../data/features/lncrna-exon1-negative-control-dataset-features-1-2-3-4-5.csv 6 && \
./add_feature_from_to.sh ../data/datasets/methylome_feature/lncrna-exon2-NC-methylome-feature.csv ../data/features/lncrna-exon2-negative-control-dataset-features-1-2-3-4-5.csv 6 && \
./add_feature_from_to.sh ../data/datasets/methylome_feature/protein-exon2-NC-methylome-feature.csv ../data/features/protein-exon2-negative-control-dataset-features-1-2-3-4-5.csv 6 && \
./add_feature_from_to.sh ../data/datasets/methylome_feature/protein-exon3-NC-methylome-feature.csv ../data/features/protein-exon3-negative-control-dataset-features-1-2-3-4-5.csv 6 && \
./add_feature_from_to.sh ../data/datasets/methylome_feature/short-ncrna-NC-methylome-feature.csv ../data/features/short-ncrna-negative-control-dataset-features-1-2-3-4-5.csv 6

# ADD Histone H3K79me2
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K79me2/H3K79me2_lncrna-exon1-histone-feature.csv ../data/features/functional-lncrna-exon1-dataset-features-1-2-3-4-5-6.csv 7 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K79me2/H3K79me2_lncrna-exon2-histone-feature.csv ../data/features/functional-lncrna-exon2-dataset-features-1-2-3-4-5-6.csv 7 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K79me2/H3K79me2_protein-exon2-histone-feature.csv ../data/features/functional-protein-exon2-dataset-features-1-2-3-4-5-6.csv 7 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K79me2/H3K79me2_protein-exon3-histone-feature.csv ../data/features/functional-protein-exon3-dataset-features-1-2-3-4-5-6.csv 7 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K79me2/H3K79me2_short-ncrna-histone-feature.csv ../data/features/functional-short-ncrna-dataset-features-1-2-3-4-5-6.csv 7 && \

./add_feature_from_to.sh ../data/datasets/histone_feature/H3K79me2/H3K79me2_lncrna-exon1-NC-histone-feature.csv ../data/features/lncrna-exon1-negative-control-dataset-features-1-2-3-4-5-6.csv 7 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K79me2/H3K79me2_lncrna-exon2-NC-histone-feature.csv ../data/features/lncrna-exon2-negative-control-dataset-features-1-2-3-4-5-6.csv 7 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K79me2/H3K79me2_protein-exon2-NC-histone-feature.csv ../data/features/protein-exon2-negative-control-dataset-features-1-2-3-4-5-6.csv 7 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K79me2/H3K79me2_protein-exon3-NC-histone-feature.csv ../data/features/protein-exon3-negative-control-dataset-features-1-2-3-4-5-6.csv 7 && \
./add_feature_from_to.sh ../data/datasets/histone_feature/H3K79me2/H3K79me2_short-ncrna-NC-histone-feature.csv ../data/features/short-ncrna-negative-control-dataset-features-1-2-3-4-5-6.csv 7

cd ../data/features/ && \
cp functional-lncrna-exon1-dataset-features-1-2-3-4-5-6-7.csv functional-lncrna-exon1-dataset-features-augmented.csv && \
cp functional-lncrna-exon2-dataset-features-1-2-3-4-5-6-7.csv functional-lncrna-exon2-dataset-features-augmented.csv && \
cp functional-protein-exon2-dataset-features-1-2-3-4-5-6-7.csv functional-protein-exon2-dataset-features-augmented.csv && \
cp functional-protein-exon3-dataset-features-1-2-3-4-5-6-7.csv functional-protein-exon3-dataset-features-augmented.csv && \
cp functional-short-ncrna-dataset-features-1-2-3-4-5-6-7.csv functional-short-ncrna-dataset-features-augmented.csv && \
cp lncrna-exon1-negative-control-dataset-features-1-2-3-4-5-6-7.csv lncrna-exon1-negative-control-dataset-features-augmented.csv && \
cp lncrna-exon2-negative-control-dataset-features-1-2-3-4-5-6-7.csv lncrna-exon2-negative-control-dataset-features-augmented.csv && \
cp protein-exon2-negative-control-dataset-features-1-2-3-4-5-6-7.csv protein-exon2-negative-control-dataset-features-augmented.csv && \
cp protein-exon3-negative-control-dataset-features-1-2-3-4-5-6-7.csv protein-exon3-negative-control-dataset-features-augmented.csv && \
cp short-ncrna-negative-control-dataset-features-1-2-3-4-5-6-7.csv short-ncrna-negative-control-dataset-features-augmented.csv
