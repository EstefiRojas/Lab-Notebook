# Configuration file for project paths

# Base directory for data
DATA_DIR <- "data"

# Base directory for results
RESULTS_DIR <- "results"

# Subdirectories within data
FEATURES_DIR <- file.path(DATA_DIR, "features")

# Subdirectories within results
DINUCLEOTIDE_FEATURES_DIR <- file.path(RESULTS_DIR, "dinucleotide_features")
EPIGENETIC_FEATURES_DIR <- file.path(RESULTS_DIR, "epigenetic_features")
EPIGENETIC_HISTONES_DIR <- file.path(EPIGENETIC_FEATURES_DIR, "histones")
EPIGENETIC_CHR_ACC_DIR <- file.path(EPIGENETIC_FEATURES_DIR, "chrm_acc")
EPIGENETIC_METHYLOME_DIR <- file.path(EPIGENETIC_FEATURES_DIR, "methylome")

ZSCORES_DIR <- file.path(RESULTS_DIR, "zscores")
EPIGENETIC_ZSCORES_DIR <- file.path(RESULTS_DIR, "epigenetic_zscores")

KS_STATS_DIR <- file.path(RESULTS_DIR, "ks_stats")

DISTANCE_EFFECT_DIR <- file.path(RESULTS_DIR, "distance_effect")
DISTANCE_EFFECT_HEATSCATTER_DIR <- file.path(DISTANCE_EFFECT_DIR, "heatscatter")
DISTANCE_EFFECT_SPEARMAN_DIR <- file.path(DISTANCE_EFFECT_DIR, "spearman")

VIOLIN_PLOTS_DIR <- file.path(RESULTS_DIR, "violinPlots")
VIOLIN_PLOTS_SUBSET_DIR <- file.path(VIOLIN_PLOTS_DIR, "Subset")

PCA_DIR <- file.path(RESULTS_DIR, "pca")

# Specific file paths
# Data files
FUNC_PROT_EXON2_FEATURES_FILE <- file.path(FEATURES_DIR, "functional-protein-exon2-dataset-features.csv")
FUNC_PROT_EXON3_FEATURES_FILE <- file.path(FEATURES_DIR, "functional-protein-exon3-dataset-features.csv")
FUNC_LNCRNA_EXON1_FEATURES_FILE <- file.path(FEATURES_DIR, "functional-lncrna-exon1-dataset-features.csv")
FUNC_LNCRNA_EXON2_FEATURES_FILE <- file.path(FEATURES_DIR, "functional-lncrna-exon2-dataset-features.csv")
FUNC_SNCRNA_FEATURES_FILE <- file.path(FEATURES_DIR, "functional-short-ncrna-dataset-features.csv")

NC_PROT_EXON2_FEATURES_FILE <- file.path(FEATURES_DIR, "protein-exon2-negative-control-dataset-features.csv")
NC_PROT_EXON3_FEATURES_FILE <- file.path(FEATURES_DIR, "protein-exon3-negative-control-dataset-features.csv")
NC_LNCRNA_EXON1_FEATURES_FILE <- file.path(FEATURES_DIR, "lncrna-exon1-negative-control-dataset-features.csv")
NC_LNCRNA_EXON2_FEATURES_FILE <- file.path(FEATURES_DIR, "lncrna-exon2-negative-control-dataset-features.csv")
NC_SNCRNA_FEATURES_FILE <- file.path(FEATURES_DIR, "short-ncrna-negative-control-dataset-features.csv")

# Dinucleotide feature files
PROT_EXON2_DIN_FEATURES_FILE <- file.path(DINUCLEOTIDE_FEATURES_DIR, "protein-exon2-dinucleotide-features.csv")
PROT_EXON3_DIN_FEATURES_FILE <- file.path(DINUCLEOTIDE_FEATURES_DIR, "protein-exon3-dinucleotide-features.csv")
LNCRNA_EXON1_DIN_FEATURES_FILE <- file.path(DINUCLEOTIDE_FEATURES_DIR, "lncrna-exon1-dinucleotide-features.csv")
LNCRNA_EXON2_DIN_FEATURES_FILE <- file.path(DINUCLEOTIDE_FEATURES_DIR, "lncrna-exon2-dinucleotide-features.csv")
SNCRNA_DIN_FEATURES_FILE <- file.path(DINUCLEOTIDE_FEATURES_DIR, "short-ncrna-dinucleotide-features.csv")

NC_PROT_EXON2_DIN_FEATURES_FILE <- file.path(DINUCLEOTIDE_FEATURES_DIR, "protein-exon2-NC-dinucleotide-features.csv")
NC_PROT_EXON3_DIN_FEATURES_FILE <- file.path(DINUCLEOTIDE_FEATURES_DIR, "protein-exon3-NC-dinucleotide-features.csv")
NC_LNCRNA_EXON1_DIN_FEATURES_FILE <- file.path(DINUCLEOTIDE_FEATURES_DIR, "lncrna-exon1-NC-dinucleotide-features.csv")
NC_LNCRNA_EXON2_DIN_FEATURES_FILE <- file.path(DINUCLEOTIDE_FEATURES_DIR, "lncrna-exon2-NC-dinucleotide-features.csv")
NC_SNCRNA_DIN_FEATURES_FILE <- file.path(DINUCLEOTIDE_FEATURES_DIR, "short-ncrna-NC-dinucleotide-features.csv")

# Epigenetic histone files
COMBINED_PROT_EXON2_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_protein-exon2.csv")
COMBINED_PROT_EXON3_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_protein-exon3.csv")
COMBINED_LNCRNA_EXON1_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_lncrna-exon1.csv")
COMBINED_LNCRNA_EXON2_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_lncrna-exon2.csv")
COMBINED_SNCRNA_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_short-ncrna.csv")

COMBINED_NC_PROT_EXON2_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_protein-exon2-NC.csv")
COMBINED_NC_PROT_EXON3_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_protein-exon3-NC.csv")
COMBINED_NC_LNCRNA_EXON1_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_lncrna-exon1-NC.csv")
COMBINED_NC_LNCRNA_EXON2_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_lncrna-exon2-NC.csv")
COMBINED_NC_SNCRNA_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_short-ncrna-NC.csv")

# Epigenetic chromatin accessibility files
CHR_ACC_PROT_EXON2_FILE <- file.path(EPIGENETIC_CHR_ACC_DIR, "chrm_acc_protein-exon2-chrm_acc-feature.csv")
CHR_ACC_PROT_EXON3_FILE <- file.path(EPIGENETIC_CHR_ACC_DIR, "chrm_acc_protein-exon3-chrm_acc-feature.csv")
CHR_ACC_LNCRNA_EXON1_FILE <- file.path(EPIGENETIC_CHR_ACC_DIR, "chrm_acc_lncrna-exon1-chrm_acc-feature.csv")
CHR_ACC_LNCRNA_EXON2_FILE <- file.path(EPIGENETIC_CHR_ACC_DIR, "chrm_acc_lncrna-exon2-chrm_acc-feature.csv")
CHR_ACC_SNCRNA_FILE <- file.path(EPIGENETIC_CHR_ACC_DIR, "chrm_acc_short-ncrna-chrm_acc-feature.csv")

NC_CHR_ACC_PROT_EXON2_FILE <- file.path(EPIGENETIC_CHR_ACC_DIR, "chrm_acc_protein-exon2-NC-chrm_acc-feature.csv")
NC_CHR_ACC_PROT_EXON3_FILE <- file.path(EPIGENETIC_CHR_ACC_DIR, "chrm_acc_protein-exon3-NC-chrm_acc-feature.csv")
NC_CHR_ACC_LNCRNA_EXON1_FILE <- file.path(EPIGENETIC_CHR_ACC_DIR, "chrm_acc_lncrna-exon1-NC-chrm_acc-feature.csv")
NC_CHR_ACC_LNCRNA_EXON2_FILE <- file.path(EPIGENETIC_CHR_ACC_DIR, "chrm_acc_lncrna-exon2-NC-chrm_acc-feature.csv")
NC_CHR_ACC_SNCRNA_FILE <- file.path(EPIGENETIC_CHR_ACC_DIR, "chrm_acc_short-ncrna-NC-chrm_acc-feature.csv")

# Epigenetic methylome files
METHYLOME_PROT_EXON2_FILE <- file.path(EPIGENETIC_METHYLOME_DIR, "protein-exon2-methylome-feature.csv")
METHYLOME_PROT_EXON3_FILE <- file.path(EPIGENETIC_METHYLOME_DIR, "protein-exon3-methylome-feature.csv")
METHYLOME_LNCRNA_EXON1_FILE <- file.path(EPIGENETIC_METHYLOME_DIR, "lncrna-exon1-methylome-feature.csv")
METHYLOME_LNCRNA_EXON2_FILE <- file.path(EPIGENETIC_METHYLOME_DIR, "lncrna-exon2-methylome-feature.csv")
METHYLOME_SNCRNA_FILE <- file.path(EPIGENETIC_METHYLOME_DIR, "short-ncrna-methylome-feature.csv")

NC_METHYLOME_PROT_EXON2_FILE <- file.path(EPIGENETIC_METHYLOME_DIR, "protein-exon2-NC-methylome-feature.csv")
NC_METHYLOME_PROT_EXON3_FILE <- file.path(EPIGENETIC_METHYLOME_DIR, "protein-exon3-NC-methylome-feature.csv")
NC_METHYLOME_LNCRNA_EXON1_FILE <- file.path(EPIGENETIC_METHYLOME_DIR, "lncrna-exon1-NC-methylome-feature.csv")
NC_METHYLOME_LNCRNA_EXON2_FILE <- file.path(EPIGENETIC_METHYLOME_DIR, "lncrna-exon2-NC-methylome-feature.csv")
NC_METHYLOME_SNCRNA_FILE <- file.path(EPIGENETIC_METHYLOME_DIR, "short-ncrna-NC-methylome-feature.csv")

# Z-score files
MRNA_ZSCORES_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "mrna_z_scores.csv")
LNCRNA_ZSCORES_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "lncrna_z_scores.csv")
SNCRNA_ZSCORES_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "sncrna_z_scores.csv")

MRNA_NC_ZSCORES_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "mrna_negative_control_z_scores.csv")
LNCRNA_NC_ZSCORES_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "lncrna_negative_control_z_scores.csv")
SNCRNA_NC_ZSCORES_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "sncrna_negative_control_z_scores.csv")

MRNA_ZSCORES_SUBSET_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "mrna_z_scores_subset.csv")
LNCRNA_ZSCORES_SUBSET_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "lncrna_z_scores_subset.csv")
SNCRNA_ZSCORES_SUBSET_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "sncrna_z_scores_subset.csv")

MRNA_NC_ZSCORES_SUBSET_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "mrna_negative_control_z_scores_subset.csv")
LNCRNA_NC_ZSCORES_SUBSET_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "lncrna_negative_control_z_scores_subset.csv")
SNCRNA_NC_ZSCORES_SUBSET_FILE <- file.path(EPIGENETIC_ZSCORES_DIR, "sncrna_negative_control_z_scores_subset.csv")

# KS Stats files
EPIGENETIC_MRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "epigenetic_mrna_ks_stat_sample.csv")
EPIGENETIC_LNCRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "epigenetic_lncrna_ks_stat_sample.csv")
EPIGENETIC_SNCRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "epigenetic_sncrna_ks_stat_sample.csv")
EPIGENETIC_MEAN_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "epigenetic_mean_ks_stat_sample.csv")
INTRINSIC_MRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "intrinsic_mrna_ks_stat_sample.csv")
INTRINSIC_LNCRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "intrinsic_lncrna_ks_stat_sample.csv")
INTRINSIC_SNCRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "intrnisic_sncrna_ks_stat_sample.csv")
INTRINSIC_MEAN_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "intrinsic_mean_ks_stat_sample.csv")
GENE_FUNCTIONALITY_MRNA_KS_STAT_FILE <- file.path(KS_STATS_DIR, "gene_functionality_mrna_ks_stat_sample.csv")
GENE_FUNCTIONALITY_LNCRNA_KS_STAT_FILE <- file.path(KS_STATS_DIR, "gene_functionality_lncrna_ks_stat_sample.csv")
GENE_FUNCTIONALITY_SNCRNA_KS_STAT_FILE <- file.path(KS_STATS_DIR, "gene_functionality_sncrna_ks_stat_sample.csv")
GENE_FUNCTIONALITY_MEAN_KS_STAT_FILE <- file.path(KS_STATS_DIR, "gene_functionality_mean_ks_stat_sample.csv")

# Distance effect files
EFFECT_ON_DISTANCE_JOINED_FILE <- file.path(DISTANCE_EFFECT_HEATSCATTER_DIR, "effectOnDistanceJoined5k_5M.png")
CORR_MATRIX_PROTEIN_FILE <- file.path(DISTANCE_EFFECT_SPEARMAN_DIR, "corr_matrix_protein.csv")
CORR_MATRIX_SNCRNA_FILE <- file.path(DISTANCE_EFFECT_SPEARMAN_DIR, "corr_matrix_sncrna.csv")
CORR_MATRIX_LNCRNA_FILE <- file.path(DISTANCE_EFFECT_SPEARMAN_DIR, "corr_matrix_lncrna.csv")

# Plot files
EPIGENETIC_ZSCORES_JITTER_PLOT_FILE <- file.path(RESULTS_DIR, "epigenetic_zscores_jitter.png")
INTRINSIC_JITTER_PLOTS_DIR <- file.path(RESULTS_DIR, "jitterPlots", "intrinsic")

# Violin plot files
EPIGENETIC_P1_1_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p1_1.png")
EPIGENETIC_P1_2_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p1_2.png")
EPIGENETIC_P2_1_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p2_1.png")
EPIGENETIC_P2_2_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p2_2.png")
EPIGENETIC_P3_1_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p3_1.png")
EPIGENETIC_P4_1_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p4_1.png")
EPIGENETIC_P5_1_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p5_1.png")

CONSERVATION_REVIEW_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "conservation_review.png")
EXPRESSION_REVIEW_TEST_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "expression_review_test.png")
INTRINSIC2_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "intrinsic2.png")
EPIGENETIC_REVIEW_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "epigenetic_review.png")
EPIGENETIC_REVIEW2_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "epigenetic_review2.png")
EPIGENETIC_REVIEW3_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "epigenetic_review3.png")
EPIGENETIC_MAXSCALEDSIGNAL_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "epigenetic_maxScaledSignal.png")
EPIGENETIC_MAXSCALEDSIGNAL2_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "epigenetic_maxScaledSignal2.png")
EPIGENETIC_MAXSCALEDSIGNAL3_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "epigenetic_maxSacaledSignal3.png")
ASSOCIATION_REVIEW_TRIMTRUE_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "association_review_trimtrue.png")
ASSOCIATION_REVIEW2_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "association_review2.png")
STRUCTURE_REVIEW_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "structure_review.png")
STRUCTURE_REVIEW2_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "structure_review2.png")
STRUCTURE_REVIEW3_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "structure_review3.png")
POP_VAR_REVIEW_PLOT_FILE <- file.path(VIOLIN_PLOTS_DIR, "popvar_review.png")

# PCA plot files
PCA_PROTEIN_20_FEATURES_PLOT_FILE <- file.path(PCA_DIR, "pca_protein_20_features.png")
PCA_LNCRNA_20_FEATURES_PLOT_FILE <- file.path(PCA_DIR, "pca_lncrna_20_features.png")
PCA_SNCRNA_20_FEATURES_PLOT_FILE <- file.path(PCA_DIR, "pca_sncrna_20_features.png")
PCA_20_FEATURES_JOINED_PLOT_FILE <- file.path(PCA_DIR, "pca_20_features_joined.png")
PCA_PROTEIN_20_FEATURES_LOADINGS_PLOT_FILE <- file.path(PCA_DIR, "pca_protein_20_features_loadings.png")
PCA_LNCRNA_20_FEATURES_LOADINGS_PLOT_FILE <- file.path(PCA_DIR, "pca_lncrna_20_features_loadings.png")
PCA_SNCRNA_20_FEATURES_LOADINGS_PLOT_FILE <- file.path(PCA_DIR, "pca_sncrna_20_features_loadings.png")
PCA_20_FEATURES_LOADINGS_JOINED_PLOT_FILE <- file.path(PCA_DIR, "pca_20_features_loadings_joined.png")

# Latest log scatter plots directory
LATEST_LOG_SCATTER_PLOTS_DIR <- file.path(RESULTS_DIR, "latest1000all", "logScatterPlots", "Paper")
