# Master script to run the entire analysis pipeline

# Load configuration
source("scripts/config.R")

# Load data loading functions
source("scripts/load_gene_functionality_features.R")
source("scripts/load_dinucleotide_features.R")
source("scripts/load_epigenetic_features.R")
source("scripts/load_gene_functionality_zscores.R")
source("scripts/load_epigenetic_zscores.R")
source("scripts/load_ks_stats_utils.R")

# Load and compute z-scores
feature_matrix <- load_gene_functionality_features()
feature_matrix_epigenetic <- load_epigenetic_features()
source("scripts/compute_robust_zscores.R")

# Run analysis scripts
source("scripts/analize_distance_effect.R")
source("scripts/mds_analisys.R")
source("scripts/pca_analysis.R")

# Generate plots
source("scripts/generate_epigenetic_heatmap.R")
source("scripts/generate_intrinsic_heatmap.R")
source("scripts/generate_epigenetic_zscore_jitter_plots.R")
source("scripts/generate_epigenetic_zscore_violin_plots.R")
source("scripts/generate_gene_functionality_zscore_violin_plots.R")

print("Analysis pipeline completed successfully!")
