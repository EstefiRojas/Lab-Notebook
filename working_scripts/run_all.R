# This R script outlines the complete analysis pipeline. It loads data, performs calculations, runs various analyses, and generates plots.

## Step 1: Load Configuration and Libraries

# This chunk loads the main configuration file and all required data loading functions.

# Load configuration
source("scripts/config.R")

# Load data loading functions
source("scripts/load_gene_functionality_features.R")
source("scripts/load_dinucleotide_features.R")
source("scripts/load_epigenetic_features.R")
source("scripts/load_gene_functionality_zscores.R")
source("scripts/load_epigenetic_zscores.R")
source("scripts/utils.R")


## Step 2: Load Data and Compute Z-Scores

# This chunk loads the feature matrices and computes the robust z-scores.

# Load and compute z-scores
feature_matrix <- load_gene_functionality_features()
feature_matrix_epigenetic <- load_epigenetic_features()
source("scripts/compute_robust_zscores.R")
source("scripts/compute_ks_stats.R")

## Step 3: Run Analysis Scripts

# This chunk executes the main analysis scripts, including distance effect, MDS, and PCA analysis.

# Run analysis scripts
source("scripts/distance_effect_analysis.R")
source("scripts/mds_analisys.R")
source("scripts/pca_analysis.R")

## Step 4: Generate Plots

# This chunk generates all the plots for the analysis, including heatmaps and violin plots.

# Generate plots
source("scripts/generate_epigenetic_heatmap.R")
source("scripts/generate_intrinsic_heatmap.R")
source("scripts/generate_gene_functionality_heatmap.R")
source("scripts/generate_epigenetic_zscore_jitter_plots.R")
source("scripts/generate_epigenetic_zscore_violin_plots.R")
source("scripts/generate_gene_functionality_zscore_violin_plots.R")

## Conclusion

# The analysis pipeline has been completed successfully. The results can be found in the `results` directory.
print("Analysis pipeline completed successfully!")
