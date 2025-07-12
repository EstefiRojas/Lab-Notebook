# LOESS Curve Analysis: Delta Z-scores based on distance to gene

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(robustbase)  # For robust statistics functions

# Define constants
SELECT_FEATURES <- c("GC_percentage","CG")
DISTANCE_BINS <- list(
  "1-5k" = c(0, 5000),
  "5k-50k" = c(5000, 50000),
  "50k-500k" = c(50000, 500000),
  "500k-5M" = c(500000, 5000000)
)

# Helper functions
load_dataset <- function(name, file_path) {
  data.frame(Dataset = name, read.csv(file_path, header = TRUE))
}

subset_by_distance_to_gene <- function(dataset) {
  lapply(DISTANCE_BINS, function(range) {
    dataset[dataset$DistanceGene > range[1] & dataset$DistanceGene <= range[2], ]
  })
}

get_gene_coords <- function(dataset) {
  dataset %>%
    mutate(
      GeneUpstreamEnd = pmax(1, Start - DistanceGene - 1),
      GeneDownstreamStart = End + DistanceGene + 1
    )
}

get_upstream_and_downstream_sets <- function(positive_df, negative_df, selected_features) {
  bind_rows(
    inner_join(negative_df, positive_df, by = join_by(Chromosome, GeneDownstreamStart == Start)),
    inner_join(negative_df, positive_df, by = join_by(Chromosome, GeneUpstreamEnd == End))
  ) %>%
    select(all_of(c(paste0(selected_features, c(".x", ".y")), "DistanceGene")))
}

get_robust_zscores <- function(subsets_df, selected_features, is_negative = FALSE) {
  suffix <- if(is_negative) ".y" else ".x"
  
  lapply(selected_features, function(col) {
    col_name <- paste0(col, suffix)
    positive_col <- subsets_df[[col_name]]
    negative_col <- subsets_df[[paste0(col, ".y")]]
    
    estimators <- list(
      MAD = mad(negative_col, na.rm = TRUE),
      Sn = Sn(negative_col, na.rm = TRUE),
      Qn = Qn(negative_col, na.rm = TRUE),
      tau = scaleTau2(negative_col, na.rm = TRUE),
      meanAD = mean(abs(negative_col - mean(negative_col, na.rm = TRUE)), na.rm = TRUE)
    )
    
    non_zero_estimators <- estimators[estimators != 0]
    
    if (length(non_zero_estimators) > 0) {
      chosen_estimator <- non_zero_estimators[[1]]
      method <- names(non_zero_estimators)[1]
      (positive_col - median(negative_col, na.rm = TRUE)) / chosen_estimator
    } else {
      (positive_col - mean(negative_col, na.rm = TRUE)) / sd(negative_col, na.rm = TRUE)
    }
  })
}

compute_zscores_all_bins <- function(subsets_df, selected_features) {
  lapply(subsets_df, function(bin) {
    list(
      positive = get_robust_zscores(bin, selected_features),
      negative = get_robust_zscores(bin, selected_features, is_negative = TRUE)
    )
  })
}

calculate_delta_zscores <- function(zscores_df1, zscores_df2) {
  Map(`-`, zscores_df2, zscores_df1)
}

# Main script
main <- function() {
  # Load datasets
  datasets <- list(
    funcProtExon2 = load_dataset("protein-coding-exon2", "../data/latest1000all/functional-protein-exon2-dataset-features-w-fic.csv"),
    funcProtExon3 = load_dataset("protein-coding-exon3", "../data/latest1000all/functional-protein-exon3-dataset-features-w-fic.csv"),
    funcLncrnaExon1 = load_dataset("lncrna-exon1", "../data/latest1000all/functional-lncrna-exon1-dataset-features-w-fic.csv"),
    funcLncrnaExon2 = load_dataset("lncrna-exon2", "../data/latest1000all/functional-lncrna-exon2-dataset-features-w-fic.csv"),
    funcSncRNA = load_dataset("short-ncrna", "../data/latest1000all/functional-short-ncrna-dataset-features-w-fic.csv"),
    protExon2NC = load_dataset("protein-exon2-negative-control", "../data/latest1000all/protein-exon2-negative-control-dataset-features-w-fic.csv"),
    protExon3NC = load_dataset("protein-exon3-negative-control", "../data/latest1000all/protein-exon3-negative-control-dataset-features-w-fic.csv"),
    lncrnaExon1NC = load_dataset("lncrna-exon1-negative-control", "../data/latest1000all/lncrna-exon1-negative-control-dataset-features-w-fic.csv"),
    lncrnaExon2NC = load_dataset("lncrna-exon2-negative-control", "../data/latest1000all/lncrna-exon2-negative-control-dataset-features-w-fic.csv"),
    sncrnaNC = load_dataset("short-ncrna-negative-control", "../data/latest1000all/short-ncrna-negative-control-dataset-features-w-fic.csv")
  )
  
  # Process negative datasets
  negative_datasets <- lapply(datasets[6:10], get_gene_coords)
  
  # Get upstream and downstream sets
  supersets <- Map(get_upstream_and_downstream_sets, 
                   datasets[1:5], negative_datasets, 
                   MoreArgs = list(selected_features = SELECT_FEATURES))
  
  # Subset data into corresponding bins
  subsets <- lapply(supersets, subset_by_distance_to_gene)
  
  # Compute z-scores
  zscores_bins <- lapply(subsets, compute_zscores_all_bins, SELECT_FEATURES)
  
  # Compute deltas
  delta_zscores <- lapply(zscores_bins, function(bins) {
    lapply(bins, function(bin) {
      calculate_delta_zscores(as.data.frame(bin$positive), as.data.frame(bin$negative))
    })
  })
  
  # Prepare data for plotting
  plot_data <- bind_rows(lapply(names(delta_zscores$funcLncrnaExon1), function(bin) {
    data.frame(
      bin = bin,
      GC_percentage = delta_zscores$funcLncrnaExon1[[bin]]$GC_percentage,
      distanceGene = subsets$funcLncrnaExon1[[bin]]$DistanceGene
    )
  }))
  
  # Create plot
  ggplot(plot_data, aes(x = bin, y = GC_percentage, fill = bin)) +
    geom_violin(scale = "width", show.legend = TRUE) +
    geom_boxplot(alpha = 0.0, outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.2) +
    facet_wrap(~ "GC_content", scales = "free") + 
    labs(title = "Delta Z-Scores for lincRNA Exon1 GC Content", x = "Distance to Gene", y = "Delta Z-score") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_brewer(palette = "Blues") +
    coord_cartesian(ylim = c(-3, 8))
}

# Run the main function
main()
