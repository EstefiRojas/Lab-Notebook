# Install necessary packages (if not already installed)
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

# Convert columns to numeric
convert_to_numeric <- function(x) {
  as.numeric(x)
}

# Remove outliers using IQR method
remove_outliers_IQR <- function(df, col) {
  df_col <- df %>% select({{col}}) %>% na.omit()
  Q_col <- quantile(df_col[[1]], probs = c(0.25, 0.75))
  iqr_col <- IQR(df_col[[1]])
  lower_bound <- Q_col[1] - 1.5 * iqr_col
  upper_bound <- Q_col[2] + 1.5 * iqr_col
  df_filtered <- df_col %>% filter(df_col[[1]] > lower_bound & df_col[[1]] < upper_bound)
  return(df_filtered[[1]])
}

# Get robust z-scores for selected features
get_robust_zscores <- function(functional_df, negative_df, selected_features) {
  list_with_zscores <- list()
  list_with_method <- list()
  for (col in selected_features) {
    #positive_col <- remove_outliers_IQR(functional_df, col)
    #negative_col <- remove_outliers_IQR(negative_df, col)
    positive_col <- functional_df[[col]]
    negative_col <- negative_df[[col]]
    if(length(positive_col) == 0) {
      positive_col <- functional_df[[col]]
    }
    if(length(negative_col) == 0) {
      negative_col <- negative_df[[col]]
    }
    mad_value <- mad(negative_col, na.rm = TRUE)
    median_value <- median(negative_col, na.rm = TRUE)
    Sn_value <- Sn(negative_col, na.rm = TRUE)
    Qn_value <- Qn(negative_col, na.rm = TRUE)
    tau_value <- scaleTau2(negative_col, na.rm = TRUE)
    mean_value <- mean(negative_col, na.rm = TRUE)
    meanAD_value <- mean(abs(negative_col - mean_value), na.rm = TRUE)
    
    if(mad_value != 0) { # MAD method
      list_with_zscores[[col]] <- abs((positive_col - median_value) / mad_value)
      list_with_method[[col]] <- "MAD"
    } else if(Sn_value != 0) { # Sn estimator
      list_with_zscores[[col]] <- abs((positive_col - median_value) / Sn_value)
      list_with_method[[col]] <- "Sn estimator"
    } else if(Qn_value != 0) { # # Qn estimator
      list_with_zscores[[col]] <- abs((positive_col - median_value) / Qn_value)
      list_with_method[[col]] <- "Qn estimator"
    } else if(tau_value != 0) { # tau estimator
      list_with_zscores[[col]] <- abs((positive_col - median_value) / tau_value)
      list_with_method[[col]] <- "tau estimator"
    } else if(meanAD_value != 0) { # meanAD estimator
      list_with_zscores[[col]] <- abs((positive_col - median_value) / (1.2533 * meanAD_value))
      list_with_method[[col]] <- "meanAD"
    } else { # regular z-score
      list_with_zscores[[col]] <- abs((positive_col - mean(negative_col)) / sd(negative_col, na.rm = TRUE))
      list_with_method[[col]] <- "regular z-score"
    }
    #list_with_zscores[[col]] <- 1.1926 * (positive_col - median_value) / Sn_value
  }
  return(list(zscores=list_with_zscores, method=list_with_method))
}

# Function to compute regular z-scores for each column in `selected_features` for functional and non functional datasets
get_zscores <-function(functional_df, negative_df, selected_features) {
  list_with_zscores <- list()
  for (col in selected_features) {
    #positive_col <- remove_outliers_IQR(functional_df, col)
    #negative_col <- remove_outliers_IQR(negative_df, col)
    positive_col <- functional_df[[col]]
    negative_col <- negative_df[[col]]
    
    mean_col <- mean(negative_col, na.rm = TRUE)
    sd_col <- sd(negative_col, na.rm = TRUE)
    
    list_with_zscores[[col]] <- abs((positive_col - mean_col) / sd_col)
  }
  return(list_with_zscores)
}

# Normalize robust z-scores data
normalize_data <- function(prot_df, lncrna_df, sncrna_df, prot_neg_df, lncrna_neg_df, sncrna_neg_df, selected_features) {
  feature_matrix <- bind_rows(prot_df, lncrna_df, sncrna_df)
  feature_matrix_neg <- bind_rows(prot_neg_df, lncrna_neg_df, sncrna_neg_df)
  
  zscores <- get_robust_zscores(feature_matrix, feature_matrix_neg, selected_features)$zscores
  #zscores_neg <- get_robust_zscores(feature_matrix_neg, feature_matrix_neg, selected_features)$zscores
  
  data_normalized_all <- as.data.frame(zscores) %>%
    mutate(across(everything(), convert_to_numeric))
  
  return(data_normalized_all)
}

# Normalize regular z-scores data
normalize_regular_data <- function(prot_df, lncrna_df, sncrna_df, prot_neg_df, lncrna_neg_df, sncrna_neg_df, selected_features) {
  feature_matrix <- bind_rows(prot_df, lncrna_df, sncrna_df)
  feature_matrix_neg <- bind_rows(prot_neg_df, lncrna_neg_df, sncrna_neg_df)
  
  zscores <- get_zscores(feature_matrix, feature_matrix_neg, selected_features)
  #zscores_neg <- get_robust_zscores(feature_matrix_neg, feature_matrix_neg, selected_features)$zscores
  
  data_normalized_all <- as.data.frame(zscores) %>%
    mutate(across(everything(), convert_to_numeric))
  
  return(data_normalized_all)
}

# Perform PCA by Feature
perform_pca <- function(data_normalized, scale=FALSE) {
  data_t_normalized <- t(data_normalized)
  pca_result <- prcomp(data_t_normalized, scale. = scale)
  return(pca_result)
}

# Plot PCA by Feature
plot_pca <- function(pca_result, selected_features) {
  pc_scores <- as.data.frame(pca_result$x[, 1:2])
  pca_df <- pc_scores %>%
    mutate(feature_name = selected_features) %>%
    mutate(feature_name = factor(feature_name, levels = selected_features))
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = feature_name, shape = feature_name)) +
    geom_point(size = 5, position = position_dodge(width = 0.1)) +
    scale_shape_manual(values = c(10, 15, 18, 18, 18, 18, 18, 7, 16, 16, 5, 5, 8, 8, 14, 6, 6, 6, 6, 11, 11, 13, 13, 13, 13, 13)) +
    scale_color_manual(values = c("red", "green", "lightpink1", "royalblue", "lightcyan4", "lawngreen", "indianred4", "darkturquoise", "dodgerblue4", "forestgreen", "purple", "orange3", "firebrick", "gold3", "sienna2", "limegreen", "khaki4", "orangered", "maroon", "purple3", "darkgreen", "lightcoral", "darkorchid4", "tomato4", "chartreuse4", "cornflowerblue")) +
    labs(title = "Regular Z-scores Principal Component Analysis by Feature", x = "PC1 (72.93%)", y = "PC2 (27.07%)") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 26),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 24),
      legend.text = element_text(size = 20),
      legend.key.size = unit(2, "lines")
    )
}

################
# Main Execution
################
selected_features <- c("Random number", "GC content", "CpG","GA","GG","TA","TC",
                     "low_complexity_density","phyloP max_241w","phyloP max_100w",
                     "RPKM_tissue","RPKM_primary cell","Copy number","Repeat free","RNAcode","Max covariance",
                     "MFE","RNAalifold","Interaction_ave","gnomAD_SNP_density","gnomAD_MAF",
                     "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome")

# Robust Z-scores analysis
data_normalized <- normalize_data(
  feature_matrix_numeric_prot, feature_matrix_numeric_lncrna, feature_matrix_numeric_sncrna,
  feature_matrix_numeric_protein_neg, feature_matrix_numeric_lncrna_neg, feature_matrix_numeric_sncrna_neg,
  selected_features
)

summary(data_normalized)
data_normalized <- na.omit(data_normalized)
pca_result <- perform_pca(data_normalized, TRUE)
summary(pca_result)
plot_pca(pca_result, selected_features)


# Regular Z-scores analysis
regular_data_normalized <- normalize_regular_data(
  feature_matrix_numeric_prot, feature_matrix_numeric_lncrna, feature_matrix_numeric_sncrna,
  feature_matrix_numeric_protein_neg, feature_matrix_numeric_lncrna_neg, feature_matrix_numeric_sncrna_neg,
  selected_features
)

summary(regular_data_normalized)
regular_data_normalized <- na.omit(regular_data_normalized)
pca_result <- perform_pca(regular_data_normalized, FALSE)
summary(pca_result)
plot_pca(pca_result, selected_features)
