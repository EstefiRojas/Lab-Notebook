library(dplyr)
library(ggplot2)
# Install LambertW package
install.packages("LambertW")
library(LambertW)

# Compute z-scores for all gene functionality features
# 1. Load gene functionality features
source("load_gene_functionality_features.R")

# 2. Compute robust z-scores for each gene type
# Function to apply transformation
apply_transformation <- function(df, method = "log", lambda = NULL) {
  transformed_data <- df
  
  #if (length(data) == 0) return(NA)  # Return NA if no data left

  if (method == "log") {
    # Log transformation (ensure no zero values)
    transformed_data <- log(df + 1)
  } else if (method == "sqrt") {
    # Square root transformation
    transformed_data <- sqrt(df + 1)
  } else if (method == "boxcox" && !is.null(lambda)) {
    # Box-Cox transformation (requires positive values)
    transformed_data <- ifelse(lambda == 0, log(df), (df^lambda - 1) / lambda)
  } else if (method == "lambertw") {
    # LambertW transformation
    transformed_data[!is.na(df)] <- Gaussianize(df[!is.na(df)], type = "h")
    transformed_data[is.na(df)] <- NA
  } else if (method == "inverse") {
    # Inverse transformation
    transformed_data <- 1 / (df + 1)
  } else {
    stop("Invalid transformation method or parameters.")
  }
  return(transformed_data)
}

# Function to remove outliers
remove_outliers_vectorized <- function(df, method = "IQR", threshold = 3) {
  cleaned_data <- df
  
  if (method == "IQR") {
    for (col in colnames(df)) {
      column_data <- as.numeric(df[[col]])
      column_data <- na.omit(column_data)
      
      # Remove outliers using IQR method
      Q <- quantile(column_data, probs = c(0.25, 0.75))
      IQR_value <- IQR(column_data)
      lower_bound <- Q[1] - 1 * IQR_value
      upper_bound <- Q[2] + 1 * IQR_value
      
      # Assign NA to outliers
      cleaned_data[[col]][df[[col]] < lower_bound | df[[col]] > upper_bound] <- NA
    }
  } else if (method == "mad") {
    for (col in colnames(df)) {
      column_data <- as.numeric(df[[col]])
      column_data <- na.omit(column_data)
      
      #if (length(column_data) == 0) return(NA)  # Return NA if no data left
      
      median_val <- median(column_data, na.rm = TRUE)
      mad_val <- mad(column_data, na.rm = TRUE)
      
      if (mad_val != 0) {
        deviation <- abs(column_data - median_val)
        # Identify points that are more than `threshold` times MAD away from median
        #outliers <- column_data[deviation > threshold * mad_val]
        cleaned_data[[col]][deviation > threshold * mad_val] <- NA
      #  return(column_data)  # No outliers if MAD is zero
      }
    }
  } else {
    stop("Invalid transformation method or parameters.")
  }
  
  return(cleaned_data)
}

# Updated function to compute robust z-scores
compute_robust_zscores_optimized <- function(functional_df, negative_df, outlier_method = "mad") {
  # Remove outliers from both datasets (vectorized)
  negative_df_clean <- remove_outliers_vectorized(negative_df, method = outlier_method)
  functional_df_clean <- remove_outliers_vectorized(functional_df, method = outlier_method)
  
  # Calculate metrics for each column without looping
  mad_values <- apply(negative_df_clean, 2, function(col) mad(col, na.rm = TRUE))
  median_values <- apply(negative_df_clean, 2, function(col) median(col, na.rm = TRUE))
  mean_values <- apply(negative_df_clean, 2, function(col) mean(col, na.rm = TRUE))
  sd_values <- apply(negative_df_clean, 2, function(col) sd(col, na.rm = TRUE))
  meanAD_values <- apply(negative_df_clean, 2, function(col) mean(abs(col - mean(col, na.rm = TRUE)), na.rm = TRUE))
  
  # Calculate z-scores
  list_with_zscores <- list()
  list_with_method <- list()
  
  for (col in colnames(functional_df_clean)) {
    positive_col <- functional_df_clean[[col]]
    # Attempt to convert columns to numeric
    positive_col_numeric <- as.numeric(positive_col)
    
    # Identify non-convertible values
    non_numeric_pos_indices <- is.na(positive_col_numeric) & !is.na(positive_col)
    
    # Replace non-convertible values with NA
    positive_col_numeric[non_numeric_pos_indices] <- NA
    
    # Get statistics
    mad_value <- mad_values[col]
    median_value <- median_values[col]
    mean_value <- mean_values[col]
    sd_value <- sd_values[col]
    meanAD_value <- meanAD_values[col]
    
    # Initialize z-scores with NA for all elements
    zscores_col <- rep(NA, length(positive_col_numeric))
    
    # Compute z-scores
    if (!is.na(mad_value)) {
      if (!is.na(mad_value) && mad_value != 0) {
        zscores_col[!non_numeric_pos_indices] <- (positive_col_numeric[!non_numeric_pos_indices] - median_value) / mad_value
        list_with_method[[col]] <- "MAD"
      } else if (meanAD_value != 0) {
        zscores_col[!non_numeric_pos_indices] <- (positive_col_numeric[!non_numeric_pos_indices] - median_value) / (1.2533 * meanAD_value)
        list_with_method[[col]] <- "meanAD"
      } else if (sd_value != 0) {
        zscores_col[!non_numeric_pos_indices] <- (positive_col_numeric[!non_numeric_pos_indices] - mean_value) / sd_value
        list_with_method[[col]] <- "regular z-score"
      } else {
        list_with_method[[col]] <- "constant value (no variability)"
      }
    }  else {
      list_with_zscores[[col]] <- NA
      list_with_method[[col]] <- "NA"
    }
    
    list_with_zscores[[col]] <- zscores_col
  }
  
  return(list(zscores = list_with_zscores, method = list_with_method))
}

# Function to visualize data distributions
visualize_data_hist <- function(original_df, cleaned_df, col_name) {
  # Create a combined dataframe to facilitate visualization
  original_data <- data.frame(Value = original_df[[col_name]], Status = "Original")
  cleaned_data <- data.frame(Value = cleaned_df[[col_name]], Status = "Cleaned")
  combined_data <- rbind(original_data, cleaned_data)
  
  # Remove NA values for visualization
  combined_data <- na.omit(combined_data)
  
  # Create histogram for original and cleaned data
  ggplot(combined_data, aes(x = Value, fill = Status)) +
    geom_histogram(alpha = 0.6, position = "identity", bins = 100) +
    ggtitle(paste("Histogram for", col_name, "Before and After Outlier Removal")) +
    theme_minimal()
}

visualize_data_dist <- function(original_df, cleaned_df, col_name) {
  # Create a combined dataframe to facilitate visualization
  original_data <- data.frame(Value = original_df[[col_name]], Status = "Original")
  cleaned_data <- data.frame(Value = cleaned_df[[col_name]], Status = "Transformed")
  combined_data <- rbind(original_data, cleaned_data)
  
  # Remove NA values for visualization
  combined_data <- na.omit(combined_data)
  
  # Create histogram for original and cleaned data
  ggplot(combined_data, aes(x = Value, fill = Status)) +
    geom_density(alpha = 0.6) +
    ggtitle(paste("Density Distribution for", col_name, "Before and After LambertW transformation")) +
    theme_minimal()
}

visualize_data_box <- function(original_df, cleaned_df, col_name) {
  # Create a combined dataframe to facilitate visualization
  original_data <- data.frame(Value = original_df[[col_name]], Status = "Original")
  cleaned_data <- data.frame(Value = cleaned_df[[col_name]], Status = "Cleaned")
  combined_data <- rbind(original_data, cleaned_data)
  
  # Remove NA values for visualization
  combined_data <- na.omit(combined_data)
  
  # Create boxplot for original and cleaned data
  ggplot(combined_data, aes(x = Status, y = Value, fill = Status)) +
    geom_boxplot() +
    ggtitle(paste("Boxplot for", col_name, "Before and After Outlier Removal")) +
    theme_minimal()
}

################
# Filter only protein coding rows
protein_positive_feature_matrix <- feature_matrix_numeric %>% 
  filter(Dataset == "protein-coding-exon2" | Dataset == "protein-coding-exon3") %>%
  dplyr::select(-Functional, -Dataset)
summary(protein_positive_feature_matrix)
protein_negative_feature_matrix <- feature_matrix_numeric %>% 
  filter(Dataset == "protein-exon2-negative-control" | Dataset == "protein-exon3-negative-control") %>%
  dplyr::select(-Functional, -Dataset)

# Transform non-Gaussian features
# Positive
transformed_features <- protein_positive_feature_matrix %>%
  select(H3K27ac,H3K79me2) %>%
  apply_transformation(method = "lambertw")
protein_positive_feature_matrix[,colnames(transformed_features)] <- transformed_features

transformed_features <- protein_positive_feature_matrix %>%
  select(Expression,`RPKM_primary cell`,`Copy number`,`Repeat free`,H3K36me3,chromatin_acc) %>%
  apply_transformation(method = "log")
protein_positive_feature_matrix[,colnames(transformed_features)] <- transformed_features

# Negative controls
transformed_features_neg <- protein_negative_feature_matrix %>%
  select(H3K27ac,H3K79me2) %>%
  apply_transformation(method = "lambertw")
protein_negative_feature_matrix[,colnames(transformed_features)] <- transformed_features_neg

transformed_features_neg <- protein_negative_feature_matrix %>%
  select(Expression,`RPKM_primary cell`,`Copy number`,`Repeat free`,H3K36me3,chromatin_acc) %>%
  apply_transformation(method = "log")
protein_negative_feature_matrix[,colnames(transformed_features)] <- transformed_features_neg

# Compute zscores for protein coding:
# Compute z-scores
protein_functional_z_scores <- as.data.frame(
  compute_robust_zscores_optimized(protein_positive_feature_matrix, 
                         protein_negative_feature_matrix)$zscores)
summary(protein_functional_z_scores)
protein_negative_control_z_scores <- as.data.frame(
  compute_robust_zscores_optimized(protein_negative_feature_matrix,
                         protein_negative_feature_matrix)$zscores)
summary(protein_negative_control_z_scores)

# Replace Dataset column
protein_functional_z_scores$Dataset <- (feature_matrix_numeric %>% 
  filter(Dataset == "protein-coding-exon2" |
           Dataset == "protein-coding-exon3"))$Dataset

protein_negative_control_z_scores$Dataset <- (feature_matrix_numeric %>% 
  filter(Dataset == "protein-exon2-negative-control" |
           Dataset == "protein-exon3-negative-control"))$Dataset

# Split z-scores in their respective dataset
functional_protein_exon2_dataset_zscores <- protein_functional_z_scores %>%
  filter(Dataset == "protein-coding-exon2")
functional_protein_exon3_dataset_zscores <- protein_functional_z_scores %>%
  filter(Dataset == "protein-coding-exon3")

protein_exon2_negative_control_dataset_zscores <- protein_negative_control_z_scores %>%
  filter(Dataset == "protein-exon2-negative-control")
protein_exon3_negative_control_dataset_zscores <- protein_negative_control_z_scores %>%
  filter(Dataset == "protein-exon3-negative-control")

# 3. Save resulting z-scores in a csv file
write.csv(functional_protein_exon2_dataset_zscores, "../data/z-scores/functional-protein-exon2-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
write.csv(functional_protein_exon3_dataset_zscores, "../data/z-scores/functional-protein-exon3-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
write.csv(protein_exon2_negative_control_dataset_zscores, "../data/z-scores/protein-exon2-negative-control-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
write.csv(protein_exon3_negative_control_dataset_zscores, "../data/z-scores/protein-exon3-negative-control-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)

#############################
# Compute zscores for lncRNA:
# Filter only lncRNA rows
lncrna_positive_feature_matrix <- feature_matrix_numeric %>% 
  filter(Dataset == "lncrna-exon1" |
           Dataset == "lncrna-exon2" ) %>%
  dplyr::select(-Functional, -Dataset)
lncrna_negative_feature_matrix <- feature_matrix_numeric %>% 
  filter(Dataset == "lncrna-exon1-negative-control" |
           Dataset == "lncrna-exon2-negative-control") %>%
  dplyr::select(-Functional, -Dataset)

# Transform non-Gaussian features
# Positive
transformed_features <- lncrna_positive_feature_matrix %>%
  select(H3K27ac,H3K79me2) %>%
  apply_transformation(method = "lambertw")
lncrna_positive_feature_matrix[,colnames(transformed_features)] <- transformed_features

transformed_features <- lncrna_positive_feature_matrix %>%
  select(Expression,`RPKM_primary cell`,`Copy number`,`Repeat free`,H3K36me3,chromatin_acc) %>%
  apply_transformation(method = "log")
lncrna_positive_feature_matrix[,colnames(transformed_features)] <- transformed_features

# Negative controls
transformed_features_neg <- lncrna_negative_feature_matrix %>%
  select(H3K27ac,H3K79me2) %>%
  apply_transformation(method = "lambertw")
lncrna_negative_feature_matrix[,colnames(transformed_features)] <- transformed_features_neg

transformed_features_neg <- lncrna_negative_feature_matrix %>%
  select(Expression,`RPKM_primary cell`,`Copy number`,`Repeat free`,H3K36me3,chromatin_acc) %>%
  apply_transformation(method = "log")
lncrna_negative_feature_matrix[,colnames(transformed_features)] <- transformed_features_neg

# Compute z-scores
lncrna_functional_z_scores <- as.data.frame(
  compute_robust_zscores_optimized(lncrna_positive_feature_matrix,
                         lncrna_negative_feature_matrix)$zscores)
summary(lncrna_functional_z_scores)
lncrna_negative_control_z_scores <- as.data.frame(
  compute_robust_zscores_optimized(lncrna_negative_feature_matrix,
                         lncrna_negative_feature_matrix)$zscores)
summary(lncrna_negative_control_z_scores)

# Replace Dataset column
lncrna_functional_z_scores$Dataset <- (feature_matrix_numeric %>% 
                                         filter(Dataset == "lncrna-exon1" |
                                                  Dataset == "lncrna-exon2" ))$Dataset

lncrna_negative_control_z_scores$Dataset <- (feature_matrix_numeric %>% 
                                               filter(Dataset == "lncrna-exon1-negative-control" |
                                                        Dataset == "lncrna-exon2-negative-control"))$Dataset

# Split z-scores in their respective dataset
functional_lncrna_exon1_dataset_zscores <- lncrna_functional_z_scores %>%
  filter(Dataset == "lncrna-exon1")
functional_lncrna_exon2_dataset_zscores <- lncrna_functional_z_scores %>%
  filter(Dataset == "lncrna-exon2")

lncrna_exon1_negative_control_dataset_zscores <- lncrna_negative_control_z_scores %>%
  filter(Dataset == "lncrna-exon1-negative-control")
lncrna_exon2_negative_control_dataset_zscores <- lncrna_negative_control_z_scores %>%
  filter(Dataset == "lncrna-exon2-negative-control")

# 3. Save resulting z-scores in a csv file
write.csv(functional_lncrna_exon1_dataset_zscores, "../data/z-scores/functional-lncrna-exon1-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
write.csv(functional_lncrna_exon2_dataset_zscores, "../data/z-scores/functional-lncrna-exon2-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
write.csv(lncrna_exon1_negative_control_dataset_zscores, "../data/z-scores/lncrna-exon1-negative-control-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
write.csv(lncrna_exon2_negative_control_dataset_zscores, "../data/z-scores/lncrna-exon2-negative-control-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)


#############################
# Compute zscores for sncRNA:
# Filter only lncRNA rows
sncrna_positive_feature_matrix <- feature_matrix_numeric %>% 
  filter(Dataset == "short-ncrna" ) %>%
  dplyr::select(-Functional, -Dataset)
sncrna_negative_feature_matrix <- feature_matrix_numeric %>% 
  filter(Dataset == "short-ncrna-negative-control" ) %>%
  dplyr::select(-Functional, -Dataset)

# Transform non-Gaussian features
# Positive
transformed_features <- sncrna_positive_feature_matrix %>%
  select(H3K27ac,H3K79me2) %>%
  apply_transformation(method = "lambertw")
sncrna_positive_feature_matrix[,colnames(transformed_features)] <- transformed_features

transformed_features <- sncrna_positive_feature_matrix %>%
  select(Expression,`RPKM_primary cell`,`Copy number`,`Repeat free`,H3K36me3,chromatin_acc) %>%
  apply_transformation(method = "log")
sncrna_positive_feature_matrix[,colnames(transformed_features)] <- transformed_features

# Negative controls
transformed_features_neg <- sncrna_negative_feature_matrix %>%
  select(H3K27ac,H3K79me2) %>%
  apply_transformation(method = "lambertw")
sncrna_negative_feature_matrix[,colnames(transformed_features)] <- transformed_features_neg

transformed_features_neg <- sncrna_negative_feature_matrix %>%
  select(Expression,`RPKM_primary cell`,`Copy number`,`Repeat free`,H3K36me3,chromatin_acc) %>%
  apply_transformation(method = "log")
sncrna_negative_feature_matrix[,colnames(transformed_features)] <- transformed_features_neg

# Compute z-scores
sncrna_functional_z_scores <- as.data.frame(
  compute_robust_zscores_optimized(sncrna_positive_feature_matrix,
                         sncrna_negative_feature_matrix)$zscores)
summary(sncrna_functional_z_scores)
sncrna_negative_control_z_scores <- as.data.frame(
  compute_robust_zscores_optimized(sncrna_negative_feature_matrix,
                         sncrna_negative_feature_matrix)$zscores)
summary(sncrna_negative_control_z_scores)

# Replace Dataset column
sncrna_functional_z_scores$Dataset <- (feature_matrix_numeric %>% 
                                         filter(Dataset == "short-ncrna" ))$Dataset

sncrna_negative_control_z_scores$Dataset <- (feature_matrix_numeric %>% 
                                               filter(Dataset == "short-ncrna-negative-control" ))$Dataset

# Split z-scores in their respective dataset
functional_short_ncrna_dataset_zscores <- sncrna_functional_z_scores %>%
  filter(Dataset == "short-ncrna")
summary(functional_short_ncrna_dataset_zscores)
short_ncrna_negative_control_dataset_zscores <- sncrna_negative_control_z_scores %>%
  filter(Dataset == "short-ncrna-negative-control")
summary(short_ncrna_negative_control_dataset_zscores)
# 3. Save resulting z-scores in a csv file
write.csv(functional_short_ncrna_dataset_zscores, "../data/z-scores/functional-short-ncrna-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
write.csv(short_ncrna_negative_control_dataset_zscores, "../data/z-scores/short-ncrna-negative-control-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)

zscores_data <- rbind(protein_functional_z_scores, protein_negative_control_z_scores,
                      lncrna_functional_z_scores, lncrna_negative_control_z_scores,
                      sncrna_functional_z_scores, sncrna_negative_control_z_scores)
write.csv(zscores_data, "../data/results/zscores/gene-functionality-zscores-transformed.csv", row.names = FALSE)
