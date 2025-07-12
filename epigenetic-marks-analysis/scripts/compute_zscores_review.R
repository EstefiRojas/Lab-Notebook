###############
# Script to Compute z-scores for all gene functionality features separated by Gene Type.
# The aim is to explore the resulting z-scores and determine the reason for some 
# extremely big values obtained for some features.
# Can be executed from the command line using:
#
#  Rscript compute_zscores_review.R /path/to/features/gene_functionality_features.csv
#
###############
# Load necessary libraries
check.packages <- function (pkg) {
  print("Installing required packages, please wait...")
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, library, character.only = TRUE)
}
check.packages(c("dplyr","ggplot2"))

main <- function() {
  ###############
  # Load csv file containing gene functionality features, passed as argument in the command line
  file <- commandArgs(trailingOnly = TRUE)[1]
  feature_matrix <- read.csv(file,
                             header = TRUE,
                             check.names = FALSE)
  #source("load_gene_functionality_features.R")
  
  # Convert all features to numeric values
  feature_matrix_numeric_pipeline <- feature_matrix_numeric %>%
    dplyr::select(-Dataset, -Functional) %>%
    sapply(function(feature) as.numeric(as.character(feature))) %>%
    as.data.frame()
  
  # Print summary stats of raw features
  cat("All features joined",sep = "\n")
  print(summary(feature_matrix_numeric_pipeline))
  cat("",sep = "\n")
  
  
  ##############################
  # MAIN EXECUTION:
  # Using the methods defined at the end of this script, we will compute, visualize, and save z-scores 
  # for protein coding, lncRNA, and sncRNA datasets separately.

  # Protein Coding Datasets:
  # FILTER only protein coding rows (TODO: better not to duplicate datasets, still thinking how to avoid this)
  protein_positive_feature_matrix <- feature_matrix_numeric_pipeline %>% 
    filter(Dataset == "protein-coding-exon2" | Dataset == "protein-coding-exon3") %>%
    dplyr::select(-Functional, -Dataset)
  
  protein_negative_feature_matrix <- feature_matrix_numeric_pipeline %>% 
    filter(Dataset == "protein-exon2-negative-control" | Dataset == "protein-exon3-negative-control") %>%
    dplyr::select(-Functional, -Dataset)
  
  #PRINT summary stats without outlier removal
  cat("mRNA(+) features",sep = "\n")
  print(summary(protein_positive_feature_matrix))
  cat("",sep = "\n")
  cat("mRNA(-) features",sep = "\n")
  print(summary(protein_negative_feature_matrix))
  cat("",sep = "\n")
  
  # ALSO print MAD, MeanAD, absolute deviations
  mad_values <- apply(protein_negative_feature_matrix, 2, function(col) median(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE))
  median_values <- apply(protein_negative_feature_matrix, 2, function(col) median(col, na.rm = TRUE))
  meanAD_values <- apply(protein_negative_feature_matrix, 2, function(col) mean(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE))
  
  cat("MAD mRNA(-)",sep = "\n")
  print(mad_values)
  cat("",sep = "\n")
  cat("meanAD mRNA(-)",sep = "\n")
  print(meanAD_values)
  cat("",sep = "\n")
  cat("median mRNA(-)",sep = "\n")
  print(median_values)
  cat("",sep = "\n")
  
  # COMPUTE zscores for protein coding:
  # Compute z-scores for functional set
  protein_functional_z_scores <- as.data.frame(
    compute_robust_zscores_optimized(protein_positive_feature_matrix, 
                           protein_negative_feature_matrix)$zscores)
  # Compute z-scores for negative control set
  protein_negative_control_z_scores <- as.data.frame(
    compute_robust_zscores_optimized(protein_negative_feature_matrix,
                           protein_negative_feature_matrix)$zscores)
  
  # PRINT SUMMARY STATS for Z_SCORES
  cat("ZSCORES mRNA(+)",sep = "\n")
  print(summary(protein_functional_z_scores))
  cat("",sep = "\n")
  cat("ZSCORES mRNA(-)",sep = "\n")
  print(summary(protein_negative_control_z_scores))
  cat("",sep = "\n")
  
  # Generate plots to visualize z-scores distributions for some features:
  # GC%
  visualize_data_hist(protein_functional_z_scores, "GC.content", "mRNA(+)")
  visualize_data_dist(protein_functional_z_scores, "GC.content", "mRNA(+)")
  visualize_data_box(protein_functional_z_scores, "GC.content", "mRNA(+)")
  
  visualize_data_hist(protein_negative_control_z_scores, "GC.content", "mRNA(-)")
  visualize_data_dist(protein_negative_control_z_scores, "GC.content", "mRNA(-)")
  visualize_data_box(protein_negative_control_z_scores, "GC.content", "mRNA(-)")
  
  # EXPRESSION (RPKM max TISSUE)
  visualize_data_hist(protein_functional_z_scores, "Expression", "mRNA(+)")
  visualize_data_dist(protein_functional_z_scores, "Expression", "mRNA(+)")
  visualize_data_box(protein_functional_z_scores, "Expression", "mRNA(+)")
  
  visualize_data_hist(protein_negative_control_z_scores, "Expression", "mRNA(-)")
  visualize_data_dist(protein_negative_control_z_scores, "Expression", "mRNA(-)")
  visualize_data_box(protein_negative_control_z_scores, "Expression", "mRNA(-)")
  
  
  ###############
  # lncRNA datasets:
  # Filter only lncRNA rows
  lncrna_positive_feature_matrix <- feature_matrix_numeric_pipeline %>% 
    filter(Dataset == "lncrna-exon1" |
             Dataset == "lncrna-exon2" ) %>%
    dplyr::select(-Functional, -Dataset)
  lncrna_negative_feature_matrix <- feature_matrix_numeric_pipeline %>% 
    filter(Dataset == "lncrna-exon1-negative-control" |
             Dataset == "lncrna-exon2-negative-control") %>%
    dplyr::select(-Functional, -Dataset)
  
  #PRINT summary stats without outlier removal
  cat("lncRNA(+) features",sep = "\n")
  print(summary(lncrna_positive_feature_matrix))
  cat("",sep = "\n")
  cat("lncRNA(-) features",sep = "\n")
  print(summary(lncrna_negative_feature_matrix))
  cat("",sep = "\n")
  
  # ALSO print MAD, MeanAD, absolute deviations
  mad_values <- apply(lncrna_negative_feature_matrix, 2, function(col) median(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE))
  median_values <- apply(lncrna_negative_feature_matrix, 2, function(col) median(col, na.rm = TRUE))
  meanAD_values <- apply(lncrna_negative_feature_matrix, 2, function(col) mean(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE))
  
  cat("MAD lncRNA(-)",sep = "\n")
  print(mad_values)
  cat("",sep = "\n")
  cat("meanAD lncRNA(-)",sep = "\n")
  print(meanAD_values)
  cat("",sep = "\n")
  cat("median lncRNA(-)",sep = "\n")
  print(median_values)
  cat("",sep = "\n")
  
  # Compute z-scores:
  lncrna_functional_z_scores <- as.data.frame(
    compute_robust_zscores_optimized(lncrna_positive_feature_matrix,
                                     lncrna_negative_feature_matrix)$zscores)
  
  lncrna_negative_control_z_scores <- as.data.frame(
    compute_robust_zscores_optimized(lncrna_negative_feature_matrix,
                                     lncrna_negative_feature_matrix)$zscores)
  
  # PRINT SUMMARY STATS for Z_SCORES
  cat("ZSCORES lncRNA(+)",sep = "\n")
  print(summary(lncrna_functional_z_scores))
  cat("",sep = "\n")
  cat("ZSCORES lncRNA(-)",sep = "\n")
  print(summary(lncrna_negative_control_z_scores))
  cat("",sep = "\n")
  
  # Generate plots to visualize z-scores distributions for some features:
  # GC%
  visualize_data_hist(lncrna_functional_z_scores, "GC.content", "lncRNA(+)")
  visualize_data_dist(lncrna_functional_z_scores, "GC.content", "lncRNA(+)")
  visualize_data_box(lncrna_functional_z_scores, "GC.content", "lncRNA(+)")
  
  visualize_data_hist(lncrna_negative_control_z_scores, "GC.content", "lncRNA(-)")
  visualize_data_dist(lncrna_negative_control_z_scores, "GC.content", "lncRNA(-)")
  visualize_data_box(lncrna_negative_control_z_scores, "GC.content", "lncRNA(-)")
  
  # EXXPRESSION
  visualize_data_hist(lncrna_functional_z_scores, "Expression", "lncRNA(+)")
  visualize_data_dist(lncrna_functional_z_scores, "Expression", "lncRNA(+)")
  visualize_data_box(lncrna_functional_z_scores, "Expression", "lncRNA(+)")
  
  visualize_data_hist(lncrna_negative_control_z_scores, "Expression", "lncRNA(-)")
  visualize_data_dist(lncrna_negative_control_z_scores, "Expression", "lncRNA(-)")
  visualize_data_box(lncrna_negative_control_z_scores, "Expression", "lncRNA(-)")
  
  ###############
  
  # sncRNA datasets:
  # Filter only sncRNA rows
  sncrna_positive_feature_matrix <- feature_matrix_numeric_pipeline %>% 
    filter(Dataset == "short-ncrna" ) %>%
    dplyr::select(-Functional, -Dataset)
  sncrna_negative_feature_matrix <- feature_matrix_numeric_pipeline %>% 
    filter(Dataset == "short-ncrna-negative-control" ) %>%
    dplyr::select(-Functional, -Dataset)
  
  #PRINT summary stats without outlier removal
  cat("sncRNA(+) features",sep = "\n")
  print(summary(sncrna_positive_feature_matrix))
  cat("",sep = "\n")
  cat("sncRNA(-) features",sep = "\n")
  print(summary(sncrna_negative_feature_matrix))
  cat("",sep = "\n")
  
  # ALSO print MAD, MeanAD, absolute deviations
  mad_values <- apply(sncrna_negative_feature_matrix, 2, function(col) median(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE))
  median_values <- apply(sncrna_negative_feature_matrix, 2, function(col) median(col, na.rm = TRUE))
  meanAD_values <- apply(sncrna_negative_feature_matrix, 2, function(col) mean(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE))
  
  cat("MAD sncRNA(-)",sep = "\n")
  print(mad_values)
  cat("",sep = "\n")
  cat("meanAD sncRNA(-)",sep = "\n")
  print(meanAD_values)
  cat("",sep = "\n")
  cat("median sncRNA()",sep = "\n")
  print(median_values)
  cat("",sep = "\n")
  
  # Compute z-scores
  sncrna_functional_z_scores <- as.data.frame(
    compute_robust_zscores_optimized(sncrna_positive_feature_matrix,
                                     sncrna_negative_feature_matrix)$zscores)
  sncrna_negative_control_z_scores <- as.data.frame(
    compute_robust_zscores_optimized(sncrna_negative_feature_matrix,
                                     sncrna_negative_feature_matrix)$zscores)
  
  # PRINT SUMMARY STATS for Z_SCORES
  cat("ZSCORES sncRNA(+)",sep = "\n")
  print(summary(sncrna_functional_z_scores))
  cat("ZSCORES sncRNA(-)",sep = "\n")
  print(summary(sncrna_negative_control_z_scores))
  # Generate plots to visualize z-scores distributions for some features:
  # GC%
  visualize_data_hist(sncrna_functional_z_scores, "GC.content", "sncRNA(+)")
  visualize_data_dist(sncrna_functional_z_scores, "GC.content", "sncRNA(+)")
  visualize_data_box(sncrna_functional_z_scores, "GC.content", "sncRNA(+)")
  
  visualize_data_hist(sncrna_negative_control_z_scores, "GC.content", "sncRNA(-)")
  visualize_data_dist(sncrna_negative_control_z_scores, "GC.content", "sncRNA(-)")
  visualize_data_box(sncrna_negative_control_z_scores, "GC.content", "sncRNA(-)")
  
  # EXPRESSION
  visualize_data_hist(sncrna_functional_z_scores, "Expression", "sncRNA(+)")
  visualize_data_dist(sncrna_functional_z_scores, "Expression", "sncRNA(+)")
  visualize_data_box(sncrna_functional_z_scores, "Expression", "sncRNA(+)")
  
  visualize_data_hist(sncrna_negative_control_z_scores, "Expression", "sncRNA(-)")
  visualize_data_dist(sncrna_negative_control_z_scores, "Expression", "sncRNA(-)")
  visualize_data_box(sncrna_negative_control_z_scores, "Expression", "sncRNA(-)")
  
  
  ########################
  # SAVE RESULTING Z-SCORES TO CSV FILE
  # Add Dataset column to corresponding z-score dataframes
  protein_functional_z_scores$Dataset <- (feature_matrix_numeric_pipeline %>% 
    filter(Dataset == "protein-coding-exon2" |
             Dataset == "protein-coding-exon3"))$Dataset
  
  protein_negative_control_z_scores$Dataset <- (feature_matrix_numeric_pipeline %>% 
    filter(Dataset == "protein-exon2-negative-control" |
             Dataset == "protein-exon3-negative-control"))$Dataset
  
  # Split z-scores in their respective exon dataset
  functional_protein_exon2_dataset_zscores <- protein_functional_z_scores %>%
    filter(Dataset == "protein-coding-exon2")
  functional_protein_exon3_dataset_zscores <- protein_functional_z_scores %>%
    filter(Dataset == "protein-coding-exon3")
  
  protein_exon2_negative_control_dataset_zscores <- protein_negative_control_z_scores %>%
    filter(Dataset == "protein-exon2-negative-control")
  protein_exon3_negative_control_dataset_zscores <- protein_negative_control_z_scores %>%
    filter(Dataset == "protein-exon3-negative-control")
  
  # SAVE resulting z-scores in a csv file
  write.csv(functional_protein_exon2_dataset_zscores, "../data/z-scores/functional-protein-exon2-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
  write.csv(functional_protein_exon3_dataset_zscores, "../data/z-scores/functional-protein-exon3-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
  write.csv(protein_exon2_negative_control_dataset_zscores, "../data/z-scores/protein-exon2-negative-control-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
  write.csv(protein_exon3_negative_control_dataset_zscores, "../data/z-scores/protein-exon3-negative-control-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
  
  
  #############################
  # Replace Dataset column
  lncrna_functional_z_scores$Dataset <- (feature_matrix_numeric_pipeline %>% 
                                           filter(Dataset == "lncrna-exon1" |
                                                    Dataset == "lncrna-exon2" ))$Dataset
  
  lncrna_negative_control_z_scores$Dataset <- (feature_matrix_numeric_pipeline %>% 
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
  
  # SAVE resulting z-scores in a csv file
  write.csv(functional_lncrna_exon1_dataset_zscores, "../data/z-scores/functional-lncrna-exon1-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
  write.csv(functional_lncrna_exon2_dataset_zscores, "../data/z-scores/functional-lncrna-exon2-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
  write.csv(lncrna_exon1_negative_control_dataset_zscores, "../data/z-scores/lncrna-exon1-negative-control-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
  write.csv(lncrna_exon2_negative_control_dataset_zscores, "../data/z-scores/lncrna-exon2-negative-control-dataset-zscores-transformed.csv", quote=FALSE, row.names=FALSE)
  
  
  # Replace Dataset column
  sncrna_functional_z_scores$Dataset <- (feature_matrix_numeric_pipeline %>% 
                                           filter(Dataset == "short-ncrna" ))$Dataset
  
  sncrna_negative_control_z_scores$Dataset <- (feature_matrix_numeric_pipeline %>% 
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
  
  # Join all z-scores in one dataframe and save that as well.
  zscores_data <- rbind(protein_functional_z_scores, protein_negative_control_z_scores,
                        lncrna_functional_z_scores, lncrna_negative_control_z_scores,
                        sncrna_functional_z_scores, sncrna_negative_control_z_scores)
  write.csv(zscores_data, "../data/results/zscores/gene-functionality-zscores-transformed.csv", row.names = FALSE)
}

##############
# FUNCTION DEFINITIONS:
# Function to COMPUTE robust z-scores
compute_robust_zscores_optimized <- function(functional_df, negative_df) {
  # Remove outliers from both datasets (ACORDING TO MEETING, DONT DO THIS PART, MIGHT BE ADDING MORE DIFFERENCES)
  #negative_df_clean <- remove_outliers_vectorized(negative_df, method = "IQR")
  #functional_df_clean <- remove_outliers_vectorized(functional_df, method = "IQR")
  
  # Calculate Stats ONLY for NEGATIVE set for each feature (MAKE THIS MORE CLEAR)
  mad_values    <- apply(negative_df, 2, function(col) median(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE))
  median_values <- apply(negative_df, 2, function(col) median(col, na.rm = TRUE))
  meanAD_values <- apply(negative_df, 2, function(col) mean(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE))
  
  # Initialize empty lists for z-scores and stat used
  list_with_zscores <- list()
  list_with_method <- list()
  
  # Process each feature one by one
  for (col in colnames(functional_df)) {
    positive_col <- functional_df[[col]]
    # Attempt to convert columns to numeric
    positive_col_numeric <- as.numeric(positive_col)
    
    # Identify non-convertible values (RENAMED THIS VARIABLE TO numeric_indices AND INVERT BOOLEANS)
    is_numeric_indices <- !( is.na(positive_col_numeric) & !is.na(positive_col) )
    
    # Overwrite non-convertible values with NA
    positive_col_numeric[ !is_numeric_indices ] <- NA
    
    # Retrieve statistics for negative controls
    mad_value <- mad_values[col]
    median_value <- median_values[col]
    meanAD_value <- meanAD_values[col]
    
    # Initialize z-scores with NA
    zscores_col <- rep(NA, length(positive_col_numeric))
    
    # Compute robust z-scores
    if (mad_value > 1e-6) { #think about very small mad values, treat them as 0, e.g. 10^-6, print warning to check that feature
      zscores_col[is_numeric_indices] <- (positive_col_numeric[is_numeric_indices] - median_value) / (1.4826 * mad_value)
      list_with_method[[col]] <- "MAD"
    } else if (meanAD_value > 1e-6) {
      zscores_col[is_numeric_indices] <- (positive_col_numeric[is_numeric_indices] - median_value) / (1.2533 * meanAD_value)
      list_with_method[[col]] <- "meanAD"
    } else {
      warning(paste(col,"has tiny variance. MAD =", mad_value, "meanAD = ", meanAD_value))
      list_with_zscores[[col]] <- NA
      list_with_method[[col]] <- "NA"
    }
    
    list_with_zscores[[col]] <- zscores_col
  }
  
  return(list(zscores = list_with_zscores, method = list_with_method))
}

# Function to remove outliers (NOT CALLED ANY MORE)
remove_outliers_vectorized <- function(df, method = "IQR") {
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
  } else {
    stop("Invalid transformation method or parameters.")
  }
  
  return(cleaned_data)
}

# Function to visualize data distributions
visualize_data_hist <- function(original_df, col_name, title) {
  # Create a combined dataframe to facilitate visualization
  original_data <- data.frame(Value = original_df[[col_name]], Status = "Original")

  # Remove NA values for visualization
  original_data <- na.omit(original_data)
  
  # Create histogram for original and cleaned data
  return(ggplot(original_data, aes(x = Value, fill = Status)) +
    geom_histogram(alpha = 0.6, position = "identity", bins = 100) +
    labs(x = "", y = "Robust z-score") +
    ggtitle(paste(title, col_name)) +
    theme_minimal() +
    theme(legend.position = "none"))
}

visualize_data_dist <- function(original_df, col_name, title) {
  # Create a combined dataframe to facilitate visualization
  original_data <- data.frame(Value = original_df[[col_name]], Status = "Original")

  # Remove NA values for visualization
  original_data <- na.omit(original_data)
  
  # Create histogram for original and cleaned data
  return(ggplot(original_data, aes(x = Value, fill = Status)) +
    geom_density(alpha = 0.6) +
    labs(x = "", y = "Robust z-score") +
    ggtitle(paste(title, col_name)) +
    theme_minimal() +
    theme(legend.position = "none"))
}

visualize_data_box <- function(original_df, col_name, title) {
  # Create a combined dataframe to facilitate visualization
  original_data <- data.frame(Value = original_df[[col_name]], Status = "Original")
 
  # Remove NA values for visualization
  original_data <- na.omit(original_data)
  
  # Create boxplot for original and cleaned data
  return(ggplot(original_data, aes(x = Status, y = Value, fill = Status)) +
    geom_boxplot() +
    labs(x = "", y = "Robust z-score") +
    ggtitle(paste(title, col_name)) +
    theme_minimal() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(),
          strip.text = element_blank()))
}

main()
