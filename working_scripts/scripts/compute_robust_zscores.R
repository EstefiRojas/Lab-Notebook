###############
# Utility script to Compute robust z-scores for all numeric features separated by Gene Type.
###############
#options(repos = c(CRAN = "https://cran.r-project.org"))
library(dplyr)
library(ggplot2)
library(tidyr)
#install.packages("patchwork", dependencies = TRUE)
library(patchwork)
library(scales)

source("scripts/load_gene_functionality_features.R")
feature_matrix <- load_gene_functionality_features()

# Load epigenetic features
source("scripts/load_epigenetic_features.R")
feature_matrix_epigenetic <- load_epigenetic_features()

# Function to compute robust z-scores
compute_robust_zscores_optimized <- function(positive_matrix, negative_matrix, verbose = FALSE) {
  # Use NEGATIVE CONTROLS median absolute deviations (MAD) and medians to compute robust z-scores
  mad_values <- apply(negative_matrix, 2, mad, na.rm = TRUE)  # Calculate MAD for each column in the negative control set
  median_values <- apply(negative_matrix, 2, median, na.rm = TRUE)  # Calculate median for each column in the negative control set
  meanAD_values <- apply(negative_matrix, 2, function(col) mean(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE))
  
  # Handle case where MAD is 0 or very small by using meanAD instead
  mad_values <- apply(negative_matrix, 2, function(col) {
    mad_value <- mad(col, na.rm = TRUE)
    if (mad_value <= 1e-6) {  # If MAD is 0 or very small, calculate mean absolute deviation (meanAD)
      meanAD <- mean(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE)
      warning(sprintf("Column %s with very small MAD of %f encountered: Using meanAD value of %f.\n", col, mad_value, meanAD))
      return(1.2533 * meanAD)  # Use meanAD if it is greater than the threshold
    } else {
      return(mad_value)  # Use MAD if it is sufficient
    }
  })
  
  # Calculate z-scores: subtract the median and divide by MAD for each feature
  zscores <- sweep(positive_matrix, 2, median_values, "-")  # Center each column by subtracting the median
  zscores <- sweep(zscores, 2, mad_values, "/")  # Scale each column by dividing by MAD
  
  #print("\n\nMedian values:\n")
  #print(median_values)
  #print("MAD values:\n")
  #print(mad_values)
  #print("meanAD values:\n")
  #print(meanAD_values)
  
  zscores
}

# Run z-score computation function. Receives a dataframe and an output directory
get_robust_zscores <- function(data, output_dir) {
  # Capture command line arguments
  #args <- commandArgs(trailingOnly = TRUE)
  #if (length(args) < 2) {
  #  stop("Please provide the path to the data file as an argument and output folder.")  # Ensure a file path is provided
  #}
  #names(feature_matrix)[1] <- "row"
  # Load and preprocess data
  #data_file <- args[1]  # Get the file path from the command line arguments
  #data <- read.csv(data_file, header = TRUE, check.names = TRUE)# Load the data from the specified CSV file
  #data <- feature_matrix %>% dplyr::select(-Dataset,-ID,-Functional,-Chromosome,-Start,-End,-Sequence)
  #data <- data_numeric_epigenetic_sample
  
  # Check if subset option is provided
  #subset_value <- as.numeric(gsub("--subset=", "", args[grep("--subset=", args)]))
  #if (is.na(subset_value)) {
  #subset_value <- nrow(data)  # Use nrow to indicate no subsetting by default
  #subset_value <- 100
  #}
  
  # Get output directory
  #output_dir <- args[2]
  #output_dir <- "../results/epigenetic_zscores"
  
  # Convert to numeric
  data_numeric <- data %>% 
    dplyr::select(-Dataset) %>%
    sapply(function(feature) as.numeric(as.character(feature))) %>%
    as.data.frame()
  #summary(data_numeric)
  # Restore Dataset column
  data_numeric$Dataset <- data$Dataset
  #unique(data_numeric$Dataset)
  #table(data_numeric$Dataset)
  #Protein coding
  # Separate the data into positive and negative feature matrices
  protein_positive_feature_matrix <- data_numeric %>% 
    filter(Dataset == "protein-coding-exon2" | Dataset == "protein-coding-exon3") #%>% 
  #dplyr::select(-Dataset)
  protein_negative_feature_matrix <- data_numeric %>% 
    filter(Dataset == "protein-exon2-negative-control" | Dataset == "protein-exon3-negative-control")# %>% 
  #dplyr::select(-Dataset)
  
  # Select random samples from positive and negative datasets
  #set.seed(2025)  # Set seed for reproducibility
  #protein_positive_feature_matrix <- protein_positive_feature_matrix %>% sample_n(min(2 * subset_value, nrow(protein_positive_feature_matrix)))  # Subset rows from positive data
  #protein_negative_feature_matrix <- protein_negative_feature_matrix %>% sample_n(min(2 * 10 * subset_value, nrow(protein_negative_feature_matrix)))  # Subset rows from negative data
  
  #print(summary(protein_positive_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #print(summary(protein_negative_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))  
  #print(apply(protein_positive_feature_matrix, 2, sd, na.rm = TRUE))
  #apply(protein_negative_feature_matrix, 2, sd, na.rm = TRUE)
  #print(sd(protein_negative_feature_matrix$RPKM_tissue, na.rm = TRUE))
  # Compute z-scores
  protein_functional_z_scores <- compute_robust_zscores_optimized(protein_positive_feature_matrix %>% dplyr::select(-Dataset), 
                                                                  protein_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                                  verbose)
  protein_negative_z_scores <- compute_robust_zscores_optimized(protein_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                                protein_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                                verbose)
  
  #print(summary(protein_functional_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #print(count(protein_positive_feature_matrix))
  #print(summary(protein_negative_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #print(summary(protein_negative_feature_matrix))
  #print(sd(protein_functional_z_scores$RPKM_tissue, na.rm = TRUE))
  #print(sd(protein_negative_z_scores$RPKM_tissue, na.rm = TRUE))
  
  #Restore Dataset
  protein_functional_z_scores$Dataset <- (protein_positive_feature_matrix %>% 
                                            filter(Dataset == "protein-coding-exon2" | Dataset == "protein-coding-exon3") %>% 
                                            dplyr::select(Dataset))$Dataset
  protein_negative_z_scores$Dataset <- (protein_negative_feature_matrix %>% 
                                          filter(Dataset == "protein-exon2-negative-control" | Dataset == "protein-exon3-negative-control") %>% 
                                          dplyr::select(Dataset))$Dataset
  #Join in one dataframe
  protein_zscores_all <- rbind(protein_functional_z_scores,protein_negative_z_scores)
  # Save z-scores to CSV
  z_scores_file <- file.path(output_dir, "mrna_z_scores.csv")
  write.csv(protein_functional_z_scores, z_scores_file, row.names = FALSE)
  
  z_scores_file <- file.path(output_dir, "mrna_negative_control_z_scores.csv")
  write.csv(protein_negative_z_scores, z_scores_file, row.names = FALSE)
  
  # Save raw data subset
  subset_file <- file.path(output_dir, "mrna_features.csv")
  write.csv(protein_positive_feature_matrix, subset_file, row.names = FALSE)
  
  subset_file <- file.path(output_dir, "mrna_negative_control_features.csv")
  write.csv(protein_negative_feature_matrix, subset_file, row.names = FALSE)
  
  #lncRNA
  # Separate the data into positive and negative feature matrices
  lncrna_positive_feature_matrix <- data_numeric %>% 
    filter(Dataset == "lncrna-exon1" | Dataset == "lncrna-exon2") #%>% 
  #dplyr::select(-Dataset)
  lncrna_negative_feature_matrix <- data_numeric %>% 
    filter(Dataset == "lncrna-exon1-negative-control" | Dataset == "lncrna-exon2-negative-control") #%>% 
  #dplyr::select(-Dataset)
  
  # Select random samples from positive and negative datasets
  #set.seed(2025)  # Set seed for reproducibility
  #lncrna_positive_feature_matrix <- lncrna_positive_feature_matrix %>% sample_n(min(2 * subset_value, nrow(lncrna_positive_feature_matrix)))  # Subset rows from positive data
  #lncrna_negative_feature_matrix <- lncrna_negative_feature_matrix %>% sample_n(min(2 * 10 * subset_value, nrow(lncrna_negative_feature_matrix)))  # Subset rows from negative data
  
  # Compute z-scores
  lncrna_functional_z_scores <- compute_robust_zscores_optimized(lncrna_positive_feature_matrix %>% dplyr::select(-Dataset), 
                                                                 lncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                                 verbose)
  lncrna_negative_z_scores <- compute_robust_zscores_optimized(lncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                               lncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                               verbose)
  
  #print(summary(lncrna_positive_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #count(lncrna_positive_feature_matrix)
  #print(summary(lncrna_negative_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #count(lncrna_negative_feature_matrix)
  
  #print(summary(lncrna_functional_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #print(summary(lncrna_negative_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  
  #Restore Dataset
  lncrna_functional_z_scores$Dataset <- (lncrna_positive_feature_matrix %>% 
                                           filter(Dataset == "lncrna-exon1" | Dataset == "lncrna-exon2") %>% 
                                           dplyr::select(Dataset))$Dataset
  lncrna_negative_z_scores$Dataset <- (lncrna_negative_feature_matrix %>% 
                                         filter(Dataset == "lncrna-exon1-negative-control" | Dataset == "lncrna-exon2-negative-control") %>% 
                                         dplyr::select(Dataset))$Dataset
  
  # Save z-scores to CSV
  z_scores_file <- file.path(output_dir, "lncrna_z_scores.csv")
  write.csv(lncrna_functional_z_scores, z_scores_file, row.names = FALSE)
  
  z_scores_file <- file.path(output_dir, "lncrna_negative_control_z_scores.csv")
  write.csv(lncrna_negative_z_scores, z_scores_file, row.names = FALSE)
  
  # Save raw subset feature to CSV
  subset_file <- file.path(output_dir, "lncrna_features.csv")
  write.csv(lncrna_positive_feature_matrix, subset_file, row.names = FALSE)
  
  subset_file <- file.path(output_dir, "lncrna_negative_control_features.csv")
  write.csv(lncrna_negative_feature_matrix, subset_file, row.names = FALSE)
  
  #sncRNA
  # Separate the data into positive and negative feature matrices
  sncrna_positive_feature_matrix <- data_numeric %>% 
    filter(Dataset == "short-ncrna" | Dataset == "short-ncrna") #%>% 
  #dplyr::select(-Dataset)
  sncrna_negative_feature_matrix <- data_numeric %>% 
    filter(Dataset == "short-ncrna-negative-control" | Dataset == "short-ncrna-negative-control") #%>% 
  #dplyr::select(-Dataset)
  
  # Select random samples from positive and negative datasets
  #set.seed(2025)  # Set seed for reproducibility
  #sncrna_positive_feature_matrix <- sncrna_positive_feature_matrix %>% sample_n(min(subset_value, nrow(sncrna_positive_feature_matrix)))  # Subset rows from positive data
  #sncrna_negative_feature_matrix <- sncrna_negative_feature_matrix %>% sample_n(min(10 * subset_value, nrow(sncrna_negative_feature_matrix)))  # Subset rows from negative data
  
  # Compute z-scores
  sncrna_functional_z_scores <- compute_robust_zscores_optimized(sncrna_positive_feature_matrix %>% dplyr::select(-Dataset), 
                                                                 sncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                                 verbose)
  sncrna_negative_z_scores <- compute_robust_zscores_optimized(sncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                               sncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                               verbose)
  
  #print(summary(sncrna_positive_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #count(sncrna_positive_feature_matrix)
  #print(summary(sncrna_negative_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #count(sncrna_negative_feature_matrix)
  
  #print(summary(sncrna_functional_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #print(summary(sncrna_negative_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  
  
  #Restore Dataset
  sncrna_functional_z_scores$Dataset <- (sncrna_positive_feature_matrix %>% 
                                           filter(Dataset ==  "short-ncrna" | Dataset == "short-ncrna") %>% 
                                           dplyr::select(Dataset))$Dataset
  sncrna_negative_z_scores$Dataset <- (sncrna_negative_feature_matrix %>% 
                                         filter(Dataset == "short-ncrna-negative-control" | Dataset == "short-ncrna-negative-control") %>% 
                                         dplyr::select(Dataset))$Dataset
  # Save z-scores to CSV
  z_scores_file <- file.path(output_dir, "sncrna_z_scores.csv")
  write.csv(sncrna_functional_z_scores, z_scores_file, row.names = FALSE)
  
  z_scores_file <- file.path(output_dir, "sncrna_negative_control_z_scores.csv")
  write.csv(sncrna_negative_z_scores, z_scores_file, row.names = FALSE)
  
  # Save raw subset to CSV
  subset_file <- file.path(output_dir, "sncrna_features.csv")
  write.csv(sncrna_positive_feature_matrix, subset_file, row.names = FALSE)
  
  subset_file <- file.path(output_dir, "sncrna_negative_control_features.csv")
  write.csv(sncrna_negative_feature_matrix, subset_file, row.names = FALSE)
}

# Compute gene functionality z-scores and save to disk:
get_robust_zscores(feature_matrix, ZSCORES_DIR)

# Compute epigenetic z-scores and save to disk:
get_robust_zscores(feature_matrix_epigenetic, EPIGENETIC_ZSCORES_DIR)
