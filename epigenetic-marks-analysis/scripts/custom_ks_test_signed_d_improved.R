# Install necessary packages
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
install.packages("randomForest")
library(randomForest)

# Load required library
library(dplyr)

# Define file paths
input_file_path <- "../data/datasets/gene_functionality_features.csv"
output_file_path <- "../data/features/gene_functionality_features.csv"

# Read data
data <- read.csv(input_file_path, header = TRUE, sep = ",")

# Extract relevant columns
functional_column <- data$Functional
dataset_column <- data$Dataset

# Select features
select_features <- c("Random_number", "GC_content",
                     "CpG","GA","GG","TA","TC",
                     "low_complexity_density","phyloP max_241w","phyloP max_100w",
                     "RPKM_tissue","RPKM_primary cell","Copy number","Repeat free","RNAcode","Max covariance",
                     "MFE","RNAalifold","Interaction_ave","gnomAD_SNP_density","gnomAD_MAF",
                     "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome")

# Subset data to include only selected features
feature_matrix <- data %>%
  select(all_of(select_features))

# Add 'Functional' and 'Dataset' columns
feature_matrix <- feature_matrix %>%
  mutate(
    Functional = factor(functional_column, levels = c("No", "Yes")),
    Dataset = dataset_column
  )

# Convert 'Functional' to numeric (0 for 'No', 1 for 'Yes')
feature_matrix$Functional <- as.numeric(feature_matrix$Functional) - 1

# Write the feature matrix to a CSV file
write.csv(feature_matrix, output_file_path, row.names = FALSE)

# Create negative control and positive data sets
negative_control <- feature_matrix %>% filter(Functional == 0)
positive_data <- feature_matrix %>% filter(Functional == 1)

# Function to perform K-S test
run_ks_tests <- function(dataN, dataP) {
  n <- ncol(dataN)
  results <- matrix(nrow = 4, ncol = n)
  colnames(results) <- colnames(dataN)
  rownames(results) <- c("signed_D", "max", "min", "p.val")
  
  for (i in colnames(dataN)) {
    positive_col <- dataP[[i]]
    negative_col <- dataN[[i]]
    
    ks_test <- custom_ks_test(negative_col, positive_col)
    ks_test_p <- ks.test(negative_col, positive_col)
    
    results[1, i] <- ks_test$signed_D
    results[2, i] <- ks_test$max_diff
    results[3, i] <- ks_test$min_diff
    results[4, i] <- ks_test_p$p.value
  }
  return(results)
}

# Custom K-S test function
custom_ks_test <- function(x, y) {
  # Sort data
  x <- sort(x)
  y <- sort(y)
  
  # Calculate empirical cumulative distribution functions (ECDF)
  ecdf_x <- ecdf(x)
  ecdf_y <- ecdf(y)
  
  # Get the unique values from both samples
  unique_vals <- sort(unique(c(x, y)))
  
  # Calculate the raw differences between ECDFs
  diffs <- ecdf_x(unique_vals) - ecdf_y(unique_vals)
  
  # Find the maximum difference (positive or negative)
  max_diff <- max(diffs)
  min_diff <- min(diffs)
  
  signed_D <- if (max_diff > abs(min_diff)) max_diff else min_diff
  
  # Return the maximum and minimum differences
  return(list(signed_D = signed_D, max_diff = max_diff, min_diff = min_diff))
}

# Run the K-S tests
ks_results <- run_ks_tests(negative_control, positive_data)

# Print the results
print(ks_results)
