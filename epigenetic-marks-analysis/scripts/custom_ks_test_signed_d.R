all_data_file <- "../data/datasets/gene_functionality_features.csv"

data <- read.csv(all_data_file, header = TRUE, sep = ",")

functional_column <- data$Functional
dataset_column <- data$Dataset
select_features <- c("GC_percent","GA","CpG","GG","TA","TC",
                     "X100w_PP_max","X241w_PP_max","RPKM_primary.cell","RPKM_tissue","copy_number","Dfam_sum",
                     "Interaction_ave","Fickett_score","RNAcode_score","RNAalifold_score","Accessibility","Max_covariance",
                     "MFE","gnomAD_aveMAF","gnomAD_SNP_density","chrm_acc_AvgSignal","H3K27ac_AvgSignal","H3K36me3_AvgSignal","H3K79me2_AvgSignal",
                     "methylome","random")
feature_matrix <- data[, names(data) %in% select_features]
feature_matrix$Functional <- functional_column
feature_matrix$Functional <- factor(feature_matrix$Functional, levels=c("No","Yes"))
functional_numeric <- as.numeric(factor(feature_matrix$Functional)) - 1
feature_matrix$Functional <- functional_numeric
feature_matrix$Dataset <- dataset_column

write.csv(feature_matrix,"../data/features/gene_functionality_features.csv", row.names = FALSE)

negative_control <- feature_matrix[feature_matrix$Functional == 0,]
positive_data <- feature_matrix[feature_matrix$Functional == 1,]

# Function to perform K-S test
run_ks_tests <- function(dataN, dataP) {
  n <- ncol(dataN)
  results <- matrix(nrow = 4, ncol = n)
  colnames(results) <- colnames(dataN)
  rownames(results) <- c("signed_D","max","min","p.val")
  
  for (i in colnames(dataN)) {
    #positive_col <- remove_outliers_IQR(dataP, i)
    #negative_col <- remove_outliers_IQR(dataN, i)
    positive_col <- dataP[[i]]
    negative_col <- dataN[[i]]
    
    ks_test <- custom_ks_test(negative_col, positive_col)
    ks_test_p <- ks.test(negative_col, positive_col)
    
    results[1,i] <- ks_test$signed_D
    results[2,i] <- ks_test$max_diff
    results[3,i] <- ks_test$min_diff
    results[4,i] <- ks_test_p$p.value
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
  
  signed_D <- 0
  if(max_diff > abs(min_diff)) {
    signed_D <- max_diff
  } else {
    signed_D <- min_diff
  }
  
  # Return the maximum and minimum differences
  return(list(signed_D = signed_D, max_diff = max_diff, min_diff = min_diff))
}

# Run the K-S tests
ks_results <- run_ks_tests(negative_control, positive_data)
ks_results
