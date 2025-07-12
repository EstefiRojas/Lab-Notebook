#load_ks_stats_utils.R

################################
# Function to perform K-S test #
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

############################
# Custom K-S test function #
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
