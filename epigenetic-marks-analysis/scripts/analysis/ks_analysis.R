# KS test analysis
library(randomForest)

# Load data
source("load_gene_functionality_features.R")

# Define features to analyze
ks_select_features <- c("GC content", "low_complexity_density", 
                        "CpG", "GA", "GG", "GT", "TA", "AC", "CC",
                        "phyloP max_241w", "phyloP max_100w", 
                        "GERP_91_mammals_max", "GERP_63_amniotes_max", 
                        "RPKM_tissue", "RPKM_primary cell", 
                        "H3K27ac", "H3K36me3", "H3K79me2", 
                        "chromatin_acc", "methylome", 
                        "Repeat free", "Copy number", "RNAcode", 
                        "Fickett_score", "Max covariance", "MFE",
                        "Interaction_ave", "gnomAD_SNP_density", "gnomAD_MAF"
)

# Select just desired features
funcProtDataSelect <- funcProtExon2Data %>% rbind(funcProtExon3Data) %>% select(all_of(ks_select_features))
funcLncrnaDataSelect <- funcLncrnaExon1Data %>% rbind(funcLncrnaExon2Data) %>% select(all_of(ks_select_features))
funcSncrnaDatasetSelect <- funcSncrnaDataset %>% select(all_of(ks_select_features))

protNCDataSelect <- protExon2NCData %>% rbind(protExon3NCData) %>% select(all_of(ks_select_features))
lncrnaNCDataSelect <- lncrnaExon1NCData %>% rbind(lncrnaExon2NCData) %>% select(all_of(ks_select_features))
sncrnaNCDataSelect <- sncrnaNCData %>% select(all_of(ks_select_features))

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

# Compute KS stats:
ks_results_prot <- run_ks_tests(protNCDataSelect, funcProtDataSelect)
ks_results_lncrna <- run_ks_tests(lncrnaNCDataSelect, funcLncrnaDataSelect)
ks_results_sncrna <- run_ks_tests(sncrnaNCDataSelect, funcSncrnaDatasetSelect)

# Use BH/FDR correction to account for multiple tests.
print(ks_results_prot["p.val",])

fdrs<-p.adjust(ks_results_prot["p.val",], method="BH")
print(fdrs)
ks_results_prot["p.val",] <- fdrs

print(ks_results_lncrna["p.val",])
fdrs<-p.adjust(ks_results_lncrna["p.val",], method="BH")
print(fdrs)
ks_results_lncrna["p.val",] <- fdrs

print(ks_results_sncrna["p.val",])
fdrs<-p.adjust(ks_results_sncrna["p.val",], method="BH")
print(fdrs)
ks_results_sncrna["p.val",] <- fdrs

write.csv(ks_results_prot,"../results/ks_results_prot.csv", row.names = TRUE)
write.csv(ks_results_lncrna,"../results/ks_results_lncrna.csv", row.names = TRUE)
write.csv(ks_results_sncrna,"../results/ks_results_sncrna.csv", row.names = TRUE)

