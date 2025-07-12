library(randomForest)

all_data_file <- "../data/datasets/gene_functionality_features.csv"

data <- read.csv(all_data_file, header = TRUE, sep = ",")

functional_column <- data$Functional
dataset_column <- data$Dataset
colnames(feature_matrix)
select_features <- c("GC_percentage","GA","CpG","GG",
                     "TA","phyloP_max_241w","phyloP_max_100w","RPKM_tissue","RPKM_primary.cell",
                     "copy_number","repeat_distance","Interaction_ave","coding_potential","RNAalifold_score",
                     "Max_covariance","MFE","SNP_density","MAF_avg","H3K27ac_AvgSignal",
                     "H3K36me3_AvgSignal","H3K79me2_AvgSignal","chrm_acc_AvgSignal","Random","accessibility",
                     "fickett","GERP_91_mammals_max","GERP_63_amniotes_max","AA","AC",
                     "AG","AT","CA","CC","CT",
                     "GC","GT","TC","TG","TT")
feature_matrix <- data[, names(data) %in% select_features]
feature_matrix$Functional <- functional_column
feature_matrix_numeric <- feature_matrix[,names(feature_matrix) %in% c("accessibility","H3K79me2")]


feature_matrix_numeric$Functional <- factor(feature_matrix_numeric$Functional, levels=c("No","Yes"))
functional_numeric <- as.numeric(factor(feature_matrix_numeric$Functional)) - 1
feature_matrix_numeric$Functional <- functional_numeric
feature_matrix$Dataset <- dataset_column
feature_matrix$GA <- data$GA
feature_matrix$CpG <- data$CpG
feature_matrix$TA <- data$TA
feature_matrix$TC <- data$TC
feature_matrix$GG <- data$GG
feature_matrix_numeric_1 <- feature_matrix_numeric[, sapply(feature_matrix_numeric, is.numeric)]
feature_matrix_numeric_imputed <- feature_matrix_numeric_1
summary(feature_matrix_numeric_imputed)
help("na.roughfix")

negative_control <- feature_matrix_numeric_imputed[feature_matrix_numeric$Functional == 0,]
positive_data <- feature_matrix_numeric_imputed[feature_matrix_numeric$Functional == 1,]

prot_positive_data <- feature_matrix_numeric_imputed[feature_matrix$Functional == 1 &
                                                       (feature_matrix$Dataset=="protein-coding-exon2" | feature_matrix$Dataset=="protein-coding-exon3"),]
prot_positive_data <- prot_positive_data[, select_features]

prot_negative_data <- feature_matrix_numeric_imputed[feature_matrix$Functional == 0 &
                                       (feature_matrix$Dataset=="protein-exon2-negative-control" | feature_matrix$Dataset=="protein-exon3-negative-control"),]
prot_negative_data <- prot_negative_data[, select_features]

lncrna_positive_data <- feature_matrix_numeric_imputed[feature_matrix$Functional == 1 &
                                         (feature_matrix$Dataset=="lncrna-exon1" | feature_matrix$Dataset=="lncrna-exon2"),]
lncrna_positive_data <- lncrna_positive_data[, select_features]

lncrna_negative_data <- feature_matrix_numeric_imputed[feature_matrix$Functional == 0 &
                                         (feature_matrix$Dataset=="lncrna-exon1-negative-control" | feature_matrix$Dataset=="lncrna-exon2-negative-control"),]
lncrna_negative_data <- lncrna_negative_data[, select_features]

sncrna_positive_data <- feature_matrix_numeric_imputed[feature_matrix$Functional == 1 &
                                                         feature_matrix$Dataset=="short-ncrna",]
sncrna_positive_data <- sncrna_positive_data[, select_features]

sncrna_negative_data <- feature_matrix_numeric_imputed[feature_matrix$Functional == 0 &
                                         (feature_matrix$Dataset=="short-ncrna-negative-control"),]
sncrna_negative_data <- sncrna_negative_data[, select_features]
help("ks.test")
# Function to perform K-S test
run_ks_tests <- function(data) {
  n <- ncol(data)
  results <- matrix(nrow = n, ncol = n)
  
  diag(results) <- 0 # Diagonals will be 0 (difference of a variable with itself)
  
  colnames(results) <- rownames(results) <- colnames(data)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      ks_test <- ks.test(data[, i], data[, j])
      
      results[i,j] <- ks_test$statistic
      #results[j,i] <- ks_test$max_diff
      #results[j,i] <- ks_test$min_diff
      results[j,i] <- ks_test$statistic
    }
  }
  
  return(results)
}

# Function to perform K-S test
run_ks_pval <- function(data) {
  n <- ncol(data)
  results <- matrix(nrow = n, ncol = n)
  
  diag(results) <- 0 # Diagonals will be 0 (difference of a variable with itself)
  
  colnames(results) <- rownames(results) <- colnames(data)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      ks_test_p <- ks.test(data[, i], data[, j])
      
      results[i,j] <- ks_test_p$p.value
      #results[j,i] <- ks_test$max_diff
      #results[j,i] <- ks_test$min_diff
      results[j,i] <- ks_test_p$p.value
    }
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

ks_results_prot <- run_ks_tests(prot_negative_data, prot_positive_data)
#ks_pval_prot <- run_ks_pval(rbind(prot_negative_data, prot_positive_data))

lncrna_positive_data$accessibility <- accessibility_column_lncrna$accessibility
lncrna_negative_data$accessibility <- accessibility_column_lncrna_nc$accessibility
summary(lncrna_positive_data)
lncrna_positive_data <- as.data.frame(sapply(lncrna_positive_data, convert_to_numeric))
lncrna_negative_data <- as.data.frame(sapply(lncrna_negative_data, convert_to_numeric))
ks_results_lncrna <- run_ks_tests(lncrna_negative_data, lncrna_positive_data)
#ks_pval_lncrna <- run_ks_pval(rbind(lncrna_negative_data, lncrna_positive_data))

ks_results_sncrna <- run_ks_tests(sncrna_negative_data, sncrna_positive_data)
#ks_pval_sncrna <- run_ks_pval(rbind(sncrna_negative_data, sncrna_positive_data))

# Use BH/FDR p-value correction to account for multiple tests.
print(ks_results_prot["p.val",])
help("p.adjust")
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
write.csv(ks_results_lncrna,"../results/ks_results_lincrna.csv", row.names = TRUE)
write.csv(ks_results_sncrna,"../results/ks_results_sincrna.csv", row.names = TRUE)


# Heatmap
library(ggplot2)
library(reshape2)
# Melt the correlation matrix into long format
# Reshape to suit ggplot, remove NAs, and sort the labels
ks_melted <- na.omit(melt(ks_results_d))
pval_melted <- na.omit(melt(pval_matrix))

combined_matrix <- merge(ks_melted, pval_melted, by = c("Var1","Var2")) 

heatmap(ks_results, Rowv = NA, Colv = NA)

# Inspiration from: http://pseudofish.com/trianglTRUE# Inspiration from: http://pseudofish.com/triangle-heatmaps-in-r-using-ggplot.html
# Triangle heatmap
ggplot(combined_matrix, aes(Var1, Var2)) +
  ggtitle('Selected Features K-S Test - short ncrna') +
  theme_bw() +
  geom_tile(aes(fill = value.x), color='white') +
  geom_text(aes(Var2, Var1, label = ifelse(value.y < corrected_alpha, 
                                           round(value.x, 2), "X"), 
                                           color = value.x), 
            nudge_x = 0, nudge_y = 0, fontface = "bold") + 
  scale_color_gradient2(low = "dodgerblue", mid = "bisque", high = "darkorchid4", midpoint = 0.5, limits = c(0, 1)) +
  scale_fill_gradient2(low = "white", high = "white", midpoint = 0, guide = "none") +  # For background
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='white'),
        axis.text.y = element_text(hjust = 1) ) +
  labs(x = "", y = "")


