library(dplyr)
install.packages("corrr")
library('corrr')
install.packages("ggcorrplot")
library(ggcorrplot)
install.packages("FactoMineR")
library("FactoMineR")


#### Spearman Correlation ANALYSIS ####

# Load Data #
file1a <- data.frame(Dataset = "protein-coding", read.csv("../data/features/functional-protein-exon2-dataset-features-augmented.csv", header=TRUE))
file1b <- data.frame(Dataset = "protein-coding",read.csv("../data/features/functional-protein-exon3-dataset-features-augmented.csv", header=TRUE))
file2a <- data.frame(Dataset = "lncrna",read.csv("../data/features/functional-lncrna-exon1-dataset-features-augmented.csv", header=TRUE))
file2b <- data.frame(Dataset = "lncrna",read.csv("../data/features/functional-lncrna-exon2-dataset-features-augmented.csv", header=TRUE))
file3a <- data.frame(Dataset = "short-ncrna",read.csv("../data/features/functional-short-ncrna-dataset-features-augmented.csv", header=TRUE))

filen1a <- data.frame(Dataset = "negative-control",read.csv("../data/features/protein-exon2-negative-control-dataset-features-augmented.csv", header=TRUE))
filen1b <- data.frame(Dataset = "negative-control",read.csv("../data/features/protein-exon3-negative-control-dataset-features-augmented.csv", header=TRUE))
filen2a <- data.frame(Dataset = "negative-control",read.csv("../data/features/lncrna-exon1-negative-control-dataset-features-augmented.csv", header=TRUE))
filen2b <- data.frame(Dataset = "negative-control",read.csv("../data/features/lncrna-exon2-negative-control-dataset-features-augmented.csv", header=TRUE))
filen3a <- data.frame(Dataset = "negative-control",read.csv("../data/features/short-ncrna-negative-control-dataset-features-augmented.csv", header=TRUE))

# Drop Distance column from negative controls #
filen1a <- filen1a[, -which(names(filen1a) == "Distance")]
filen1b <- filen1b[, -which(names(filen1b) == "Distance")]
filen2a <- filen2a[, -which(names(filen2a) == "Distance")]
filen2b <- filen2b[, -which(names(filen2b) == "Distance")]
filen3a <- filen3a[, -which(names(filen3a) == "Distance")]

# Load dinucleotide frequencies feature
din_file1a <- data.frame(read.csv("../data/datasets/dinucleotide_feature/protein-exon2-dinucleotide-features.csv", header=TRUE))
din_file1b <- data.frame(read.csv("../data/datasets/dinucleotide_feature/protein-exon3-dinucleotide-features.csv", header=TRUE))
din_file2a <- data.frame(read.csv("../data/datasets/dinucleotide_feature/lncrna-exon1-dinucleotide-features.csv", header=TRUE))
din_file2b <- data.frame(read.csv("../data/datasets/dinucleotide_feature/lncrna-exon2-dinucleotide-features.csv", header=TRUE))
din_file3a <- data.frame(read.csv("../data/datasets/dinucleotide_feature/short-ncrna-dinucleotide-features.csv", header=TRUE))

din_filen1a <- data.frame(read.csv("../data/datasets/dinucleotide_feature/protein-exon2-NC-dinucleotide-features.csv", header=TRUE))
din_filen1b <- data.frame(read.csv("../data/datasets/dinucleotide_feature/protein-exon3-NC-dinucleotide-features.csv", header=TRUE))
din_filen2a <- data.frame(read.csv("../data/datasets/dinucleotide_feature/lncrna-exon1-NC-dinucleotide-features.csv", header=TRUE))
din_filen2b <- data.frame(read.csv("../data/datasets/dinucleotide_feature/lncrna-exon2-NC-dinucleotide-features.csv", header=TRUE))
din_filen3a <- data.frame(read.csv("../data/datasets/dinucleotide_feature/short-ncrna-NC-dinucleotide-features.csv", header=TRUE))

# Drop Functional column from dinucleotide features #
din_file1a <- din_file1a[, -which(names(din_file1a) == "Functional")]
din_file1b <- din_file1b[, -which(names(din_file1b) == "Functional")]
din_file2a <- din_file2a[, -which(names(din_file2a) == "Functional")]
din_file2b <- din_file2b[, -which(names(din_file2b) == "Functional")]
din_file3a <- din_file3a[, -which(names(din_file3a) == "Functional")]
din_filen1a <- din_filen1a[, -which(names(din_filen1a) == "Functional")]
din_filen1b <- din_filen1b[, -which(names(din_filen1b) == "Functional")]
din_filen2a <- din_filen2a[, -which(names(din_filen2a) == "Functional")]
din_filen2b <- din_filen2b[, -which(names(din_filen2b) == "Functional")]
din_filen3a <- din_filen3a[, -which(names(din_filen3a) == "Functional")]

# Join dinucleotide features to other features #
file1a <- cbind(file1a, din_file1a)
file1b <- cbind(file1b, din_file1b)
file2a <- cbind(file2a, din_file2a)
file2b <- cbind(file2b, din_file2b)
file3a <- cbind(file3a, din_file3a)
filen1a <- cbind(filen1a, din_filen1a)
filen1b <- cbind(filen1b, din_filen1b)
filen2a <- cbind(filen2a, din_filen2a)
filen2b <- cbind(filen2b, din_filen2b)
filen3a <- cbind(filen3a, din_filen3a)

# Join data sets #
matrix_protein <- rbind(file1a, file1b, filen1a, filen1b)
matrix_lncrna <- rbind(file2a, file2b, filen2a, filen2b)
matrix_sncrna <- rbind(file3a,filen3a)
#matrix_negative_control <- rbind(,filen3a)

#all_matrix <- rbind(matrix_protein, matrix_lncrna, matrix_short_ncrna, matrix_negative_control)

# Check for nulls
colSums(is.na(matrix_protein))
colSums(is.na(matrix_lncrna))
colSums(is.na(matrix_sncrna))
# Remove nulls
matrix_protein <- na.omit(matrix_protein)
matrix_lncrna <- na.omit(matrix_lncrna)
matrix_sncrna <- na.omit(matrix_sncrna)

# Store functional and database (sequence type) columns
functional_column_protein <- matrix_protein$Functional
functional_column_lncrna <- matrix_lncrna$Functional
functional_column_sncrna <- matrix_sncrna$Functional
dataset_column_protein <- matrix_protein$Dataset
dataset_column_lncrna <- matrix_lncrna$Dataset
dataset_column_sncrna <- matrix_sncrna$Dataset
# Keep only variables of interest
list_of_features <- c("Start","copy_number","GC.","CG","GA","GG","TA","TC","X241w_PP_max",
                      "X100w_PP_max","RPKM_tissue","RPKM_primary.cell","Dfam_sum","RNAcode_score",
                      "Fickett_score","Max_covariance","MFE","Accessibility","RNAalifold_score",
                      "Interaction_ave","gnomAD_SNP_density","gnomAD_aveMAF","H3K36me3_AvgSignal",
                      "H3K27ac_AvgSignal","H3K79me2_AvgSignal","chrm_acc_AvgSignal","methylome")
feature_matrix_protein <- matrix_protein[, names(matrix_protein) %in% list_of_features]
feature_matrix_protein_numeric <- feature_matrix_protein[, sapply(feature_matrix_protein, is.numeric)]
functional_protein_numeric <- as.numeric(factor(functional_column_protein)) - 1

feature_matrix_lncrna <- matrix_lncrna[, names(matrix_lncrna) %in% list_of_features]
feature_matrix_lncrna_numeric <- feature_matrix_lncrna[, sapply(feature_matrix_lncrna, is.numeric)]
functional_lncrna_numeric <- as.numeric(factor(functional_column_lncrna)) - 1

feature_matrix_sncrna <- matrix_sncrna[, names(matrix_sncrna) %in% list_of_features]
feature_matrix_sncrna_numeric <- feature_matrix_sncrna[, sapply(feature_matrix_sncrna, is.numeric)]
functional_sncrna_numeric <- as.numeric(factor(functional_column_sncrna)) - 1


### Spearman corr calculation ###
install.packages("RVAideMemoire")
library(RVAideMemoire)
install.packages("randomForest")
library(randomForest)
# Function to calculate correlation against function from dataframe
spearman_correlation <- function(functional_numeric, numeric_data) {

  # Calculate the Spearman correlation between "functional" column and each numeric column
  correlations <- sapply(numeric_data, function(col) {
    spearman <- spearman.ci(functional_numeric, col, conf.level = 0.95, nrep = 1000)
    cor_matrix <- cor.test(functional_numeric, col, method = "spearman", alpha = 0.05)
    rho <- spearman$estimate[[1]]
    ci_inf <- spearman$conf.int[[1]] 
    ci_sup <- spearman$conf.int[[2]]
    pval <- cor_matrix$p.value
    return(c( cor = rho, ci_inf = ci_inf, ci_sup = ci_sup, pval = pval))
  })
  
  # Return correlation coefficients and confidence intervals as a named list with variable name as the name
  return(as.data.frame(correlations))
}

# Run correlation calculation
corr_protein_df <- spearman_correlation(functional_protein_numeric, feature_matrix_protein_numeric)
corr_lncrna_df <- spearman_correlation(functional_lncrna_numeric, feature_matrix_lncrna_numeric)
corr_sncrna_df <- spearman_correlation(functional_sncrna_numeric, feature_matrix_sncrna_numeric)

# Prepare datasets for ploting #
corr_protein_df <- as.data.frame(t(corr_protein_df))
corr_protein_df$feature_name <- rownames(corr_protein_df)
corr_protein_df$sequence_type <- "protein"
rownames(corr_protein_df) <- NULL

corr_lncrna_df <- as.data.frame(t(corr_lncrna_df))
corr_lncrna_df$feature_name <- rownames(corr_lncrna_df)
corr_lncrna_df$sequence_type <- "lncrna"
rownames(corr_lncrna_df) <- NULL

corr_sncrna_df <- as.data.frame(t(corr_sncrna_df))
corr_sncrna_df$feature_name <- rownames(corr_sncrna_df)
corr_sncrna_df$sequence_type <- "short_ncrna"
rownames(corr_sncrna_df) <- NULL

# Create dataset for ploting
df <- rbind(corr_protein_df, corr_lncrna_df, corr_sncrna_df)

sorted_df <- df[order(-df$cor), ]
# Create a factor for feature_type for colors
#df$feature_name <- c(list_of_features)
sorted_df$feature_name <- factor(sorted_df$feature_name, 
                          levels = unique(sorted_df$feature_name))

ggplot(sorted_df, aes(x = feature_name, y = cor, color = feature_name, shape = sequence_type, alpha = ifelse(pval > 0.05, 0.4, 1))) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ci_inf, ymax = ci_sup), width = 0.05, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(19, 15, 17)) +  #c(19, 17, 15))
  scale_color_manual(values = c("blue","green","pink4","orange4","brown","purple","red","darkgreen","darkorchid4",
                                "gold","lightcoral","lightblue4","orange3","maroon4","purple3","khaki4","limegreen","hotpink4",
                                "gold4","forestgreen","firebrick","dodgerblue4","deeppink4","darkviolet","darkslategrey","darksalmon","darkolivegreen4","darkcyan",
                                "yellowgreen","wheat4","violetred4","turquoise4","tomato4","chartreuse4","cadetblue4")) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),  # Rotate x-axis labels
        legend.position = "right",
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(title = "Spearman correlation Analysis of Selected Features",
       x = "",
       y = "Spearman correlation",
       color = "Feature",
       shape = "Sequence Type") + 
  ylim(-1, 1) + 
  scale_alpha_continuous(range = c(0.4, 1),
                         guide = "none")
