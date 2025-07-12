# Spearman Correlation Heatmap plot
# Load necessary libraries
library(corrr)
library(Hmisc) # For rcorr() function
library(ggcorrplot)
library(dplyr)
library(reshape2)

# Load z-scores data
source("load_gene_functionality_zscores.R")

# load Utility functions like feature renaming
source("utils.R")

# Set features to analyze
HEATMAP_SELECT_FEATURES <- c("Random number","GC content", "low_complexity_density", 
                             "CpG", "GA", "GG", "GT", "TA", "AC", "CC",
                             "phyloP max_241w", "phyloP max_100w", 
                             "GERP_91_mammals_max", "GERP_63_amniotes_max", 
                             "RPKM_tissue", "RPKM_primary cell", 
                             "H3K27ac", "H3K36me3", "H3K79me2", 
                             "chromatin_acc", "methylome", 
                             "Repeat free", "Copy number", "RNAcode", 
                             "Fickett_score", "Max covariance", "MFE",
                             "Interaction_ave", "gnomAD_SNP_density", "gnomAD_MAF")


# Set features to analyze
colnames(feature_matrix)
HEATMAP_SELECT_FEATURES_OLD_NAMES <- c("Random","GC_percentage", "lowComplexity_density", 
                             "CpG", "GA", "GG", "GT", "TA", "AC", "CC",
                             "phyloP_max_241w", "phyloP_max_100w", 
                             "GERP_91_mammals_max", "GERP_63_amniotes_max", 
                             "RPKM_tissue", "RPKM_primary.cell", 
                             "H3K9ac_MaxScaledSignal", "H3K79me2_MaxScaledSignal", 
                             "chrm_acc_MaxScaledSignal", "methylome", 
                             "repeat_distance", "copy_number", "coding_potential", 
                             "fickett", "Max_covariance", "MFE",
                             "Interaction_ave", "SNP_density", "MAF_avg")


HEATMAP_SELECT_FEATURES_LABELS <- c("Random number", "GC%", "Complexity",
                                    "CpG",  "GA", "GG", "GT", "TA", "AC", "CC",
                                    "PhyloP-mammals", "PhyloP-vertebrates",
                                    "GERP-mammals", "GERP-vertebrates",
                                    "Tissue RPKM", "Primary cell RPKM", 
                                    "H3K9ac", "H3K79me2",
                                    "Chromatin Accessibility", "Methylome",
                                    "Repeat free", "Copies", "RNAcode",
                                    "Fickett", "Covariance", "MFE",
                                    "Interactions", "SNPs", "MAF")




help(rcorr)
# Spearman Correlation Heat map
unique(feature_matrix$Dataset)
prot_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3",
                            "protein-exon2-negative-control", "protein-exon3-negative-control")) %>%
      dplyr::select(all_of(HEATMAP_SELECT_FEATURES_OLD_NAMES))
      #select(all_of(SELECTED_HISTONE_FEATURES))
    ),
  type = "spearman"
)


lncrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2",
                            "lncrna-exon1-negative-control", "lncrna-exon2-negative-control")) %>%
      dplyr::select(all_of(HEATMAP_SELECT_FEATURES_OLD_NAMES))
      #select(all_of(SELECTED_HISTONE_FEATURES))
  ),
  type = "spearman"
)

sncrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset %in% c("short-ncrna", "short-ncrna-negative-control")) %>%
      dplyr::select(all_of(HEATMAP_SELECT_FEATURES_OLD_NAMES))
      #dplyr::select(all_of(SELECTED_HISTONE_FEATURES))
  ),
  type = "spearman"
)


data_numeric_histones <- feature_matrix_histones %>%
  dplyr::select(-Dataset, ends_with("_MaxScaledSignal")) %>% 
  #rename_with(~ sub("_MaxScaledSignal$", "", .x), .cols = ends_with("_MaxScaledSignal"))
  mutate(across(everything(), .fns = ~as.numeric(as.character(.)))) %>%
  as.data.frame()
summary(data_numeric_histones)
# Restore Dataset column
data_numeric_histones$Dataset <- feature_matrix_histones$Dataset
allrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric_histones %>%
      filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3",
                            "lncrna-exon1", "lncrna-exon2",
                            "short-ncrna")) %>%
      dplyr::select(ends_with("_MaxScaledSignal"))
      #select(all_of(SELECTED_HISTONE_FEATURES))
  ),
  type = "spearman"
)

#prot_corr_obj <- rcorr(as.matrix(protein_positive_feature_matrix), type = "spearman")
#lncrna_corr_obj <- rcorr(as.matrix(lncrna_positive_feature_matrix), type = "spearman")
#sncrna_corr_obj <- rcorr(as.matrix(sncrna_positive_feature_matrix), type = "spearman")

# Protein coding #
prot_corr_matrix <- prot_corr_obj[["r"]]
prot_corr_matrix[is.na(prot_corr_matrix)] <- NA
prot_pvalues <- prot_corr_obj[["P"]]

# Number of features
n <- ncol(prot_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
prot_p_values_adjusted <- p.adjust(prot_pvalues, method="BH")

prot_p_values_adjusted <- matrix(prot_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(prot_p_values_adjusted) <- diag

rownames(prot_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
#rownames(prot_corr_matrix) <- SELECTED_HISTONE_LABELS
colnames(prot_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
#colnames(prot_corr_matrix) <- SELECTED_HISTONE_LABELS

rownames(prot_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(prot_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
#rownames(prot_p_values_adjusted) <- SELECTED_HISTONE_LABELS
#colnames(prot_p_values_adjusted) <- SELECTED_HISTONE_LABELS

help(ggcorrplot)
ggcorrplot(prot_corr_matrix, type = "lower", method = "square", lab = TRUE, 
           lab_col = "black", lab_size = 3, ggtheme = theme_void, 
           title = "Spearman correlation Heatmap - mRNA", show.diag = TRUE, 
           p.mat = prot_p_values_adjusted, sig.level = corrected_alpha, insig = "pch"
) +
  theme(
    plot.title = element_text(size = 24),  # Increase title size
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    
  )

ggsave("spearman_corr_protein_p&n.png",path = "../results/latest1000all/spearman/", scale = 3,  width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

write.csv(prot_corr_matrix, "../data/spearman_correlation/protein-correlation-matrix_p&n.csv",)
write.csv(prot_p_values_adjusted, "../data/spearman_correlation/protein-correlation-p-adj_p&n.csv")

protein_corr_matrix_melted <- melt(prot_corr_matrix)
protein_p_values_adjusted_melted <- melt(prot_p_values_adjusted)

write.csv(protein_corr_matrix_melted, "../data/spearman_correlation/protein-correlation-matrix_melted_p&n.csv")
write.csv(protein_p_values_adjusted_melted, "../data/spearman_correlation/protein-correlation-p-adj_melt_p&n.csv")


# lncRNA #
lncrna_corr_matrix <- lncrna_corr_obj[["r"]]
lncrna_pvalues <- lncrna_corr_obj[["P"]]

# Number of features
n <- ncol(lncrna_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
lncrna_p_values_adjusted <- p.adjust(lncrna_pvalues, method="BH")

lncrna_p_values_adjusted <- matrix(lncrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(lncrna_p_values_adjusted) <- diag

rownames(lncrna_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(lncrna_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
#rownames(lncrna_corr_matrix) <- SELECTED_HISTONE_LABELS
#colnames(lncrna_corr_matrix) <- SELECTED_HISTONE_LABELS

rownames(lncrna_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(lncrna_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
#rownames(lncrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS
#colnames(lncrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS

ggcorrplot(lncrna_corr_matrix, type = "lower", method = "square", lab = TRUE, 
           lab_col = "black", lab_size = 3, ggtheme = theme_void, 
           title = "Spearman correlation Heatmap - lncRNA", show.diag = TRUE, 
           p.mat = lncrna_p_values_adjusted, sig.level = corrected_alpha, insig = "pch"
) +
  theme(
    plot.title = element_text(size = 24),  # Increase title size
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    
  )
ggsave("spearman_corr_lncrna_p&n.png",path = "../results/latest1000all/spearman/", scale = 3,  width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

write.csv(lncrna_corr_matrix, "../data/spearman_correlation/lncrna-correlation-matrix_p&n.csv")
write.csv(lncrna_p_values_adjusted, "../data/spearman_correlation/lncrna-correlation-p-adj_p&n.csv")

lncrna_corr_matrix_melted <- melt(lncrna_corr_matrix)
lncrna_p_values_adjusted_melted <- melt(lncrna_p_values_adjusted)

write.csv(lncrna_corr_matrix_melted, "../data/spearman_correlation/lncrna-correlation-matrix_melted_p&n.csv")
write.csv(lncrna_p_values_adjusted_melted, "../data/spearman_correlation/lncrna-correlation-p-adj_melt_p&n.csv")


# sncRNA #
sncrna_corr_matrix <- sncrna_corr_obj[["r"]]
sncrna_pvalues <- sncrna_corr_obj[["P"]]

# Number of features
n <- ncol(sncrna_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
sncrna_p_values_adjusted <- p.adjust(sncrna_pvalues, method="BH")

sncrna_p_values_adjusted <- matrix(sncrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(sncrna_p_values_adjusted) <- diag

rownames(sncrna_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(sncrna_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
#rownames(sncrna_corr_matrix) <- SELECTED_HISTONE_LABELS
#colnames(sncrna_corr_matrix) <- SELECTED_HISTONE_LABELS

rownames(sncrna_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(sncrna_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
#rownames(sncrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS
#colnames(sncrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS

ggcorrplot(sncrna_corr_matrix, type = "lower", method = "square", lab = TRUE, 
           lab_col = "black", lab_size = 3, ggtheme = theme_void, 
           title = "Spearman correlation Heatmap - sncRNA", show.diag = TRUE, 
           p.mat = sncrna_p_values_adjusted, sig.level = corrected_alpha, insig = "pch"
) +
  theme(
    plot.title = element_text(size = 24),  # Increase title size
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    
  )
ggsave("spearman_corr_sncrna_p&n.png",path = "../results/latest1000all/spearman/", scale = 3,  width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

write.csv(sncrna_corr_matrix, "../data/spearman_correlation/sncrna-correlation-matrix_p&n.csv")
write.csv(sncrna_p_values_adjusted, "../data/spearman_correlation/sncrna-correlation-p-adj_p&n.csv")

sncrna_corr_matrix_melted <- melt(sncrna_corr_matrix)
sncrna_p_values_adjusted_melted <- melt(sncrna_p_values_adjusted)

write.csv(sncrna_corr_matrix_melted, "../data/spearman_correlation/sncrna-correlation-matrix_melted_p&n.csv")
write.csv(sncrna_p_values_adjusted_melted, "../data/spearman_correlation/sncrna-correlation-p-adj_melt_p&n.csv")



##################
# allRNA #
allrna_corr_matrix <- allrna_corr_obj[["r"]]
allrna_corr_matrix[is.na(allrna_corr_matrix)] <- NA
allrna_pvalues <- allrna_corr_obj[["P"]]


# Number of features
n <- ncol(allrna_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
allrna_p_values_adjusted <- p.adjust(allrna_pvalues, method="BH")

allrna_p_values_adjusted <- matrix(allrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(allrna_p_values_adjusted) <- diag

#rownames(prot_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
SELECTED_HISTONE_LABELS <- colnames(data_numeric_histones %>% dplyr::select(ends_with("_MaxScaledSignal")) %>% rename_with(~ sub("_MaxScaledSignal$", "", .x), .cols = ends_with("_MaxScaledSignal")))
rownames(allrna_corr_matrix) <- SELECTED_HISTONE_LABELS
#colnames(prot_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(allrna_corr_matrix) <- SELECTED_HISTONE_LABELS

#rownames(prot_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
#colnames(prot_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
rownames(allrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS
colnames(allrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS

help(ggcorrplot)
ggcorrplot(allrna_corr_matrix, type = "lower", method = "square", lab = TRUE, 
           lab_col = "black", lab_size = 3, ggtheme = theme_void(), 
           title = "Spearman correlation Heatmap - All RNA types", show.diag = TRUE, 
           p.mat = prot_p_values_adjusted, sig.level = corrected_alpha, insig = "pch"
) +
  theme(
    plot.title = element_text(size = 24),  # Increase title size
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    
  )
help(heatmap)
heatmap(allrna_corr_matrix)
ggsave("spearman_corr_allrna_plus.png",path = "../results/latest1000all/spearman/histones/", scale = 3,  width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

write.csv(prot_corr_matrix, "../data/spearman_correlation/protein-correlation-matrix_plus_histones.csv",)
write.csv(prot_p_values_adjusted, "../data/spearman_correlation/protein-correlation-p-adj_plus_histones.csv")

protein_corr_matrix_melted <- melt(prot_corr_matrix)
protein_p_values_adjusted_melted <- melt(prot_p_values_adjusted)

write.csv(protein_corr_matrix_melted, "../data/spearman_correlation/protein-correlation-matrix_melted_histones.csv")
write.csv(protein_p_values_adjusted_melted, "../data/spearman_correlation/protein-correlation-p-adj_melt_histones.csv")

################

########################################################
#Obtain correlation between all features
# Spearman Correlation Heat map
unique(feature_matrix$Dataset)
prot_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset == "protein-coding-exon2" | Dataset == "protein-coding-exon3") %>%
      select(-Dataset)
  ),
  type = "spearman"
)


lncrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset == "lncrna-exon1" | Dataset == "lncrna-exon2") %>%
      select(-Dataset)
  ),
  type = "spearman"
)

sncrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset == "short-ncrna") %>%
      select(-Dataset)
  ),
  type = "spearman"
)

# Protein coding #
prot_corr_matrix <- prot_corr_obj[["r"]]
prot_pvalues <- prot_corr_obj[["P"]]

# Number of features
n <- ncol(prot_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
prot_p_values_adjusted <- p.adjust(prot_pvalues, method="BH")

prot_p_values_adjusted <- matrix(prot_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(prot_p_values_adjusted) <- diag

# Melt
protein_corr_matrix_melted <- melt(prot_corr_matrix)
protein_p_values_adjusted_melted <- melt(prot_p_values_adjusted)

write.csv(protein_corr_matrix_melted, "../data/spearman_correlation/protein-correlation-matrix_melted.csv")
write.csv(protein_p_values_adjusted_melted, "../data/spearman_correlation/protein-correlation-p-adj_melt.csv")


# lncRNA #
lncrna_corr_matrix <- lncrna_corr_obj[["r"]]
lncrna_pvalues <- lncrna_corr_obj[["P"]]

# Number of features
n <- ncol(lncrna_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
lncrna_p_values_adjusted <- p.adjust(lncrna_pvalues, method="BH")

lncrna_p_values_adjusted <- matrix(lncrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(lncrna_p_values_adjusted) <- diag

#Melt
lncrna_corr_matrix_melted <- melt(lncrna_corr_matrix)
lncrna_p_values_adjusted_melted <- melt(lncrna_p_values_adjusted)

write.csv(lncrna_corr_matrix_melted, "../data/spearman_correlation/lncrna-correlation-matrix_melted.csv")
write.csv(lncrna_p_values_adjusted_melted, "../data/spearman_correlation/lncrna-correlation-p-adj_melt.csv")


# sncRNA #
sncrna_corr_matrix <- sncrna_corr_obj[["r"]]
sncrna_pvalues <- sncrna_corr_obj[["P"]]

# Number of features
n <- ncol(sncrna_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
sncrna_p_values_adjusted <- p.adjust(sncrna_pvalues, method="BH")

sncrna_p_values_adjusted <- matrix(sncrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(sncrna_p_values_adjusted) <- diag

sncrna_corr_matrix_melted <- melt(sncrna_corr_matrix)
sncrna_p_values_adjusted_melted <- melt(sncrna_p_values_adjusted)

write.csv(sncrna_corr_matrix_melted, "../data/spearman_correlation/sncrna-correlation-matrix_melted.csv")
write.csv(sncrna_p_values_adjusted_melted, "../data/spearman_correlation/sncrna-correlation-p-adj_melt.csv")

