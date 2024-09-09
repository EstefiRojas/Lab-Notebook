# Declare required libraries
library(magrittr)
library(dplyr)
library(ggpubr)

# Load data
source("load_gene_functionality_features.R", chdir = TRUE)

# Define features to analyze
MDS_SELECT_FEATURES <- c("Random number",
                         "GC content",
                         "CpG", "GA", "GG", "TA", "TC",
                         "low_complexity_density", "Conservation",
                         "phyloP max_100w", "Expression", "RPKM_primary cell",
                         "Copy number", "Repeat free", "RNAcode",
                         "RNA Covariance", "MFE", "RNAalifold",
                         "RNA-RNA interactions", "SNP density", "gnomAD_MAF",
                         "H3K27ac", "H3K36me3", "H3K79me2", "chromatin_acc",
                         "methylome")
feature_matrix <- na.omit(feature_matrix)
feature_matrix_numeric_no_nas <-
  get_numeric_features(MDS_SELECT_FEATURES)$feature_matrix_numeric_no_nas


functional_z_scores <- as.data.frame(get_mod_zscores_equal_cols(feature_matrix_numeric_no_nas %>% 
                                                    filter(Dataset == "protein-coding-exon2" |
                                                             Dataset == "protein-coding-exon3" |
                                                             Dataset == "lncrna-exon1" |
                                                             Dataset == "lncrna-exon2" |
                                                             Dataset == "short-ncrna"), 
                                                  feature_matrix_numeric_no_nas %>% 
                                                    filter(Dataset == "protein-exon2-negative-control" |
                                                             Dataset == "protein-exon3-negative-control" |
                                                             Dataset == "lncrna-exon1-negative-control" |
                                                             Dataset == "lncrna-exon2-negative-control" |
                                                             Dataset == "short-ncrna-negative-control"),
                                                  selected_features = MDS_SELECT_FEATURES)$zscores)
negative_control_z_scores <- as.data.frame(get_mod_zscores_equal_cols(feature_matrix_numeric_no_nas %>% 
                                                          filter(Dataset == "protein-exon2-negative-control" |
                                                                   Dataset == "protein-exon3-negative-control" |
                                                                   Dataset == "lncrna-exon1-negative-control" |
                                                                   Dataset == "lncrna-exon2-negative-control" |
                                                                   Dataset == "short-ncrna-negative-control"), 
                                                        feature_matrix_numeric_no_nas %>% 
                                                          filter(Dataset == "protein-exon2-negative-control" |
                                                                   Dataset == "protein-exon3-negative-control" |
                                                                   Dataset == "lncrna-exon1-negative-control" |
                                                                   Dataset == "lncrna-exon2-negative-control" |
                                                                   Dataset == "short-ncrna-negative-control"),
                                                        selected_features = MDS_SELECT_FEATURES)$zscores)
str(feature_matrix_numeric_no_nas[MDS_SELECT_FEATURES])
data_scaled <- scale(feature_matrix_numeric_no_nas[MDS_SELECT_FEATURES])
data_normalized <- as.data.frame(rbind(as.data.frame(functional_z_scores), as.data.frame(negative_control_z_scores)))


# MDS Analisys on z-scores
mds_nor <- data_normalized %>%
  dist() %>%
  cmdscale() %>%
  as_tibble()

colnames(mds_nor) <- c("Dim.1", "Dim.2")
mds_nor$gene_type <- factor(feature_matrix_numeric_no_nas$Dataset,
                            levels = c("protein-coding-exon2","protein-coding-exon3","lncrna-exon1",
                                       "lncrna-exon2","short-ncrna","protein-exon2-negative-control",
                                       "protein-exon3-negative-control","lncrna-exon1-negative-control","lncrna-exon2-negative-control",
                                       "short-ncrna-negative-control"),
                            labels = c("Protein coding", "Protein coding",
                                       "lncRNA", "lncRNA",
                                       "sncRNA",
                                       "Negative controls",
                                       "Negative controls",
                                       "Negative controls",
                                       "Negative controls",
                                       "Negative controls"))
mds_nor$Dim.1_scaled <- scale(mds_nor$Dim.1)
mds_nor$Dim.2_scaled <- scale(mds_nor$Dim.2)
ggscatter(mds_nor, x = "Dim.1", y = "Dim.2",
          size = 1,
          repel = TRUE,
          color = "gene_type"
          ) +
  labs(color = "Gene Type") +
  scale_color_manual(values = c("Protein coding" = "firebrick",
                                "sncRNA" = "limegreen",
                                "lncRNA" = "darkorange",
                                "Negative controls" = "blue"))

# Retry analysis with scaled raw data instead of z-scores.
mds_scale <- data_scaled %>%
  dist() %>%
  cmdscale() %>%
  as_tibble()

colnames(mds_scale) <- c("Dim.1", "Dim.2")
mds_scale$gene_type <- factor(feature_matrix_numeric_no_nas$Dataset,
                            levels = c("protein-coding-exon2","protein-coding-exon3","lncrna-exon1",
                                       "lncrna-exon2","short-ncrna","protein-exon2-negative-control",
                                       "protein-exon3-negative-control","lncrna-exon1-negative-control","lncrna-exon2-negative-control",
                                       "short-ncrna-negative-control"),
                            labels = c("Protein coding", "Protein coding",
                                       "lncRNA", "lncRNA",
                                       "sncRNA",
                                       "Negative controls",
                                       "Negative controls",
                                       "Negative controls",
                                       "Negative controls",
                                       "Negative controls"))

ggscatter(mds_scale, x = "Dim.1", y = "Dim.2",
          size = 1,
          repel = TRUE,
          color = "gene_type",
          facet.by = "gene_type") +
  labs(color = "Gene Type") +
  scale_color_manual(values = c("Protein coding" = "firebrick",
                                "sncRNA" = "limegreen",
                                "lncRNA" = "darkorange",
                                "Negative controls" = "blue")) +
  theme(legend.position = "right")

