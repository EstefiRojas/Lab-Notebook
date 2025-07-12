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

# compute zscores for all data
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

# Compute zscores for protein coding
protein_functional_z_scores <- as.data.frame(get_mod_zscores_equal_cols(feature_matrix_numeric_no_nas %>% 
                                                                  filter(Dataset == "protein-coding-exon2" |
                                                                           Dataset == "protein-coding-exon3"), 
                                                                feature_matrix_numeric_no_nas %>% 
                                                                  filter(Dataset == "protein-exon2-negative-control" |
                                                                           Dataset == "protein-exon3-negative-control"),
                                                                selected_features = MDS_SELECT_FEATURES)$zscores)
protein_negative_control_z_scores <- as.data.frame(get_mod_zscores_equal_cols(feature_matrix_numeric_no_nas %>% 
                                                                        filter(Dataset == "protein-exon2-negative-control" |
                                                                                 Dataset == "protein-exon3-negative-control"), 
                                                                      feature_matrix_numeric_no_nas %>% 
                                                                        filter(Dataset == "protein-exon2-negative-control" |
                                                                                 Dataset == "protein-exon3-negative-control"),
                                                                      selected_features = MDS_SELECT_FEATURES)$zscores)

protein_data_normalized <- as.data.frame(rbind(as.data.frame(protein_functional_z_scores), as.data.frame(protein_negative_control_z_scores)))

# Compute zscores for lncRNA
lncrna_functional_z_scores <- as.data.frame(get_mod_zscores_equal_cols(feature_matrix_numeric_no_nas %>% 
                                                                          filter(Dataset == "lncrna-exon1" |
                                                                                   Dataset == "lncrna-exon2" ), 
                                                                        feature_matrix_numeric_no_nas %>% 
                                                                          filter(Dataset == "lncrna-exon1-negative-control" |
                                                                                   Dataset == "lncrna-exon2-negative-control"),
                                                                        selected_features = MDS_SELECT_FEATURES)$zscores)
lncrna_negative_control_z_scores <- as.data.frame(get_mod_zscores_equal_cols(feature_matrix_numeric_no_nas %>% 
                                                                                filter(Dataset == "lncrna-exon1-negative-control" |
                                                                                         Dataset == "lncrna-exon2-negative-control"), 
                                                                              feature_matrix_numeric_no_nas %>% 
                                                                                filter(Dataset == "lncrna-exon1-negative-control" |
                                                                                         Dataset == "lncrna-exon2-negative-control"),
                                                                              selected_features = MDS_SELECT_FEATURES)$zscores)

lncrna_data_normalized <- as.data.frame(rbind(as.data.frame(lncrna_functional_z_scores), as.data.frame(lncrna_negative_control_z_scores)))


# Compute zscores for sncRNA
sncrna_functional_z_scores <- as.data.frame(get_mod_zscores_equal_cols(feature_matrix_numeric_no_nas %>% 
                                                                         filter(Dataset == "short-ncrna"), 
                                                                       feature_matrix_numeric_no_nas %>% 
                                                                         filter(Dataset == "short-ncrna-negative-control"),
                                                                       selected_features = MDS_SELECT_FEATURES)$zscores)
sncrna_negative_control_z_scores <- as.data.frame(get_mod_zscores_equal_cols(feature_matrix_numeric_no_nas %>% 
                                                                               filter(Dataset == "short-ncrna-negative-control"), 
                                                                             feature_matrix_numeric_no_nas %>% 
                                                                               filter(Dataset == "short-ncrna-negative-control"),
                                                                             selected_features = MDS_SELECT_FEATURES)$zscores)

sncrna_data_normalized <- as.data.frame(rbind(as.data.frame(sncrna_functional_z_scores), as.data.frame(sncrna_negative_control_z_scores)))



# MDS Analisys on z-scores for all data
help("cmdscale")
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

# Set up the plot with custom margins and 'cex' scaling
par(mar = c(5, 5, 4, 2) + 0.1,  # Adjust margins as needed
    cex = 1.5)                   # Adjust 'cex' scaling factor as needed

# Plot the MDS results
plot(mds_result$points, 
     xlab = "Dimension 1", 
     ylab = "Dimension 2", 
     main = "MDS Plot",        # Add a title
     pch = 20,                 # Use filled circles for points
     col = "blue")             # Set point color

# Add labels to the points (optional)
text(mds_result$points, 
     labels = rownames(mds_result$points), 
     pos = 3,                 # Position labels above the points
     cex = 1.2)               # Adjust label size as needed


# MDS Analisys on z-scores for protein coding
protein_mds_nor <- protein_data_normalized %>%
  dist() %>%
  cmdscale() %>%
  as_tibble()

colnames(protein_mds_nor) <- c("Dim.1", "Dim.2")
protein_dataset_vector <- feature_matrix_numeric_no_nas %>%
  filter(Dataset=="protein-coding-exon2" | Dataset=="protein-coding-exon3" |
           Dataset=="protein-exon2-negative-control" | Dataset=="protein-exon3-negative-control")
protein_mds_nor$gene_type <- factor(protein_dataset_vector$Dataset,
                            levels = c("protein-coding-exon2",
                                       "protein-coding-exon3",
                                       "protein-exon2-negative-control",
                                       "protein-exon3-negative-control"),
                            labels = c("Protein coding", "Protein coding",
                                       "Negative controls",
                                       "Negative controls"))

ggscatter(protein_mds_nor, x = "Dim.1", y = "Dim.2",
          size = 1,
          repel = TRUE,
          color = "gene_type"
) +
  labs(color = "Gene Type") +
  scale_color_manual(values = c("Protein coding" = "firebrick",
                                "sncRNA" = "limegreen",
                                "lncRNA" = "darkorange",
                                "Negative controls" = "blue"))


# MDS Analisys on z-scores for lncRNA
lncrna_mds_nor <- lncrna_data_normalized %>%
  dist() %>%
  cmdscale() %>%
  as_tibble()

colnames(lncrna_mds_nor) <- c("Dim.1", "Dim.2")
lncrna_dataset_vector <- feature_matrix_numeric_no_nas %>%
  filter(Dataset=="lncrna-exon1" | Dataset=="lncrna-exon2" |
           Dataset=="lncrna-exon1-negative-control" | Dataset=="lncrna-exon2-negative-control")
lncrna_mds_nor$gene_type <- factor(lncrna_dataset_vector$Dataset,
                                    levels = c("lncrna-exon1", "lncrna-exon2",
                                               "lncrna-exon1-negative-control","lncrna-exon2-negative-control"),
                                    labels = c("lncRNA", "lncRNA",
                                               "Negative controls",
                                               "Negative controls"))

ggscatter(lncrna_mds_nor, x = "Dim.1", y = "Dim.2",
          size = 1,
          repel = TRUE,
          color = "gene_type"
) +
  labs(color = "Gene Type") +
  scale_color_manual(values = c("Protein coding" = "firebrick",
                                "sncRNA" = "limegreen",
                                "lncRNA" = "darkorange",
                                "Negative controls" = "blue"))

# MDS Analisys on z-scores for sncRNA
sncrna_mds_nor <- sncrna_data_normalized %>%
  dist() %>%
  cmdscale() %>%
  as_tibble()

colnames(sncrna_mds_nor) <- c("Dim.1", "Dim.2")
sncrna_dataset_vector <- feature_matrix_numeric_no_nas %>%
  filter(Dataset=="short-ncrna" |
           Dataset=="short-ncrna-negative-control")
sncrna_mds_nor$gene_type <- factor(sncrna_dataset_vector$Dataset,
                                   levels = c("short-ncrna",
                                              "short-ncrna-negative-control"),
                                   labels = c("sncRNA",
                                              "Negative controls"))

ggscatter(sncrna_mds_nor, x = "Dim.1", y = "Dim.2",
          size = 1,
          repel = TRUE,
          color = "gene_type"
) +
  labs(color = "Gene Type") +
  scale_color_manual(values = c("Protein coding" = "firebrick",
                                "sncRNA" = "limegreen",
                                "lncRNA" = "darkorange",
                                "Negative controls" = "blue"))
