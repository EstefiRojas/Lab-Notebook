# Declare required libraries
library(magrittr)
library(dplyr)
install.packages("ggpubr", dependencies = TRUE)
library(ggpubr)
library(MASS)
library(patchwork)

# Sourcing the following R file will load into the current session previously 
# computed robust z-scores of gene functionality features:
# 
source("load_features_robust_zscores.R")
#source("load_gene_functionality_features.R", chdir = TRUE)

colnames(zscores_all)
# Define features to analyze
PCA_12_SELECT_FEATURES <- c("GC_percentage",
                         "CpG",
                         "phyloP_max_241w",
                         "RPKM_tissue",
                         "copy_number", "coding_potential",
                         "Max_covariance", "MFE",
                         "methylome",
                         "Interaction_ave", 
                         "H3K36me3_MaxScaledSignal", 
                         "chrm_acc_MaxScaledSignal")

PCA_12_SELECT_FEATURES_LABELS <- c("GC%",
                                "CpG",
                                "PhyloP-mammals",
                                "Tissue RPKM",
                                "Copies", "RNAcode",
                                "Covariance", "MFE",
                                "Accessibility",
                                "Interactions",
                                "H3K36me3",
                                "Chromatin")

PCA_20_SELECT_FEATURES <- c("GC_percentage",
                            "CpG", "TA", "GA",
                            "phyloP_max_241w", "phyloP_max_100w",
                            "GERP_91_mammals_max", "GERP_63_amniotes_max",
                            "RPKM_tissue", "RPKM_primary.cell",
                            "copy_number", "coding_potential",
                            "Max_covariance", "MFE",
                            "methylome",
                            "Interaction_ave", 
                            "H3K9ac_MaxScaledSignal", "H3K79me2_MaxScaledSignal",
                            "chrm_acc_MaxScaledSignal", "repeat_distance")

PCA_20_SELECT_FEATURES_LABELS <- c("GC%",
                                   "CpG",  "GA", "TA",
                                   "PhyloP-mammals", "PhyloP-vertebrates",
                                   "GERP-mammals", "GERP-vertebrates",
                                   "Tissue RPKM", "Primary cell RPKM", 
                                   "Copies", "RNAcode",
                                   "Covariance", "MFE",
                                   "Methylome",
                                   "Interactions",
                                   "H3K9ac", "H3K79me2",
                                   "Chromatin", "Repeat free")

#Use prot_zscores_all
prot_zscores_all <- rbind(protein_functional_z_scores,protein_negative_z_scores)
protein_data_normalized <- prot_zscores_all |> 
  dplyr::select(all_of(PCA_20_SELECT_FEATURES), Dataset) |> 
  na.omit()


# Use lncrna_zscores_all
unique(lncrna_functional_z_scores$Dataset)
lncrna_zscores_all <- rbind(lncrna_functional_z_scores,lncrna_negative_z_scores)
lncrna_data_normalized <- lncrna_zscores_all |> 
  dplyr::select(all_of(PCA_20_SELECT_FEATURES), Dataset) |> 
  na.omit()


# Use sncrna_zscores_all
sncrna_zscores_all <- rbind(sncrna_functional_z_scores,sncrna_negative_z_scores)
sncrna_data_normalized <- sncrna_zscores_all |> 
  dplyr::select(all_of(PCA_20_SELECT_FEATURES), Dataset) |> 
  na.omit()



############
# PCA Analisys on z-scores for protein coding
protein_pca_result <- prcomp(protein_data_normalized |> dplyr::select(all_of(PCA_20_SELECT_FEATURES)) |> scale(), center = TRUE, rank. = 5)
summary(protein_pca_result)
protein_pc_scores <- as.data.frame(protein_pca_result$x[,1:2])

protein_dataset_vector <- protein_data_normalized$Dataset
protein_pc_scores$`Gene type` <- factor(protein_dataset_vector,
                                      levels = c("protein-coding-exon2","protein-coding-exon3",
                                                 "protein-exon2-negative-control","protein-exon3-negative-control"),
                                      labels = c("Protein coding","Protein coding",
                                                 "Negative controls","Negative controls"))
#Steps for 2d contours
protein_pc_scores_pos <- protein_pc_scores %>%
  filter(`Gene type`=="Protein coding")
protein_pc_scores_neg <- protein_pc_scores %>%
  filter(`Gene type`=="Negative controls")
lims <- c(range(protein_pc_scores$PC1), range(protein_pc_scores$PC2))
density_func_prot <- kde2d(protein_pc_scores_pos$PC1, protein_pc_scores_pos$PC2, n = 80, lims = lims)
density_neg_prot <- kde2d(protein_pc_scores_neg$PC1, protein_pc_scores_neg$PC2, n = 80, lims = lims)

filled.contour(density_func_prot, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_prot, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")

df_density_pos_prot <- expand.grid(PC1 = density_func_prot$x, PC2 = density_func_prot$y)
df_density_pos_prot$z <- as.vector(density_func_prot$z)

df_density_neg_prot <- expand.grid(PC1 = density_neg_prot$x, PC2 = density_neg_prot$y)
df_density_neg_prot$z <- as.vector(density_neg_prot$z)

# Create a factor for sequence_type for shapes
df_density_pos_prot$`Gene type` <- c("Protein coding")
df_density_pos_prot$`Gene type` <- factor(df_density_pos_prot$`Gene type`)
df_density_neg_prot$`Gene type` <- c("Negative controls")
df_density_neg_prot$`Gene type` <- factor(df_density_neg_prot$`Gene type`)

# Plot PCA
pca_protein_plot <- ggplot(protein_pc_scores, aes(x = PC1, y = PC2, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "Negative controls"),
             aes(x = PC1, y = PC2), shape = 4, color = "#c9e3f6FF", size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_prot, aes(x = PC1, y = PC2, z = z), 
               bins= 20, color = "#72b6e7", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "Protein coding"),
             aes(x = PC1, y = PC2), shape = 15, color = "#F4A582FF", size = 4) +
  # Add semi-transparent density contours for protein coding
  geom_contour(data = df_density_pos_prot, aes(x = PC1, y = PC2, z = z), 
               bins = 20, color = "#9c3a0e", alpha = 0.6, linewidth = 1) +
  labs(subtitle="mRNA", x = "PC1 (23.27%)", y = "PC2 (12.83%)") +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "none",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-10, 5) +
  ylim(-5, 5)
pca_protein_plot
ggsave("pca_protein_20_features_scale_new_colours_new.png",path = "../results/latest1000all/pca/improved/", pca_protein_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)



# PCA Analisys on z-scores for lncRNA
lncrna_pca_result <- prcomp(lncrna_data_normalized |> dplyr::select(all_of(PCA_20_SELECT_FEATURES)) |> scale(), center = TRUE, rank. = 5)
summary(lncrna_pca_result)
lncrna_pc_scores <- as.data.frame(lncrna_pca_result$x[,1:2])
unique(lncrna_data_normalized$Dataset)
lncrna_dataset_vector <- lncrna_data_normalized$Dataset
lncrna_pc_scores$`Gene type` <- factor(lncrna_dataset_vector,
                                     levels = c("lncrna-exon1","lncrna-exon2",
                                                "lncrna-exon1-negative-control","lncrna-exon2-negative-control"),
                                     labels = c("lncRNA","lncRNA",
                                                "Negative controls","Negative controls"))
#Steps for 2d contours
lncrna_pc_scores_pos <- lncrna_pc_scores %>%
  filter(`Gene type`=="lncRNA")
lncrna_pc_scores_neg <- lncrna_pc_scores %>%
  filter(`Gene type`=="Negative controls")
lims <- c(range(lncrna_pc_scores$PC1), range(lncrna_pc_scores$PC2))
density_func_lncrna <- kde2d(lncrna_pc_scores_pos$PC1, lncrna_pc_scores_pos$PC2, n = 80, lims = lims)
density_neg_lncrna <- kde2d(lncrna_pc_scores_neg$PC1, lncrna_pc_scores_neg$PC2, n = 80, lims = lims)

filled.contour(density_func_lncrna, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_lncrna, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")

df_density_pos_lncrna <- expand.grid(PC1 = density_func_lncrna$x, PC2 = density_func_lncrna$y)
df_density_pos_lncrna$z <- as.vector(density_func_lncrna$z)

df_density_neg_lncrna <- expand.grid(PC1 = density_neg_lncrna$x, PC2 = density_neg_lncrna$y)
df_density_neg_lncrna$z <- as.vector(density_neg_lncrna$z)

# Create a factor for sequence_type for shapes
df_density_pos_lncrna$`Gene type` <- c("lncRNA")
df_density_pos_lncrna$`Gene type` <- factor(df_density_pos_lncrna$`Gene type`)
df_density_neg_lncrna$`Gene type` <- c("Negative controls")
df_density_neg_lncrna$`Gene type` <- factor(df_density_neg_lncrna$`Gene type`)

# Plot PCA
pca_lncrna_plot <- ggplot(lncrna_pc_scores, aes(x = PC1, y = PC2, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "Negative controls"),
             aes(x = PC1, y = PC2), shape = 4, color = "#53a4f5FF", size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_lncrna, aes(x = PC1, y = PC2, z = z), 
               bins= 20, color = "#0c71d6", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "lncRNA"),
             aes(x = PC1, y = PC2), shape = 19, color = "#e37b88FF", size = 4) +
  # Add semi-transparent density contours for lncRNA
  geom_contour(data = df_density_pos_lncrna, aes(x = PC1, y = PC2, z = z), 
               bins = 20, color = "#781a25", alpha = 0.6, linewidth = 1) +
  labs(subtitle = "lncRNA", x = "PC1 (22.22%)", y = "PC2 (13.57%)") +
  theme_minimal()+
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "none",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
  ) +
  xlim(-10, 5) + 
  ylim(-5, 5)
pca_lncrna_plot
ggsave("pca_lncrna_20_features_scale_new_colours_new.png",path = "../results/latest1000all/pca/improved/", pca_lncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


# PCA Analisys on z-scores for sncRNA
sncrna_pca_result <- prcomp(sncrna_data_normalized |> dplyr::select(all_of(PCA_20_SELECT_FEATURES)) |> scale(), scale. = TRUE, center = TRUE, rank. = 5)
summary(sncrna_pca_result)
sncrna_pc_scores <- as.data.frame(sncrna_pca_result$x[,1:2])

sncrna_dataset_vector <- sncrna_data_normalized$Dataset
sncrna_pc_scores$`Gene type` <- factor(sncrna_dataset_vector,
                                     levels = c("short-ncrna",
                                                "short-ncrna-negative-control"),
                                     labels = c("sncRNA",
                                                "Negative controls"))
#Steps for 2d contours
sncrna_pc_scores_pos <- sncrna_pc_scores %>%
  filter(`Gene type`=="sncRNA")
sncrna_pc_scores_neg <- sncrna_pc_scores %>%
  filter(`Gene type`=="Negative controls")
lims <- c(range(sncrna_pc_scores$PC1), range(sncrna_pc_scores$PC2))
density_func_sncrna <- kde2d(sncrna_pc_scores_pos$PC1, sncrna_pc_scores_pos$PC2, n = 80, lims = lims)
density_neg_sncrna <- kde2d(sncrna_pc_scores_neg$PC1, sncrna_pc_scores_neg$PC2, n = 80, lims = lims)

filled.contour(density_func_sncrna, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_sncrna, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")

df_density_pos_sncrna <- expand.grid(PC1 = density_func_sncrna$x, PC2 = density_func_sncrna$y)
df_density_pos_sncrna$z <- as.vector(density_func_sncrna$z)

df_density_neg_sncrna <- expand.grid(PC1 = density_neg_sncrna$x, PC2 = density_neg_sncrna$y)
df_density_neg_sncrna$z <- as.vector(density_neg_sncrna$z)

# Create a factor for sequence_type for shapes
df_density_pos_sncrna$`Gene type` <- c("sncRNA")
df_density_pos_sncrna$`Gene type` <- factor(df_density_pos_sncrna$`Gene type`)
df_density_neg_sncrna$`Gene type` <- c("Negative controls")
df_density_neg_sncrna$`Gene type` <- factor(df_density_neg_sncrna$`Gene type`)
# Plot PCA
pca_sncrna_plot <- ggplot(sncrna_pc_scores, aes(x = PC1, y = PC2, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "Negative controls"),
             aes(x = PC1, y = PC2), shape = 4, color = "#56bdfcFF", size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_sncrna, aes(x = PC1, y = PC2, z = z), 
               bins= 20, color = "#0491e8", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "sncRNA"),
             aes(x = PC1, y = PC2), shape = 17, color = "#D6604DFF", size = 4) +
  # Add semi-transparent density contours for sncRNA
  geom_contour(data = df_density_pos_sncrna, aes(x = PC1, y = PC2, z = z), 
               bins = 20, color = "#471810", alpha = 0.6, linewidth = 1) +
  labs(subtitle = "sncRNA", x = "PC1 (26.55%)", y = "PC2 (12.30%)") +
  theme_minimal()+
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "none",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-10, 5) +
  ylim(-5, 5)
pca_sncrna_plot
ggsave("pca_sncrna_20_features_scale_new_colours_new.png",path = "../results/latest1000all/pca/improved/", pca_sncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


# Join the plots in a single patchwork
pca_joined_plot <- (pca_protein_plot + pca_sncrna_plot + pca_lncrna_plot) +
  plot_annotation(title = "Principal Component Analysis",
                  tag_levels = list(c('A','B','C'))) +
  plot_layout(axis_titles = "collect", guides = "collect") &
  theme(plot.title = element_text(size = 56, hjust = 0.5, margin = ggplot2::margin(0,0,40,0)),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 34, face = "bold", hjust = 0, vjust = 0))
pca_joined_plot
ggsave("pca_20_features_scale_new_colours_new.png",path = "../results/latest1000all/pca/improved/new/", pca_joined_plot, scale = 3, width = 3840, height = 1620, units = "px", bg = "white", dpi = 600)
