# This script computes a single PCA analysis of the gene functionality features data.
# The script produces contour plots for each type of RNA: mRNA, lncRNA, and sncRNA.

# Declare required libraries
library(magrittr)
library(dplyr)
library(ggpubr)
library(MASS)
library(patchwork)
library(factoextra)

# Load data
source("scripts/config.R")
source("scripts/load_gene_functionality_zscores.R")
zscores_all <- load_gene_functionality_zscores()

# Define features to analyze
PCA_12_SELECT_FEATURES <- c("GC_percentage",
                         "CpG",
                         "phyloP_max_241w",
                         "RPKM_tissue",
                         "copy_number", "coding_potential",
                         "Max_covariance", "MFE",
                         "methylome",
                         "Interaction_ave", 
                         "H3K79me2_MaxScaledSignal", 
                         "chrm_acc_MaxScaledSignal")

PCA_12_SELECT_FEATURES_LABELS <- c("GC%",
                                "CpG",
                                "PhyloP-mammals",
                                "Tissue RPKM",
                                "Copies", "RNAcode",
                                "Covariance", "MFE",
                                "Accessibility",
                                "Interactions",
                                "H3K79me2",
                                "Chromatin")

PCA_20_SELECT_FEATURES <- c("GC_percentage",
                            "CpG", "TA", "GA",
                            "phyloP_max_241w", "phyloP_max_100w",
                            "GERP_91_mammals_max", "GERP_63_amniotes_max",
                            "RPKM_tissue", "RPKM_primary.cell",
                            "copy_number", 
                            "coding_potential",
                            "Max_covariance", "MFE",
                            "methylome",
                            "Interaction_ave", 
                            "H3K9ac_MaxScaledSignal", "H3K79me2_MaxScaledSignal", "H3K79me1_MaxScaledSignal",
                            "chrm_acc_MaxScaledSignal")#, "repeat_distance")

PCA_20_SELECT_FEATURES_LABELS <- c("GC%",
                                   "CpG",  "GA", "TA",
                                   "PhyloP-mammals", "PhyloP-vertebrates",
                                   "GERP-mammals", "GERP-vertebrates",
                                   "Tissue RPKM", "Primary cell RPKM", 
                                   "Copies", "RNAcode",
                                   "Covariance", "MFE",
                                   "Methylome",
                                   "Interactions",
                                   "H3K9ac", "H3K79me2", "H3K79me1",
                                   "Chromatin")#, "Repeat free")

# Select mRNA data
prot_zscores_all <- rbind(zscores_all %>% filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")),
                          zscores_all %>% filter(Dataset %in% c("protein-exon2-negative-control", "protein-exon3-negative-control")))
protein_data_normalized <- prot_zscores_all |> 
  dplyr::select(all_of(PCA_20_SELECT_FEATURES), Dataset) |> 
  na.omit()


# Select lncRNA data
lncrna_zscores_all <- rbind(zscores_all %>% filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2")),
                            zscores_all %>% filter(Dataset %in% c("lncrna-exon1-negative-control", "lncrna-exon2-negative-control")))
lncrna_data_normalized <- lncrna_zscores_all |> 
  dplyr::select(all_of(PCA_20_SELECT_FEATURES), Dataset) |> 
  na.omit()


# Select sncRNA data
sncrna_zscores_all <- rbind(zscores_all %>% filter(Dataset %in% c("short-ncrna")),
                            zscores_all %>% filter(Dataset %in% c("short-ncrna-negative-control")))
sncrna_data_normalized <- sncrna_zscores_all |> 
  dplyr::select(all_of(PCA_20_SELECT_FEATURES), Dataset) |> 
  na.omit()

# Select all data
allrna_data_normalized <- zscores_all |>
  dplyr::select(all_of(PCA_20_SELECT_FEATURES), Dataset) |>
  na.omit() #|>
  #mutate(across(
  #  .cols = where(is.numeric),
  #  .fns = ~as.vector(scale(.)),
  #  .names = "{.col}_scaled"
  #))

summary(allrna_data_normalized$RPKM_primary.cell)
#summary(allrna_data_normalized$RPKM_primary.cell_scaled)

# View the result
glimpse(allrna_data_normalized)


ggplot(allrna_data_normalized, aes(x = coding_potential)) +
  geom_density() + 
  coord_cartesian(xlim = c(-1, 15)) +
  stat_ecdf(geom = "step")

#########################################################
# --- PCA for all z-scores, including all RNA types --- #
allrna_pca_result <- prcomp(allrna_data_normalized |> dplyr::select(-Dataset) |> scale(), center = FALSE, rank. = 6)
summary(allrna_pca_result)
#unique(allrna_data_normalized$Dataset)
allrna_pc_scores <- as.data.frame(allrna_pca_result$x)
allrna_pc_scores$`Gene type` <- factor(allrna_data_normalized$Dataset,
                                       levels = unique(allrna_data_normalized$Dataset),
                                       labels = c("mRNA(+)","mRNA(+)",
                                                  "lncRNA(+)","lncRNA(+)",
                                                  "sncRNA(+)",
                                                  "mRNA(-)","mRNA(-)",
                                                  "lncRNA(-)","lncRNA(-)",
                                                  "sncRNA(-)"))

# Generate Scree plot, and feature variation contribution plots for combinations of PC1 vs other PCs.
# Change colnames
colnames(allrna_data_normalized) <- c(PCA_20_SELECT_FEATURES_LABELS, "Dataset")

# Compute contribution of features to variation of PCAs
allrna_pca <- princomp(allrna_data_normalized |> dplyr::select(-Dataset) |> scale())
allrna_pca_scores <- as.data.frame(allrna_pca$scores)
allrna_pca_scores$`Gene type` <- factor(allrna_data_normalized$Dataset,
                                       levels = unique(allrna_data_normalized$Dataset),
                                       labels = c("mRNA(+)","mRNA(+)",
                                                  "lncRNA(+)","lncRNA(+)",
                                                  "sncRNA(+)",
                                                  "mRNA(-)","mRNA(-)",
                                                  "lncRNA(-)","lncRNA(-)",
                                                  "sncRNA(-)"))
summary(allrna_pca)

# Inspect variable loading
allrna_pca$loadings[, 1:6]

# Scree Plot
p1 <- fviz_eig(allrna_pca, addlabels = TRUE, title = "Scree plot - all RNA")
p1_mod <- p1 + labs(x = "Principal Components") + theme(text = element_text(size = 26))
p1_mod$layers[[4]]$aes_params$size <- 8
print(p1_mod)

# Graph of the variables
fviz_pca_var(allrna_pca, col.var = "black")

# Contribution of each variable
# 1. Create the plot object
p <- fviz_cos2(allrna_pca, choice = "var", axes = 1, title = "Contribution of Features to PC 1 Variance - all RNA")

# 2. Modify the plot object's aesthetics
p_modified <- p + labs(y = "Contribution to PC1 Variance") + 
  theme(text = element_text(size = 26),
        axis.text.x = element_text(angle = 90, vjust = 0.5))

# 3. Print the modified plot
print(p_modified)


# Biplot
fviz_pca_var(allrna_pca, 
             axes = c(1, 2),
             col.var = "cos2",
             gradient.cols = c("black", "orchid", "blue"),
             repel = TRUE,
             title = "Biplot - All RNA",
             xlab = "PC1 (22.4%)",
             ylab = "PC2 (12.6%)"
) +
  theme(
    text = element_text(size = 26)
  )

fviz_pca_var(allrna_pca,
             axes = c(1, 3), # This specifies plotting PC1 vs PC3
             col.var = "cos2",
             gradient.cols = c("black", "orchid", "blue"),
             repel = TRUE,
             title = "Biplot - All RNA",
             xlab = "PC1 (22.4%)", # This remains the same
             ylab = "PC3 (7.8%)"
) +
  theme(
    text = element_text(size = 26)
  )

fviz_pca_var(allrna_pca,
             axes = c(1, 4), # This specifies plotting PC1 vs PC4
             col.var = "cos2",
             gradient.cols = c("black", "orchid", "blue"),
             repel = TRUE,
             title = "Biplot - All RNA",
             xlab = "PC1 (22.4%)", # This remains the same
             ylab = "PC4 (6.9%)"
) +
  theme(
    text = element_text(size = 26)
  )

fviz_pca_var(allrna_pca,
             axes = c(1, 5), # This specifies plotting PC1 vs PC5
             col.var = "cos2",
             gradient.cols = c("black", "orchid", "blue"),
             repel = TRUE,
             title = "Biplot - All RNA",
             xlab = "PC1 (22.4%)", # This remains the same
             ylab = "PC5 (6.1%)"
) +
  theme(
    text = element_text(size = 26)
  )

fviz_pca_var(allrna_pca,
             axes = c(1, 6), # This specifies plotting PC1 vs PC6
             col.var = "cos2",
             gradient.cols = c("black", "orchid", "blue"),
             repel = TRUE,
             title = "Biplot - All RNA",
             xlab = "PC1 (22.4%)", # This remains the same
             ylab = "PC6 (5.3%)"
) +
  theme(
    text = element_text(size = 26)
  )

# Print the loadings matrix to see the coefficients
print(allrna_pca$loadings)

#############################################
# --- PCA Analisys on z-scores for mRNA --- #
#############################################
#protein_pca_result <- prcomp(protein_data_normalized |> dplyr::select(-Dataset) |> scale(), center = TRUE, rank. = 5)
#summary(protein_pca_result)
#protein_pc_scores <- as.data.frame(allrna_pca_result$x)

#protein_dataset_vector <- protein_data_normalized$Dataset
protein_pc_scores <- allrna_pca_scores %>%
  filter(`Gene type`=="mRNA(+)" | `Gene type`=="mRNA(-)")
#Steps for 2d contours
protein_pc_scores_pos <- allrna_pca_scores %>%
  filter(`Gene type`=="mRNA(+)")
protein_pc_scores_neg <- allrna_pca_scores %>%
  filter(`Gene type`=="mRNA(-)")
lims2 <- c(range(protein_pc_scores$Comp.1), range(protein_pc_scores$Comp.2))
lims3 <- c(range(protein_pc_scores$Comp.1), range(protein_pc_scores$Comp.3))
lims4 <- c(range(protein_pc_scores$Comp.1), range(protein_pc_scores$Comp.4))
lims5 <- c(range(protein_pc_scores$Comp.1), range(protein_pc_scores$Comp.5))
lims6 <- c(range(protein_pc_scores$Comp.1), range(protein_pc_scores$Comp.6))
density_func_prot2 <- kde2d(protein_pc_scores_pos$Comp.1, protein_pc_scores_pos$Comp.2, n = 80, lims = lims2)
density_neg_prot2 <- kde2d(protein_pc_scores_neg$Comp.1, protein_pc_scores_neg$Comp.2, n = 80, lims = lims2)
density_func_prot3 <- kde2d(protein_pc_scores_pos$Comp.1, protein_pc_scores_pos$Comp.3, n = 80, lims = lims3)
density_neg_prot3 <- kde2d(protein_pc_scores_neg$Comp.1, protein_pc_scores_neg$Comp.3, n = 80, lims = lims3)
density_func_prot4 <- kde2d(protein_pc_scores_pos$Comp.1, protein_pc_scores_pos$Comp.4, n = 80, lims = lims4)
density_neg_prot4 <- kde2d(protein_pc_scores_neg$Comp.1, protein_pc_scores_neg$Comp.4, n = 80, lims = lims4)
density_func_prot5 <- kde2d(protein_pc_scores_pos$Comp.1, protein_pc_scores_pos$Comp.5, n = 80, lims = lims5)
density_neg_prot5 <- kde2d(protein_pc_scores_neg$Comp.1, protein_pc_scores_neg$Comp.5, n = 80, lims = lims5)
density_func_prot6 <- kde2d(protein_pc_scores_pos$Comp.1, protein_pc_scores_pos$Comp.6, n = 80, lims = lims6)
density_neg_prot6 <- kde2d(protein_pc_scores_neg$Comp.1, protein_pc_scores_neg$Comp.6, n = 80, lims = lims6)

filled.contour(density_func_prot2, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_prot2, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_prot3, xlab = "PC1", ylab = "PC3", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_prot3, xlab = "PC1", ylab = "PC3", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_prot4, xlab = "PC1", ylab = "PC4", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_prot4, xlab = "PC1", ylab = "PC4", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_prot5, xlab = "PC1", ylab = "PC5", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_prot5, xlab = "PC1", ylab = "PC5", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_prot6, xlab = "PC1", ylab = "PC6", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_prot6, xlab = "PC1", ylab = "PC6", main = "Filled Contour Plot of PCA Density")

df_density_pos_prot2 <- expand.grid(PC1 = density_func_prot2$x, PC2 = density_func_prot2$y)
df_density_pos_prot2$z <- as.vector(density_func_prot2$z)
df_density_neg_prot2 <- expand.grid(PC1 = density_neg_prot2$x, PC2 = density_neg_prot2$y)
df_density_neg_prot2$z <- as.vector(density_neg_prot2$z)
df_density_pos_prot3 <- expand.grid(PC1 = density_func_prot3$x, PC3 = density_func_prot3$y)
df_density_pos_prot3$z <- as.vector(density_func_prot3$z)
df_density_neg_prot3 <- expand.grid(PC1 = density_neg_prot3$x, PC3 = density_neg_prot3$y)
df_density_neg_prot3$z <- as.vector(density_neg_prot3$z)
df_density_pos_prot4 <- expand.grid(PC1 = density_func_prot4$x, PC4 = density_func_prot4$y)
df_density_pos_prot4$z <- as.vector(density_func_prot4$z)
df_density_neg_prot4 <- expand.grid(PC1 = density_neg_prot4$x, PC4 = density_neg_prot4$y)
df_density_neg_prot4$z <- as.vector(density_neg_prot4$z)
df_density_pos_prot5 <- expand.grid(PC1 = density_func_prot5$x, PC5 = density_func_prot5$y)
df_density_pos_prot5$z <- as.vector(density_func_prot5$z)
df_density_neg_prot5 <- expand.grid(PC1 = density_neg_prot5$x, PC5 = density_neg_prot5$y)
df_density_neg_prot5$z <- as.vector(density_neg_prot5$z)
df_density_pos_prot6 <- expand.grid(PC1 = density_func_prot6$x, PC6 = density_func_prot6$y)
df_density_pos_prot6$z <- as.vector(density_func_prot6$z)
df_density_neg_prot6 <- expand.grid(PC1 = density_neg_prot6$x, PC6 = density_neg_prot6$y)
df_density_neg_prot6$z <- as.vector(density_neg_prot6$z)

# Create a factor for sequence_type for shapes
df_density_pos_prot2$`Gene type` <- c("mRNA(+)")
df_density_pos_prot2$`Gene type` <- factor(df_density_pos_prot2$`Gene type`)
df_density_neg_prot2$`Gene type` <- c("mRNA(-)")
df_density_neg_prot2$`Gene type` <- factor(df_density_neg_prot2$`Gene type`)

df_density_pos_prot3$`Gene type` <- c("mRNA(+)")
df_density_pos_prot3$`Gene type` <- factor(df_density_pos_prot3$`Gene type`)
df_density_neg_prot3$`Gene type` <- c("mRNA(-)")
df_density_neg_prot3$`Gene type` <- factor(df_density_neg_prot3$`Gene type`)

df_density_pos_prot4$`Gene type` <- c("mRNA(+)")
df_density_pos_prot4$`Gene type` <- factor(df_density_pos_prot4$`Gene type`)
df_density_neg_prot4$`Gene type` <- c("mRNA(-)")
df_density_neg_prot4$`Gene type` <- factor(df_density_neg_prot4$`Gene type`)

df_density_pos_prot5$`Gene type` <- c("mRNA(+)")
df_density_pos_prot5$`Gene type` <- factor(df_density_pos_prot5$`Gene type`)
df_density_neg_prot5$`Gene type` <- c("mRNA(-)")
df_density_neg_prot5$`Gene type` <- factor(df_density_neg_prot5$`Gene type`)

df_density_pos_prot6$`Gene type` <- c("mRNA(+)")
df_density_pos_prot6$`Gene type` <- factor(df_density_pos_prot6$`Gene type`)
df_density_neg_prot6$`Gene type` <- c("mRNA(-)")
df_density_neg_prot6$`Gene type` <- factor(df_density_neg_prot6$`Gene type`)


# Define the colors and shapes you want to use
my_colors <- c("mRNA(-)" = "#c9e3f6FF", "mRNA(+)" = "#F4A582FF")
my_shapes <- c("mRNA(-)" = 4, "mRNA(+)" = 15)

# Plot PC distributions
pc1_dist_plot <- ggplot(protein_pc_scores, aes(x = Comp.1, color = `Gene type`)) +
  #geom_density(size = 1.5) + # Using a single geom_density is more efficient
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("mRNA(-)" = "#72b6e7", "mRNA(+)" = "#F4A582FF")) +
  labs(subtitle="mRNA", x = "PC1 (22.4%)", y = "ECDF") + 
  coord_cartesian(c(-5, 10)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc1_dist_plot

pc2_dist_plot <- ggplot(protein_pc_scores, aes(x = Comp.2, color = `Gene type`)) +
  #geom_density(size = 1.5) + # Using a single geom_density is more efficient
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("mRNA(-)" = "#72b6e7", "mRNA(+)" = "#F4A582FF")) +
  labs(subtitle="mRNA", x = "PC2 (12.6%)", y = "ECDF") + 
  coord_cartesian(c(-5, 5)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc2_dist_plot

pc3_dist_plot <- ggplot(protein_pc_scores, aes(x = Comp.3, color = `Gene type`)) +
  #geom_density(size = 1.5) + # Using a single geom_density is more efficient
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("mRNA(-)" = "#72b6e7", "mRNA(+)" = "#F4A582FF")) +
  labs(subtitle="mRNA", x = "PC3 (7.8%)", y = "ECDF") + 
  coord_cartesian(c(-5, 10)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc3_dist_plot

pc4_dist_plot <- ggplot(protein_pc_scores, aes(x = Comp.4, color = `Gene type`)) +
  #geom_density(size = 1.5) + # Using a single geom_density is more efficient
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("mRNA(-)" = "#72b6e7", "mRNA(+)" = "#F4A582FF")) +
  labs(subtitle="mRNA", x = "PC4 (6.9%)", y = "ECDF") + 
  coord_cartesian(c(-5, 5)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc4_dist_plot

pc5_dist_plot <- ggplot(protein_pc_scores, aes(x = Comp.5, color = `Gene type`)) +
  #geom_density(size = 1.5) + # Using a single geom_density is more efficient
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("mRNA(-)" = "#72b6e7", "mRNA(+)" = "#F4A582FF")) +
  labs(subtitle="mRNA", x = "PC5 (6.1%)", y = "ECDF") + 
  coord_cartesian(c(-5, 5)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc5_dist_plot

pc6_dist_plot <- ggplot(protein_pc_scores, aes(x = Comp.6, color = `Gene type`)) +
  #geom_density(size = 1.5) + # Using a single geom_density is more efficient
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("mRNA(-)" = "#72b6e7", "mRNA(+)" = "#F4A582FF")) +
  labs(subtitle="mRNA", x = "PC6 (5.3%)", y = "ECDF") + 
  coord_cartesian(c(-5, 10)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc6_dist_plot


# Plot PCA
pca_protein_plot <- ggplot(protein_pc_scores, aes(x = Comp.1, y = Comp.2, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "mRNA(-)"),
             aes(x = Comp.1, y = Comp.2), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_prot2, aes(x = PC1, y = PC2, z = z), 
               bins= 20, color = "#72b6e7", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "mRNA(+)"),
             aes(x = Comp.1, y = Comp.2), shape = 15, size = 4) +
  # Add semi-transparent density contours for protein coding
  geom_contour(data = df_density_pos_prot2, aes(x = PC1, y = PC2, z = z), 
               bins = 20, color = "#9c3a0e", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle="mRNA", x = "PC1 (22.4%)", y = "PC2 (12.6%)") +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-5, 10) +
  ylim(-5, 5)
pca_protein_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_protein_20_features_loadings","_pc2.png"), pca_protein_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

pca_protein_plot <- ggplot(protein_pc_scores, aes(x = Comp.1, y = Comp.3, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "mRNA(-)"),
             aes(x = Comp.1, y = Comp.3), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_prot3, aes(x = PC1, y = PC3, z = z), 
               bins= 20, color = "#72b6e7", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "mRNA(+)"),
             aes(x = Comp.1, y = Comp.3), shape = 15, size = 4) +
  # Add semi-transparent density contours for protein coding
  geom_contour(data = df_density_pos_prot3, aes(x = PC1, y = PC3, z = z), 
               bins = 20, color = "#9c3a0e", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle="mRNA", x = "PC1 (22.4%)", y = "PC3 (7.8%)") +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-5, 10) +
  ylim(-5, 5)
pca_protein_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_protein_20_features_loadings","_pc3.png"), pca_protein_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

pca_protein_plot <- ggplot(protein_pc_scores, aes(x = Comp.1, y = Comp.4, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "mRNA(-)"),
             aes(x = Comp.1, y = Comp.4), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_prot4, aes(x = PC1, y = PC4, z = z), 
               bins= 20, color = "#72b6e7", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "mRNA(+)"),
             aes(x = Comp.1, y = Comp.4), shape = 15, size = 4) +
  # Add semi-transparent density contours for protein coding
  geom_contour(data = df_density_pos_prot4, aes(x = PC1, y = PC4, z = z), 
               bins = 20, color = "#9c3a0e", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle="mRNA", x = "PC1 (22.4%)", y = "PC4 (6.9%)") +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-5, 10) +
  ylim(-5, 5)
pca_protein_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_protein_20_features_loadings","_pc4.png"), pca_protein_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

pca_protein_plot <- ggplot(protein_pc_scores, aes(x = Comp.1, y = Comp.5, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "mRNA(-)"),
             aes(x = Comp.1, y = Comp.5), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_prot5, aes(x = PC1, y = PC5, z = z), 
               bins= 20, color = "#72b6e7", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "mRNA(+)"),
             aes(x = Comp.1, y = Comp.5), shape = 15, size = 4) +
  # Add semi-transparent density contours for protein coding
  geom_contour(data = df_density_pos_prot5, aes(x = PC1, y = PC5, z = z), 
               bins = 20, color = "#9c3a0e", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle="mRNA", x = "PC1 (22.4%)", y = "PC5 (6.1%)") +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-5, 10) +
  ylim(-5, 5)
pca_protein_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_protein_20_features_loadings","_pc5.png"), pca_protein_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

pca_protein_plot <- ggplot(protein_pc_scores, aes(x = Comp.1, y = Comp.6, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "mRNA(-)"),
             aes(x = Comp.1, y = Comp.6), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_prot6, aes(x = PC1, y = PC6, z = z), 
               bins= 20, color = "#72b6e7", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(protein_pc_scores, `Gene type` == "mRNA(+)"),
             aes(x = Comp.1, y = Comp.6), shape = 15, size = 4) +
  # Add semi-transparent density contours for protein coding
  geom_contour(data = df_density_pos_prot6, aes(x = PC1, y = PC6, z = z), 
               bins = 20, color = "#9c3a0e", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle="mRNA", x = "PC1 (22.4%)", y = "PC6 (5.3%)") +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-5, 10) +
  ylim(-5, 5)
pca_protein_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_protein_20_features_loadings","_pc6.png"), pca_protein_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


# Change colnames
#colnames(protein_data_normalized) <- c(PCA_20_SELECT_FEATURES_LABELS, "Dataset")
# Compute contribution of features to variation of PCAs
#mrna_pca <- princomp(protein_data_normalized |> dplyr::select(-Dataset) |> scale())
#summary(mrna_pca)

# Inspect variable loadings
#mrna_pca$loadings[, 1:3]

# Scree Plot
#p1 <- fviz_eig(mrna_pca, addlabels = TRUE, title = "Scree plot - mRNA")
#p1_mod <- p1 + labs(x = "Principal Components") + theme(text = element_text(size = 26))
#p1_mod$layers[[4]]$aes_params$size <- 8
#print(p1_mod)

# Graph of the variables
#fviz_pca_var(mrna_pca, col.var = "black")

# Contribution of each variable
#fviz_cos2(mrna_pca, choice = "var", axes = 1, title = "Contribution of Features to PC 1 Variance - mRNA")

# 1. Create the plot object, but don't print it yet
#p <- fviz_cos2(mrna_pca, choice = "var", axes = 5, title = "Contribution of Features to PC 5 Variance - mRNA")

# 2. Modify the plot object's aesthetics
# The labels are text, so we target the 'text' aesthetic.
# We set a new, larger size for all text elements.
#p_modified <- p + labs(y = "Contribution to PC5 Variance") + theme(text = element_text(size = 26),
#                        axis.text.x = element_text(angle = 90, vjust = 0.5)) # Change 16 to your desired size

# 3. Print the modified plot
#print(p_modified)

# Biplot
#fviz_pca_var(mrna_pca, col.var = "cos2",
#             gradient.cols = c("black", "orchid", "blue"),
#             repel = TRUE,
#             title = "Biplot - mRNA",
#             xlab = "PC1 (23.3%)",
#             ylab = "PC2 (12.8%)"
#             ) +
#  theme(
#    text = element_text(size = 26)
#  )


###############################################
# --- PCA Analisys on z-scores for lncRNA --- #
###############################################
#lncrna_pca_result <- prcomp(lncrna_data_normalized |> dplyr::select(-Dataset) |> scale(), center = TRUE, rank. = 5)
#summary(lncrna_pca_result)

lncrna_pc_scores <- allrna_pca_scores %>%
  filter(`Gene type`=="lncRNA(+)" | `Gene type`=="lncRNA(-)")

#Steps for 2d contours
lncrna_pc_scores_pos <- allrna_pca_scores %>%
  filter(`Gene type`=="lncRNA(+)")
lncrna_pc_scores_neg <- allrna_pca_scores %>%
  filter(`Gene type`=="lncRNA(-)")
lims2 <- c(range(lncrna_pc_scores$Comp.1), range(lncrna_pc_scores$Comp.2))
lims3 <- c(range(lncrna_pc_scores$Comp.1), range(lncrna_pc_scores$Comp.3))
lims4 <- c(range(lncrna_pc_scores$Comp.1), range(lncrna_pc_scores$Comp.4))
lims5 <- c(range(lncrna_pc_scores$Comp.1), range(lncrna_pc_scores$Comp.5))
lims6 <- c(range(lncrna_pc_scores$Comp.1), range(lncrna_pc_scores$Comp.6))
density_func_lncrna2 <- kde2d(lncrna_pc_scores_pos$Comp.1, lncrna_pc_scores_pos$Comp.2, n = 80, lims = lims2)
density_neg_lncrna2 <- kde2d(lncrna_pc_scores_neg$Comp.1, lncrna_pc_scores_neg$Comp.2, n = 80, lims = lims2)
density_func_lncrna3 <- kde2d(lncrna_pc_scores_pos$Comp.1, lncrna_pc_scores_pos$Comp.3, n = 80, lims = lims3)
density_neg_lncrna3 <- kde2d(lncrna_pc_scores_neg$Comp.1, lncrna_pc_scores_neg$Comp.3, n = 80, lims = lims3)
density_func_lncrna4 <- kde2d(lncrna_pc_scores_pos$Comp.1, lncrna_pc_scores_pos$Comp.4, n = 80, lims = lims4)
density_neg_lncrna4 <- kde2d(lncrna_pc_scores_neg$Comp.1, lncrna_pc_scores_neg$Comp.4, n = 80, lims = lims4)
density_func_lncrna5 <- kde2d(lncrna_pc_scores_pos$Comp.1, lncrna_pc_scores_pos$Comp.5, n = 80, lims = lims5)
density_neg_lncrna5 <- kde2d(lncrna_pc_scores_neg$Comp.1, lncrna_pc_scores_neg$Comp.5, n = 80, lims = lims5)
density_func_lncrna6 <- kde2d(lncrna_pc_scores_pos$Comp.1, lncrna_pc_scores_pos$Comp.6, n = 80, lims = lims6)
density_neg_lncrna6 <- kde2d(lncrna_pc_scores_neg$Comp.1, lncrna_pc_scores_neg$Comp.6, n = 80, lims = lims6)

filled.contour(density_func_lncrna2, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_lncrna2, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_lncrna3, xlab = "PC1", ylab = "PC3", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_lncrna3, xlab = "PC1", ylab = "PC3", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_lncrna4, xlab = "PC1", ylab = "PC4", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_lncrna4, xlab = "PC1", ylab = "PC4", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_lncrna5, xlab = "PC1", ylab = "PC5", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_lncrna5, xlab = "PC1", ylab = "PC5", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_lncrna6, xlab = "PC1", ylab = "PC6", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_lncrna6, xlab = "PC1", ylab = "PC6", main = "Filled Contour Plot of PCA Density")

df_density_pos_lncrna2 <- expand.grid(PC1 = density_func_lncrna2$x, PC2 = density_func_lncrna2$y)
df_density_pos_lncrna2$z <- as.vector(density_func_lncrna2$z)
df_density_pos_lncrna3 <- expand.grid(PC1 = density_func_lncrna3$x, PC3 = density_func_lncrna3$y)
df_density_pos_lncrna3$z <- as.vector(density_func_lncrna3$z)
df_density_pos_lncrna4 <- expand.grid(PC1 = density_func_lncrna4$x, PC4 = density_func_lncrna4$y)
df_density_pos_lncrna4$z <- as.vector(density_func_lncrna4$z)
df_density_pos_lncrna5 <- expand.grid(PC1 = density_func_lncrna5$x, PC5 = density_func_lncrna5$y)
df_density_pos_lncrna5$z <- as.vector(density_func_lncrna5$z)
df_density_pos_lncrna6 <- expand.grid(PC1 = density_func_lncrna6$x, PC6 = density_func_lncrna6$y)
df_density_pos_lncrna6$z <- as.vector(density_func_lncrna6$z)

df_density_neg_lncrna2 <- expand.grid(PC1 = density_neg_lncrna2$x, PC2 = density_neg_lncrna2$y)
df_density_neg_lncrna2$z <- as.vector(density_neg_lncrna2$z)
df_density_neg_lncrna3 <- expand.grid(PC1 = density_neg_lncrna3$x, PC3 = density_neg_lncrna3$y)
df_density_neg_lncrna3$z <- as.vector(density_neg_lncrna3$z)
df_density_neg_lncrna4 <- expand.grid(PC1 = density_neg_lncrna4$x, PC4 = density_neg_lncrna4$y)
df_density_neg_lncrna4$z <- as.vector(density_neg_lncrna4$z)
df_density_neg_lncrna5 <- expand.grid(PC1 = density_neg_lncrna5$x, PC5 = density_neg_lncrna5$y)
df_density_neg_lncrna5$z <- as.vector(density_neg_lncrna5$z)
df_density_neg_lncrna6 <- expand.grid(PC1 = density_neg_lncrna6$x, PC6 = density_neg_lncrna6$y)
df_density_neg_lncrna6$z <- as.vector(density_neg_lncrna6$z)

# Create a factor for sequence_type for shapes
df_density_pos_lncrna2$`Gene type` <- c("lncRNA(+)")
df_density_pos_lncrna2$`Gene type` <- factor(df_density_pos_lncrna2$`Gene type`)
df_density_neg_lncrna2$`Gene type` <- c("lncRNA(-)")
df_density_neg_lncrna2$`Gene type` <- factor(df_density_neg_lncrna2$`Gene type`)

df_density_pos_lncrna3$`Gene type` <- c("lncRNA(+)")
df_density_pos_lncrna3$`Gene type` <- factor(df_density_pos_lncrna3$`Gene type`)
df_density_neg_lncrna3$`Gene type` <- c("lncRNA(-)")
df_density_neg_lncrna3$`Gene type` <- factor(df_density_neg_lncrna3$`Gene type`)

df_density_pos_lncrna4$`Gene type` <- c("lncRNA(+)")
df_density_pos_lncrna4$`Gene type` <- factor(df_density_pos_lncrna4$`Gene type`)
df_density_neg_lncrna4$`Gene type` <- c("lncRNA(-)")
df_density_neg_lncrna4$`Gene type` <- factor(df_density_neg_lncrna4$`Gene type`)

df_density_pos_lncrna5$`Gene type` <- c("lncRNA(+)")
df_density_pos_lncrna5$`Gene type` <- factor(df_density_pos_lncrna5$`Gene type`)
df_density_neg_lncrna5$`Gene type` <- c("lncRNA(-)")
df_density_neg_lncrna5$`Gene type` <- factor(df_density_neg_lncrna5$`Gene type`)

df_density_pos_lncrna6$`Gene type` <- c("lncRNA(+)")
df_density_pos_lncrna6$`Gene type` <- factor(df_density_pos_lncrna6$`Gene type`)
df_density_neg_lncrna6$`Gene type` <- c("lncRNA(-)")
df_density_neg_lncrna6$`Gene type` <- factor(df_density_neg_lncrna6$`Gene type`)

# Define the colors and shapes you want to use
my_colors <- c("lncRNA(-)" = "#53a4f5FF", "lncRNA(+)" = "#e37b88FF")
my_shapes <- c("lncRNA(-)" = 4, "lncRNA(+)" = 19)

# Plot PC distributions
pc1_dist_plot <- ggplot(lncrna_pc_scores, aes(x = Comp.1, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("lncRNA(-)" = "#53a4f5FF", "lncRNA(+)" = "#e37b88FF")) +
  labs(subtitle="lncRNA", x = "PC1 (22.4%)", y = "ECDF") + 
  coord_cartesian(c(-5, 10)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc1_dist_plot

pc2_dist_plot <- ggplot(lncrna_pc_scores, aes(x = Comp.2, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("lncRNA(-)" = "#53a4f5FF", "lncRNA(+)" = "#e37b88FF")) +
  labs(subtitle="lncRNA", x = "PC2 (12.6%)", y = "ECDF") + 
  coord_cartesian(c(-5, 5)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc2_dist_plot

pc3_dist_plot <- ggplot(lncrna_pc_scores, aes(x = Comp.1, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("lncRNA(-)" = "#53a4f5FF", "lncRNA(+)" = "#e37b88FF")) +
  labs(subtitle="lncRNA", x = "PC3 (7.8%)", y = "ECDF") + 
  coord_cartesian(c(-5, 10)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc3_dist_plot

pc4_dist_plot <- ggplot(lncrna_pc_scores, aes(x = Comp.4, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("lncRNA(-)" = "#53a4f5FF", "lncRNA(+)" = "#e37b88FF")) +
  labs(subtitle="lncRNA", x = "PC4 (6.9%)", y = "ECDF") + 
  coord_cartesian(c(-5, 5)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc4_dist_plot

pc5_dist_plot <- ggplot(lncrna_pc_scores, aes(x = Comp.5, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("lncRNA(-)" = "#53a4f5FF", "lncRNA(+)" = "#e37b88FF")) +
  labs(subtitle="lncRNA", x = "PC5 (6.1%)", y = "ECDF") + 
  coord_cartesian(c(-5, 5)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc5_dist_plot

pc6_dist_plot <- ggplot(lncrna_pc_scores, aes(x = Comp.6, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("lncRNA(-)" = "#53a4f5FF", "lncRNA(+)" = "#e37b88FF")) +
  labs(subtitle="lncRNA", x = "PC6 (5.3%)", y = "ECDF") + 
  coord_cartesian(c(-5, 10)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc6_dist_plot

# Plot PCA
pca_lncrna_plot <- ggplot(lncrna_pc_scores, aes(x = Comp.1, y = Comp.2, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "lncRNA(-)"),
             aes(x = Comp.1, y = Comp.2), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_lncrna2, aes(x = PC1, y = PC2, z = z), 
               bins= 20, color = "#0c71d6", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "lncRNA(+)"),
             aes(x = Comp.1, y = Comp.2), shape = 19, size = 4) +
  # Add semi-transparent density contours for lncRNA
  geom_contour(data = df_density_pos_lncrna2, aes(x = PC1, y = PC2, z = z), 
               bins = 20, color = "#781a25", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "lncRNA", x = "PC1 (22.4%)", y = "PC2 (12.6%)") +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
  ) +
  xlim(-5, 10) + 
  ylim(-5, 5)
pca_lncrna_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_lncrna_20_features_loadings","_pc2.png"), pca_lncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# Plot PCA
pca_lncrna_plot <- ggplot(lncrna_pc_scores, aes(x = Comp.1, y = Comp.3, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "lncRNA(-)"),
             aes(x = Comp.1, y = Comp.3), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_lncrna3, aes(x = PC1, y = PC3, z = z), 
               bins= 20, color = "#0c71d6", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "lncRNA(+)"),
             aes(x = Comp.1, y = Comp.3), shape = 19, size = 4) +
  # Add semi-transparent density contours for lncRNA
  geom_contour(data = df_density_pos_lncrna3, aes(x = PC1, y = PC3, z = z), 
               bins = 20, color = "#781a25", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "lncRNA", x = "PC1 (22.4%)", y = "PC3 (7.8%)") +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
  ) +
  xlim(-5, 10) + 
  ylim(-5, 5)
pca_lncrna_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_lncrna_20_features_loadings","_pc3.png"), pca_lncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# Plot PCA
pca_lncrna_plot <- ggplot(lncrna_pc_scores, aes(x = Comp.1, y = Comp.4, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "lncRNA(-)"),
             aes(x = Comp.1, y = Comp.4), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_lncrna4, aes(x = PC1, y = PC4, z = z), 
               bins= 20, color = "#0c71d6", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "lncRNA(+)"),
             aes(x = Comp.1, y = Comp.4), shape = 19, size = 4) +
  # Add semi-transparent density contours for lncRNA
  geom_contour(data = df_density_pos_lncrna4, aes(x = PC1, y = PC4, z = z), 
               bins = 20, color = "#781a25", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "lncRNA", x = "PC1 (22.4%)", y = "PC4 (6.9%)") +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
  ) +
  xlim(-5, 10) + 
  ylim(-5, 5)
pca_lncrna_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_lncrna_20_features_loadings","_pc4.png"), pca_lncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# Plot PCA
pca_lncrna_plot <- ggplot(lncrna_pc_scores, aes(x = Comp.1, y = Comp.5, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "lncRNA(-)"),
             aes(x = Comp.1, y = Comp.5), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_lncrna5, aes(x = PC1, y = PC5, z = z), 
               bins= 20, color = "#0c71d6", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "lncRNA(+)"),
             aes(x = Comp.1, y = Comp.5), shape = 19, size = 4) +
  # Add semi-transparent density contours for lncRNA
  geom_contour(data = df_density_pos_lncrna5, aes(x = PC1, y = PC5, z = z), 
               bins = 20, color = "#781a25", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "lncRNA", x = "PC1 (22.4%)", y = "PC5 (6.1%)") +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
  ) +
  xlim(-5, 10) + 
  ylim(-5, 5)
pca_lncrna_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_lncrna_20_features_loadings","_pc5.png"), pca_lncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# Plot PCA
pca_lncrna_plot <- ggplot(lncrna_pc_scores, aes(x = Comp.1, y = Comp.6, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "lncRNA(-)"),
             aes(x = Comp.1, y = Comp.6), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_lncrna6, aes(x = PC1, y = PC6, z = z), 
               bins= 20, color = "#0c71d6", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(lncrna_pc_scores, `Gene type` == "lncRNA(+)"),
             aes(x = Comp.1, y = Comp.6), shape = 19, size = 4) +
  # Add semi-transparent density contours for lncRNA
  geom_contour(data = df_density_pos_lncrna6, aes(x = PC1, y = PC6, z = z), 
               bins = 20, color = "#781a25", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "lncRNA", x = "PC1 (22.4%)", y = "PC6 (5.3%)") +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
  ) +
  xlim(-5, 10) + 
  ylim(-5, 5)
pca_lncrna_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_lncrna_20_features_loadings","_pc6.png"), pca_lncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# Change colnames
#colnames(lncrna_data_normalized) <- c(PCA_20_SELECT_FEATURES_LABELS, "Dataset")
# Compute contribution of features to variation of PCAs
#lncrna_pca <- princomp(lncrna_data_normalized |> dplyr::select(-Dataset) |> scale())
#summary(lncrna_pca)

# Inspect variable loadings
#lncrna_pca$loadings[, 1:3]

# Scree Plot
#p1 <- fviz_eig(lncrna_pca, addlabels = TRUE, title = "Scree plot - lncRNA")
#p1_mod <- p1 + labs(x = "Principal Components") + theme(text = element_text(size = 26))
#1_mod$layers[[4]]$aes_params$size <- 8
#print(p1_mod)

# Graph of the variables
#fviz_pca_var(lncrna_pca, col.var = "black")

# Contribution of each variable
#fviz_cos2(lncrna_pca, choice = "var", axes = 1, title = "Contribution of Features to PC 1 Variance - lncRNA")

# 1. Create the plot object, but don't print it yet
#p <- fviz_cos2(lncrna_pca, choice = "var", axes = 1, title = "Contribution of Features to PC 1 Variance - lncRNA")

# 2. Modify the plot object's aesthetics
# The labels are text, so we target the 'text' aesthetic.
# We set a new, larger size for all text elements.
#p_modified <- p + labs(y = "Contribution to PC1 Variance") + theme(text = element_text(size = 26),
#                        axis.text.x = element_text(angle = 90, vjust = 0.5)) # Change 16 to your desired size

# 3. Print the modified plot
#print(p_modified)

# Biplot
#fviz_pca_var(lncrna_pca, col.var = "cos2",
#             gradient.cols = c("black", "orchid", "blue"),
#             repel = TRUE,
#             title = "Biplot - lncRNA",
#             xlab = "PC1 (22.2%)",
#             ylab = "PC2 (13.6%)"
#             ) +
#  theme(
#    text = element_text(size = 26)
#  )


###############################################
# --- PCA Analisys on z-scores for sncRNA --- #
###############################################
#sncrna_pca_result <- prcomp(sncrna_data_normalized |> dplyr::select(-Dataset) |> scale(), scale. = TRUE, center = TRUE, rank. = 5)
#summary(sncrna_pca_result)
#sncrna_pc_scores <- as.data.frame(sncrna_pca_result$x)

#sncrna_dataset_vector <- sncrna_data_normalized$Dataset
#sncrna_pc_scores$`Gene type` <- factor(sncrna_dataset_vector,
#                                     levels = c("short-ncrna",
#                                                "short-ncrna-negative-control"),
#                                     labels = c("sncRNA",
#                                                "Negative controls"))

sncrna_pc_scores <- allrna_pca_scores %>%
  filter(`Gene type`=="sncRNA(+)" | `Gene type`=="sncRNA(-)")

#Steps for 2d contours
sncrna_pc_scores_pos <- allrna_pca_scores %>%
  filter(`Gene type`=="sncRNA(+)")
sncrna_pc_scores_neg <- allrna_pca_scores %>%
  filter(`Gene type`=="sncRNA(-)")
lims2 <- c(range(sncrna_pc_scores$Comp.1), range(sncrna_pc_scores$Comp.2))
lims3 <- c(range(sncrna_pc_scores$Comp.1), range(sncrna_pc_scores$Comp.3))
lims4 <- c(range(sncrna_pc_scores$Comp.1), range(sncrna_pc_scores$Comp.4))
lims5 <- c(range(sncrna_pc_scores$Comp.1), range(sncrna_pc_scores$Comp.5))
lims6 <- c(range(sncrna_pc_scores$Comp.1), range(sncrna_pc_scores$Comp.6))
density_func_sncrna2 <- kde2d(sncrna_pc_scores_pos$Comp.1, sncrna_pc_scores_pos$Comp.2, n = 80, lims = lims2)
density_neg_sncrna2 <- kde2d(sncrna_pc_scores_neg$Comp.1, sncrna_pc_scores_neg$Comp.2, n = 80, lims = lims2)
density_func_sncrna3 <- kde2d(sncrna_pc_scores_pos$Comp.1, sncrna_pc_scores_pos$Comp.3, n = 80, lims = lims3)
density_neg_sncrna3 <- kde2d(sncrna_pc_scores_neg$Comp.1, sncrna_pc_scores_neg$Comp.3, n = 80, lims = lims3)
density_func_sncrna4 <- kde2d(sncrna_pc_scores_pos$Comp.1, sncrna_pc_scores_pos$Comp.4, n = 80, lims = lims4)
density_neg_sncrna4 <- kde2d(sncrna_pc_scores_neg$Comp.1, sncrna_pc_scores_neg$Comp.4, n = 80, lims = lims4)
density_func_sncrna5 <- kde2d(sncrna_pc_scores_pos$Comp.1, sncrna_pc_scores_pos$Comp.5, n = 80, lims = lims5)
density_neg_sncrna5 <- kde2d(sncrna_pc_scores_neg$Comp.1, sncrna_pc_scores_neg$Comp.5, n = 80, lims = lims5)
density_func_sncrna6 <- kde2d(sncrna_pc_scores_pos$Comp.1, sncrna_pc_scores_pos$Comp.6, n = 80, lims = lims6)
density_neg_sncrna6 <- kde2d(sncrna_pc_scores_neg$Comp.1, sncrna_pc_scores_neg$Comp.6, n = 80, lims = lims6)

filled.contour(density_func_sncrna2, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_sncrna2, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_sncrna3, xlab = "PC1", ylab = "PC3", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_sncrna3, xlab = "PC1", ylab = "PC3", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_sncrna4, xlab = "PC1", ylab = "PC4", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_sncrna4, xlab = "PC1", ylab = "PC4", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_sncrna5, xlab = "PC1", ylab = "PC5", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_sncrna5, xlab = "PC1", ylab = "PC5", main = "Filled Contour Plot of PCA Density")
filled.contour(density_func_sncrna6, xlab = "PC1", ylab = "PC6", main = "Filled Contour Plot of PCA Density")
filled.contour(density_neg_sncrna6, xlab = "PC1", ylab = "PC6", main = "Filled Contour Plot of PCA Density")

df_density_pos_sncrna2 <- expand.grid(PC1 = density_func_sncrna2$x, PC2 = density_func_sncrna2$y)
df_density_pos_sncrna2$z <- as.vector(density_func_sncrna2$z)
df_density_neg_sncrna2 <- expand.grid(PC1 = density_neg_sncrna2$x, PC2 = density_neg_sncrna2$y)
df_density_neg_sncrna2$z <- as.vector(density_neg_sncrna2$z)
df_density_pos_sncrna3 <- expand.grid(PC1 = density_func_sncrna3$x, PC3 = density_func_sncrna3$y)
df_density_pos_sncrna3$z <- as.vector(density_func_sncrna3$z)
df_density_neg_sncrna3 <- expand.grid(PC1 = density_neg_sncrna3$x, PC3 = density_neg_sncrna3$y)
df_density_neg_sncrna3$z <- as.vector(density_neg_sncrna3$z)
df_density_pos_sncrna4 <- expand.grid(PC1 = density_func_sncrna4$x, PC4 = density_func_sncrna4$y)
df_density_pos_sncrna4$z <- as.vector(density_func_sncrna4$z)
df_density_neg_sncrna4 <- expand.grid(PC1 = density_neg_sncrna4$x, PC4 = density_neg_sncrna4$y)
df_density_neg_sncrna4$z <- as.vector(density_neg_sncrna4$z)
df_density_pos_sncrna5 <- expand.grid(PC1 = density_func_sncrna5$x, PC5 = density_func_sncrna5$y)
df_density_pos_sncrna5$z <- as.vector(density_func_sncrna5$z)
df_density_neg_sncrna5 <- expand.grid(PC1 = density_neg_sncrna5$x, PC5 = density_neg_sncrna5$y)
df_density_neg_sncrna5$z <- as.vector(density_neg_sncrna5$z)
df_density_pos_sncrna6 <- expand.grid(PC1 = density_func_sncrna6$x, PC6 = density_func_sncrna6$y)
df_density_pos_sncrna6$z <- as.vector(density_func_sncrna6$z)
df_density_neg_sncrna6 <- expand.grid(PC1 = density_neg_sncrna6$x, PC6 = density_neg_sncrna6$y)
df_density_neg_sncrna6$z <- as.vector(density_neg_sncrna6$z)

# Create a factor for sequence_type for shapes
df_density_pos_sncrna2$`Gene type` <- c("sncRNA(+)")
df_density_pos_sncrna2$`Gene type` <- factor(df_density_pos_sncrna2$`Gene type`)
df_density_neg_sncrna2$`Gene type` <- c("sncRNA(-)")
df_density_neg_sncrna2$`Gene type` <- factor(df_density_neg_sncrna2$`Gene type`)

df_density_pos_sncrna3$`Gene type` <- c("sncRNA(+)")
df_density_pos_sncrna3$`Gene type` <- factor(df_density_pos_sncrna3$`Gene type`)
df_density_neg_sncrna3$`Gene type` <- c("sncRNA(-)")
df_density_neg_sncrna3$`Gene type` <- factor(df_density_neg_sncrna3$`Gene type`)

df_density_pos_sncrna4$`Gene type` <- c("sncRNA(+)")
df_density_pos_sncrna4$`Gene type` <- factor(df_density_pos_sncrna4$`Gene type`)
df_density_neg_sncrna4$`Gene type` <- c("sncRNA(-)")
df_density_neg_sncrna4$`Gene type` <- factor(df_density_neg_sncrna4$`Gene type`)

df_density_pos_sncrna5$`Gene type` <- c("sncRNA(+)")
df_density_pos_sncrna5$`Gene type` <- factor(df_density_pos_sncrna5$`Gene type`)
df_density_neg_sncrna5$`Gene type` <- c("sncRNA(-)")
df_density_neg_sncrna5$`Gene type` <- factor(df_density_neg_sncrna5$`Gene type`)

df_density_pos_sncrna6$`Gene type` <- c("sncRNA(+)")
df_density_pos_sncrna6$`Gene type` <- factor(df_density_pos_sncrna6$`Gene type`)
df_density_neg_sncrna6$`Gene type` <- c("sncRNA(-)")
df_density_neg_sncrna6$`Gene type` <- factor(df_density_neg_sncrna6$`Gene type`)

# Define the colors and shapes you want to use
my_colors <- c("sncRNA(-)" = "#56bdfcFF", "sncRNA(+)" = "#D6604DFF")
my_shapes <- c("sncRNA(-)" = 4, "sncRNA(+)" = 17)

# Plot PC distributions
pc1_dist_plot <- ggplot(sncrna_pc_scores, aes(x = Comp.1, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("sncRNA(-)" = "#0491e8", "sncRNA(+)" = "#D6604DFF")) +
  labs(subtitle="sncRNA", x = "PC1 (22.4%)", y = "ECDF") + 
  coord_cartesian(c(-5, 10)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc1_dist_plot

pc2_dist_plot <- ggplot(sncrna_pc_scores, aes(x = Comp.2, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("sncRNA(-)" = "#0491e8", "sncRNA(+)" = "#D6604DFF")) +
  labs(subtitle="sncRNA", x = "PC2 (12.6%)", y = "ECDF") + 
  coord_cartesian(c(-5, 5)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc2_dist_plot

pc3_dist_plot <- ggplot(sncrna_pc_scores, aes(x = Comp.3, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("sncRNA(-)" = "#0491e8", "sncRNA(+)" = "#D6604DFF")) +
  labs(subtitle="sncRNA", x = "PC3 (7.8%)", y = "ECDF") + 
  coord_cartesian(c(-5, 10)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc3_dist_plot

pc4_dist_plot <- ggplot(sncrna_pc_scores, aes(x = Comp.4, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("sncRNA(-)" = "#0491e8", "sncRNA(+)" = "#D6604DFF")) +
  labs(subtitle="sncRNA", x = "PC4 (6.9%)", y = "ECDF") + 
  coord_cartesian(c(-5, 5)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc4_dist_plot

pc5_dist_plot <- ggplot(sncrna_pc_scores, aes(x = Comp.5, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("sncRNA(-)" = "#0491e8", "sncRNA(+)" = "#D6604DFF")) +
  labs(subtitle="sncRNA", x = "PC5 (6.1%)", y = "ECDF") + 
  coord_cartesian(c(-5, 5)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc5_dist_plot

pc6_dist_plot <- ggplot(sncrna_pc_scores, aes(x = Comp.6, color = `Gene type`)) +
  #geom_density(size = 1.5) +
  stat_ecdf(geom = "step", linewidth = 1.5) +
  scale_color_manual(values = c("sncRNA(-)" = "#0491e8", "sncRNA(+)" = "#D6604DFF")) +
  labs(subtitle="sncRNA", x = "PC6 (5.3%)", y = "ECDF") + 
  coord_cartesian(c(-5, 10)) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "lines")
  )

pc6_dist_plot

# Plot PCA
pca_sncrna_plot <- ggplot(sncrna_pc_scores, aes(x = Comp.1, y = Comp.2, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "sncRNA(-)"),
             aes(x = Comp.1, y = Comp.2), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_sncrna2, aes(x = PC1, y = PC2, z = z), 
               bins= 20, color = "#0491e8", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "sncRNA(+)"),
             aes(x = Comp.1, y = Comp.2), shape = 17, size = 4) +
  # Add semi-transparent density contours for sncRNA
  geom_contour(data = df_density_pos_sncrna2, aes(x = PC1, y = PC2, z = z), 
               bins = 20, color = "#471810", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "sncRNA", x = "PC1 (22.4%)", y = "PC2 (12.6%)") +
  theme_minimal()+
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-5, 10) +
  ylim(-5, 5)
pca_sncrna_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_sncrna_20_features_loadings","_pc2.png"), pca_sncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# Plot PCA
pca_sncrna_plot <- ggplot(sncrna_pc_scores, aes(x = Comp.1, y = Comp.3, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "sncRNA(-)"),
             aes(x = Comp.1, y = Comp.3), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_sncrna3, aes(x = PC1, y = PC3, z = z), 
               bins= 20, color = "#0491e8", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "sncRNA(+)"),
             aes(x = Comp.1, y = Comp.3), shape = 17, size = 4) +
  # Add semi-transparent density contours for sncRNA
  geom_contour(data = df_density_pos_sncrna3, aes(x = PC1, y = PC3, z = z), 
               bins = 20, color = "#471810", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "sncRNA", x = "PC1 (22.4%)", y = "PC3 (7.8%)") +
  theme_minimal()+
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-5, 10) +
  ylim(-5, 5)
pca_sncrna_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_sncrna_20_features_loadings","_pc3.png"), pca_sncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# Plot PCA
pca_sncrna_plot <- ggplot(sncrna_pc_scores, aes(x = Comp.1, y = Comp.4, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "sncRNA(-)"),
             aes(x = Comp.1, y = Comp.4), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_sncrna4, aes(x = PC1, y = PC4, z = z), 
               bins= 20, color = "#0491e8", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "sncRNA(+)"),
             aes(x = Comp.1, y = Comp.4), shape = 17, size = 4) +
  # Add semi-transparent density contours for sncRNA
  geom_contour(data = df_density_pos_sncrna4, aes(x = PC1, y = PC4, z = z), 
               bins = 20, color = "#471810", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "sncRNA", x = "PC1 (22.4%)", y = "PC4 (6.9%)") +
  theme_minimal()+
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-5, 10) +
  ylim(-5, 5)
pca_sncrna_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_sncrna_20_features_loadings","_pc4.png"), pca_sncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# Plot PCA
pca_sncrna_plot <- ggplot(sncrna_pc_scores, aes(x = Comp.1, y = Comp.5, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "sncRNA(-)"),
             aes(x = Comp.1, y = Comp.5), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_sncrna5, aes(x = PC1, y = PC5, z = z), 
               bins= 20, color = "#0491e8", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "sncRNA(+)"),
             aes(x = Comp.1, y = Comp.5), shape = 17, size = 4) +
  # Add semi-transparent density contours for sncRNA
  geom_contour(data = df_density_pos_sncrna5, aes(x = PC1, y = PC5, z = z), 
               bins = 20, color = "#471810", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "sncRNA", x = "PC1 (22.4%)", y = "PC5 (6.1%)") +
  theme_minimal()+
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-5, 10) +
  ylim(-5, 5)
pca_sncrna_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_sncrna_20_features_loadings","_pc5.png"), pca_sncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# Plot PCA
pca_sncrna_plot <- ggplot(sncrna_pc_scores, aes(x = Comp.1, y = Comp.6, shape = `Gene type`, color = `Gene type`)) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "sncRNA(-)"),
             aes(x = Comp.1, y = Comp.6), shape = 4, size = 4) +
  # Add semi-transparent density contours
  geom_contour(data = df_density_neg_sncrna6, aes(x = PC1, y = PC6, z = z), 
               bins= 20, color = "#0491e8", alpha = 0.6, linewidth = 1) +
  geom_point(data = subset(sncrna_pc_scores, `Gene type` == "sncRNA(+)"),
             aes(x = Comp.1, y = Comp.6), shape = 17, size = 4) +
  # Add semi-transparent density contours for sncRNA
  geom_contour(data = df_density_pos_sncrna6, aes(x = PC1, y = PC6, z = z), 
               bins = 20, color = "#471810", alpha = 0.6, linewidth = 1) +
  # Use scale_manual functions
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "sncRNA", x = "PC1 (22.4%)", y = "PC6 (5.3%)") +
  theme_minimal()+
  theme(
    plot.subtitle = element_text(size = 38, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 30),  # Increase axis title size
    axis.text = element_text(size = 26),    # Increase axis label size
    legend.position = "right",
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.text = element_text(size = 28),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  xlim(-5, 10) +
  ylim(-5, 5)
pca_sncrna_plot
ggsave(paste0("results/pca/9-Dec-2025/pca_sncrna_20_features_loadings","_pc6.png"), pca_sncrna_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# Change colnames
#colnames(sncrna_data_normalized) <- c(PCA_20_SELECT_FEATURES_LABELS, "Dataset")
# Compute contribution of features to variation of PCAs
#sncrna_pca <- princomp(sncrna_data_normalized |> dplyr::select(-Dataset) |> scale())
#summary(sncrna_pca)

# Inspect variable loadings
#sncrna_pca$loadings[, 1:3]

# Scree Plot
#p1 <- fviz_eig(sncrna_pca, addlabels = TRUE, title = "Scree plot - sncRNA")
#p1_mod <- p1 + labs(x = "Principal Components") + theme(text = element_text(size = 26))
#p1_mod$layers[[4]]$aes_params$size <- 8
#print(p1_mod)

# Graph of the variables
#fviz_pca_var(sncrna_pca, col.var = "black")

# Contribution of each variable
#fviz_cos2(sncrna_pca, choice = "var", axes = 1:2, title = "Contribution of Features to PC 1 and 2 Variance - sncRNA")

# 1. Create the plot object, but don't print it yet
#p <- fviz_cos2(sncrna_pca, choice = "var", axes = 5, title = "Contribution of Features to PC 5 Variance - sncRNA")

# 2. Modify the plot object's aesthetics
# The labels are text, so we target the 'text' aesthetic.
# We set a new, larger size for all text elements.
#p_modified <- p + labs(y = "Contribution to PC5 Variance") + theme(text = element_text(size = 26),
#                        axis.text.x = element_text(angle = 90, vjust = 0.5)) # Change 16 to your desired size

# 3. Print the modified plot
#print(p_modified)

# Biplot
#fviz_pca_var(sncrna_pca, col.var = "cos2",
#             gradient.cols = c("black", "orchid", "blue"),
#             repel = TRUE,
#             title = "Biplot - sncRNA",
#             xlab = "PC1 (26.5%)",
#             ylab = "PC2 (12.3%)"
#             ) +
#  theme(
#    text = element_text(size = 26)
#  )
  

################################################
# --- Join the plots in a single patchwork --- #
pca_joined_plot <- (pca_protein_plot + pca_sncrna_plot + pca_lncrna_plot) +
  plot_annotation(title = "Principal Component Analysis",
                  tag_levels = list(c('A','B','C'))) +
  plot_layout(axis_titles = "collect", guides = "collect") &
  theme(plot.title = element_text(size = 56, hjust = 0.5, margin = ggplot2::margin(0,0,40,0)),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 34, face = "bold", hjust = 0, vjust = 0))
pca_joined_plot
ggsave(PCA_20_FEATURES_LOADINGS_JOINED_PLOT_FILE, pca_loadings_joined_plot, scale = 3, width = 3840, height = 1620, units = "px", bg = "white", dpi = 600)





