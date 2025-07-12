library(dplyr)
install.packages("MASS")     # For the kde2d function
library(MASS)
library(ggplot2)
#### PCA PLOT ####
# Load Data #
source("load_gene_functionality_features.R")
source("load_gene_functionality_features.R")

# List selected features #
pca_select_features <- c("GC_percentage", "CpG", "GA", "TA", 
                         "phyloP_max_241w", "phyloP_max_100w", 
                         "GERP_91_mammals_max", "GERP_63_amniotes_max", 
                         "RPKM_tissue", "RPKM_primary.cell", 
                         "H3K27ac_MaxScaledSignal", "H3K36me3_MaxScaledSignal", "H3K79me2_MaxScaledSignal", 
                         "chrm_acc_MaxScaledSignal", "methylome", 
                         "copy_number", "coding_potential", 
                         "Max_covariance", "MFE", 
                         "Interaction_ave"
)

pca_select_features_labels <- c("GC%", "CpG",  "GA", "TA",
                                "PhyloP-mammals", "PhyloP-vertebrates",
                                "GERP-mammals", "GERP-vertebrates",
                                "Tissue RPKM", "Primary cell RPKM",
                                "H3K27ac", "H3K36me3", "H3K79me2",
                                "Chromatin", "Methylome",
                                "Copies", "RNAcode",
                                "Covariance", "MFE",
                                "Interactions"
)
###########

# Perform PCA by Sequence
help("scale")
normalized_fmn <- scale(data_numeric[,pca_select_features])
cor_matrix <- cor(normalized_fmn, method = "spearman")
#pca_result <- prcomp(feature_matrix_numeric_no_nas, scale. = TRUE, center = TRUE, rank. = 5)
pca_result <- prcomp(normalized_fmn, scale. = TRUE)
pca_result <- prcomp(cor_matrix, scale. = TRUE)
#pca_df <- as.data.frame(pca_result$x)
#pc_scores <- as.data.frame(pca_result$x[, 1:2])
summary(pca_result)
pc_scores <- pca_result$x
pca_df <- cbind(pc_scores, feature_matrix_numeric_no_nas)


# Create a factor for sequence_type for shapes
unique(dataset_column_sel)
pca_df$`Gene type` <- factor(feature_matrix_no_nas$Dataset,
                             levels = c("protein-coding-exon2","protein-coding-exon3",
                                        "lncrna-exon1", "lncrna-exon2",
                                        "short-ncrna",
                                        "protein-exon2-negative-control",
                                        "protein-exon3-negative-control",
                                        "lncrna-exon1-negative-control",
                                        "lncrna-exon2-negative-control",
                                        "short-ncrna-negative-control"),
                             labels = c("Protein coding","Protein coding",
                                        "lncRNA","lncRNA",
                                        "sncRNA",
                                        "Negative controls",
                                        "Negative controls",
                                        "Negative controls",
                                        "Negative controls",
                                        "Negative controls")
                              )
custom_order <- c("Protein coding", "sncRNA", "lncRNA", "Negative controls")
pca_df$`Gene type` <- factor(pca_df$`Gene type`, levels = custom_order)

# PCA regression
var_explained <- pca_result$sdev^2
pct_contribution <- var_explained / sum(var_explained) * 100

print(pct_contribution[1:15]) # Print contribution of the first 15 PCs
pca_df <- na.omit(pca_df)
str(pca_df)
model <- lm(as.numeric(pca_df$sequence_type) ~ pc_scores, data = feature_matrix_numeric)
summary(model)

install.packages("factoextra")
library(factoextra)
fviz_pca_biplot(pca_result, label = "var", habillage = pca_df$sequence_type, addEllipses = TRUE, repel = TRUE,
                ellipse.level=0.95, col.var = "darkblue", palette = c("firebrick","darkorange","limegreen","steelblue3"),
                ggtheme = theme_minimal()) +
  scale_shape_manual(values = c(15, 19, 17, 4))

help(fviz_pca_biplot)
fviz_pca_var(pca_result, label = "var")
###################

# Perform PCA by Feature
# Normalize data
data_normalized <- scale(feature_matrix_numeric_t)

pca_result_by_feature <- prcomp(data_normalized, scale. = TRUE)
pc_scores_by_feature <- as.data.frame(pca_result_by_feature$x[, 1:2])
pca_by_feature_df <- cbind(pc_scores_by_feature)

# Add Row names as factor to the PCA by Feature result data frame
pca_by_feature_df$feature_name <- rownames(feature_matrix_numeric_t)
pca_by_feature_df$feature_name <- factor(pca_by_feature_df$feature_name, levels = c(rownames(feature_matrix_numeric_t)))

###################

# Plot PCA by Sequence
ggplot(pca_df, aes(x = PC1, y = PC2, shape = sequence_type, color = sequence_type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(15, 19, 17, 4)) +  #c(19, 17, 15))
  scale_color_manual(values = c("firebrick","darkorange","limegreen","blue")) +  # Specify colors
  labs(title = "Principal Component Analysis by Sequence", x = "PC1", y = "PC2") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10),  # Increase title size
    axis.title = element_text(size = 10),  # Increase axis title size
    axis.text = element_text(size = 10),    # Increase axis label size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 18),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  # Add semi-transparent density contours
  geom_contour(data = df_density, aes(x = PC1, y = PC2, z = z), 
               bins = 10, color = "black", alpha = 1) 



# Plot PCA by Feature
ggplot(pca_by_feature_df, aes(x = PC1, y = PC2, color = feature_name, shape = feature_name)) +
  geom_point(size = 4, position = position_dodge(width = 0.1)) +
  scale_shape_manual(values = c(2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,1,2,3,4,5,6,7,8,9,10,11)) +  # Specify shapes
  scale_color_manual(values = c("green","pink4","orange4","brown","purple","darkgreen","darkorchid4",
                                "orchid4","lightcoral","lightblue4","orange3","maroon4","purple3","khaki4","limegreen","hotpink4",
                                "gold4","forestgreen","firebrick","dodgerblue4","deeppink4","darkviolet","darkslategrey","darksalmon","darkolivegreen4","darkcyan",
                                "yellowgreen","wheat4","violetred4","turquoise4","tomato4","chartreuse4","cadetblue4","brown","chocolate","cornflowerblue")) +  # Specify colors
  labs(title = "Principal Component Analysis by Feature", x = "PC1", y = "PC2") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10),  # Increase title size
    axis.title = element_text(size = 10),  # Increase axis title size
    axis.text = element_text(size = 10),    # Increase axis label size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 18),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  )

########################
# 2D contour density plot by Sequence
lims <- c(range(pca_df$PC1), range(pca_df$PC2))
density <- kde2d(pca_df$PC1, pca_df$PC2, n = 100, lims = lims)
#help(kde2d)
filled.contour(density, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")

df_density <- expand.grid(PC1 = density$x, PC2 = density$y)
df_density$z <- as.vector(density$z)

# Create a factor for sequence_type for shapes
df_density$`Gene type` <- factor(df_density$z, 
                                   levels = c("Protein coding","lncRNA","sncRNA","Negative controls"))

plot_data <- with(density, expand.grid(PC1 = density$x, PC2 = density$y))
plot_data$z <- as.vector(density$z)

ggplot(df_density, aes(x = PC1, y = PC2, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c() +
  geom_contour(aes(z = z), bins = 10) +
  labs(title = "Density Contour Plot of PCA density by Sequence", x = "PC1", y = "PC2")


ggplot(plot_data, aes(x = PC1, y = PC2, z = z)) +
  geom_contour(aes(z = z), bins = 10) +
  scale_fill_viridis_c() +
  labs(title = "2D Density Contour of PCA Scores", x = "PC1", y = "PC2") +
  theme_minimal()
###################


###################
#Protein Coding PCA
protein_feature_matrix_numeric <- data_numeric[data_numeric$Dataset=="protein-coding-exon2" | data_numeric$Dataset=="protein-coding-exon3", pca_select_features]
protein_feature_matrix_numeric <- protein_feature_matrix_numeric %>% na.omit()
# Perform PCA by Sequence
protein_pca_result <- prcomp(protein_feature_matrix_numeric, scale. = TRUE, center = TRUE, rank. = 5)
protein_pc_scores <- as.data.frame(protein_pca_result$x[, 1:2])
protein_pca_df <- cbind(protein_pc_scores, protein_feature_matrix_numeric[, !names(protein_feature_matrix_numeric) %in% c("Functional", "Dataset")])

# Create a factor for sequence_type for shapes
protein_pca_df$`Gene type` <- factor(c("Protein coding"), 
                               levels = c("Protein coding"))

########################
# 2D contour density plot by Sequence
protein_density <- kde2d(pca_df[pca_df$`Gene type` == "Protein coding",]$PC1,
                         pca_df[pca_df$`Gene type` == "Protein coding",]$PC2, 
                         n = 100, 
                         lims = lims)
filled.contour(protein_density, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")
#help("filled.contour")
protein_df_density <- expand.grid(PC1 = protein_density$x, PC2 = protein_density$y)
protein_df_density$z <- as.vector(protein_density$z)

# Create a factor for sequence_type for shapes
sort(unique(protein_df_density$z))
protein_df_density$`Gene type` <- factor(c("Protein coding"), 
                                           levels = c("Protein coding"))

ggplot(protein_df_density, aes(x = PC1, y = PC2, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c() +
  geom_contour(aes(z = z), bins = 0) +
  labs(title = "Density Contour Plot of PCA density by Sequence", x = "PC1", y = "PC2")



# Plot PCA by Sequence
ggplot(pca_df[pca_df$Dataset == "protein-coding",], aes(x = PC1, y = PC2, shape = sequence_type, color = sequence_type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(15, 19, 17, 4)) +  #c(19, 17, 15))
  scale_color_manual(values = c("firebrick","darkorange","limegreen","blue")) +  # Specify colors
    labs(title = "Principal Component Analysis by Sequence", x = "PC1", y = "PC2") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10),  # Increase title size
    axis.title = element_text(size = 10),  # Increase axis title size
    axis.text = element_text(size = 10),    # Increase axis label size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 18),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  # Add semi-transparent density contours
  geom_contour(data = protein_df_density, aes(x = PC1, y = PC2, z = z), 
               bins = 10, color = "firebrick4", alpha = 1) +
  scale_fill_viridis_c()
  


########################

########################
#Long non-coding RNA PCA
matrix_lncrna <- na.omit(matrix_lncrna)
lncrna_functional_column <- matrix_lncrna$Functional
lncrna_dataset_column <- matrix_lncrna$Dataset
lncrna_feature_matrix <- matrix_lncrna[, names(matrix_lncrna) %in% c("Start","copy_number","GC.","CG","GA","GG","TA","TC","X241w_PP_max",
                                                                     "X100w_PP_max","RPKM_tissue","RPKM_primary.cell","Dfam_sum","RNAcode_score",
                                                                     "Fickett_score","Max_covariance","MFE","Accessibility","RNAalifold_score",
                                                                     "Interaction_ave","gnomAD_SNP_density","gnomAD_aveMAF","H3K36me3_AvgSignal",
                                                                     "H3K27ac_AvgSignal","H3K79me2_AvgSignal","chrm_acc_AvgSignal","methylome")]
lncrna_feature_matrix_numeric <- lncrna_feature_matrix[, sapply(lncrna_feature_matrix, is.numeric)] # Matrix for PCA by Sequence
lncrna_feature_matrix_numeric_t <- t(lncrna_feature_matrix_numeric) # Matrix for PCA by Feature without segregation by sequence type


# Perform PCA by Sequence
lncrna_pca_result <- prcomp(lncrna_feature_matrix_numeric, scale. = TRUE, center = TRUE, rank. = 5)
lncrna_pc_scores <- as.data.frame(lncrna_pca_result$x[, 1:2])
lncrna_pca_df <- cbind(lncrna_pc_scores, matrix_lncrna[, !names(matrix_lncrna) %in% c("Functional", "Dataset")])

# Add 'Functional' and 'Dataset' columns back to the PCA by Sequence result data frame
lncrna_pca_df$Functional <- lncrna_functional_column
lncrna_pca_df$Dataset <- lncrna_dataset_column

# Create a factor for sequence_type for shapes
lncrna_pca_df$sequence_type <- factor(lncrna_pca_df$Dataset, 
                                       levels = c("lncrna"))

########################

# 2D contour density plot by Sequence
lncrna_density <- kde2d(pca_df[pca_df$`Gene type` == "lncRNA",]$PC1, 
                        pca_df[pca_df$`Gene type` == "lncRNA",]$PC2, 
                        n = 100, 
                        lims = lims)
filled.contour(lncrna_density, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")

lncrna_df_density <- expand.grid(PC1 = lncrna_density$x, PC2 = lncrna_density$y)
lncrna_df_density$z <- as.vector(lncrna_density$z)

# Create a factor for sequence_type for shapes
lncrna_df_density$`Gene type` <- factor(c("lncRNA"), 
                                           levels = c("lncRNA"))

ggplot(lncrna_df_density, aes(x = PC1, y = PC2, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c() +
  geom_contour(aes(z = z), bins = 10) +
  labs(title = "Density Contour Plot of PCA density by Sequence", x = "PC1", y = "PC2")



# Plot PCA by Sequence
ggplot(pca_df[pca_df$Dataset == "lncrna",], aes(x = PC1, y = PC2, shape = sequence_type, color = sequence_type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(19)) +  #c(19, 17, 15))
  scale_color_manual(values = c("darkorange")) +  # Specify colors
  labs(title = "Principal Component Analysis by Sequence", x = "PC1", y = "PC2") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10),  # Increase title size
    axis.title = element_text(size = 10),  # Increase axis title size
    axis.text = element_text(size = 10),    # Increase axis label size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 18),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  # Add semi-transparent density contours
  geom_contour(data = lncrna_df_density, aes(x = PC1, y = PC2, z = z), 
               bins = 10, color = "orange4", alpha = 1)


#########################

#########################
#Short non-coding RNA PCA
matrix_short_ncrna <- na.omit(matrix_short_ncrna)
short_ncrna_functional_column <- matrix_short_ncrna$Functional
short_ncrna_dataset_column <- matrix_short_ncrna$Dataset
short_ncrna_feature_matrix <- matrix_short_ncrna[, names(matrix_short_ncrna) %in% c("Start","copy_number","GC.","CG","GA","GG","TA","TC","X241w_PP_max",
                                                                                    "X100w_PP_max","RPKM_tissue","RPKM_primary.cell","Dfam_sum","RNAcode_score",
                                                                                    "Fickett_score","Max_covariance","MFE","Accessibility","RNAalifold_score",
                                                                                    "Interaction_ave","gnomAD_SNP_density","gnomAD_aveMAF","H3K36me3_AvgSignal",
                                                                                    "H3K27ac_AvgSignal","H3K79me2_AvgSignal","chrm_acc_AvgSignal","methylome")]
short_ncrna_feature_matrix_numeric <- short_ncrna_feature_matrix[, sapply(short_ncrna_feature_matrix, is.numeric)] # Matrix for PCA by Sequence
short_ncrna_feature_matrix_numeric_t <- t(short_ncrna_feature_matrix_numeric) # Matrix for PCA by Feature without segregation by sequence type


# Perform PCA by Sequence
short_ncrna_pca_result <- prcomp(short_ncrna_feature_matrix_numeric, scale. = TRUE, center = TRUE, rank. = 5)
short_ncrna_pc_scores <- as.data.frame(short_ncrna_pca_result$x[, 1:2])
short_ncrna_pca_df <- cbind(short_ncrna_pc_scores, matrix_short_ncrna[, !names(matrix_short_ncrna) %in% c("Functional", "Dataset")])

# Add 'Functional' and 'Dataset' columns back to the PCA by Sequence result data frame
short_ncrna_pca_df$Functional <- short_ncrna_functional_column
short_ncrna_pca_df$Dataset <- short_ncrna_dataset_column

# Create a factor for sequence_type for shapes
short_ncrna_pca_df$sequence_type <- factor(short_ncrna_pca_df$Dataset, 
                                       levels = c("sincrna"))

########################
# 2D contour density plot by Sequence
short_ncrna_density <- kde2d(pca_df[pca_df$`Gene type` == "sncRNA",]$PC1, 
                             pca_df[pca_df$`Gene type` == "sncRNA",]$PC2, 
                             n = 100, 
                             lims = lims)
filled.contour(short_ncrna_density, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density")

short_ncrna_df_density <- expand.grid(PC1 = short_ncrna_density$x, PC2 = short_ncrna_density$y)
short_ncrna_df_density$z <- as.vector(short_ncrna_density$z)

# Create a factor for sequence_type for shapes
short_ncrna_df_density$`Gene type` <- factor(c("sncRNA"), 
                                           levels = c("sncRNA"))

ggplot(short_ncrna_df_density, aes(x = PC1, y = PC2, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c() +
  geom_contour(aes(z = z), bins = 10) +
  labs(title = "Density Contour Plot of PCA density by Sequence", x = "PC1", y = "PC2")



# Plot PCA by Sequence
ggplot(pca_df[pca_df$Dataset == "short-ncrna",], aes(x = PC1, y = PC2, shape = sequence_type, color = sequence_type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(17)) + 
  scale_color_manual(values = c("limegreen")) +  # Specify colors
  labs(title = "Principal Component Analysis by Sequence", x = "PC1", y = "PC2") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10),  # Increase title size
    axis.title = element_text(size = 10),  # Increase axis title size
    axis.text = element_text(size = 10),    # Increase axis label size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 18),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  # Add semi-transparent density contours
  geom_contour(data = short_ncrna_df_density, aes(x = PC1, y = PC2, z = z), 
               bins = 10, color = "green4", alpha = 1)

###############

#########################
#Negative controls PCA
matrix_negative_control <- na.omit(matrix_negative_control)
negative_control_functional_column <- matrix_negative_control$Functional
negative_control_dataset_column <- matrix_negative_control$Dataset
negative_control_feature_matrix <- matrix_negative_control[, names(matrix_negative_control) %in% c("Start","copy_number","GC.","CG","GA","GG","TA","TC","X241w_PP_max",
                                                                                                   "X100w_PP_max","RPKM_tissue","RPKM_primary.cell","Dfam_sum","RNAcode_score",
                                                                                                   "Fickett_score","Max_covariance","MFE","Accessibility","RNAalifold_score",
                                                                                                   "Interaction_ave","gnomAD_SNP_density","gnomAD_aveMAF","H3K36me3_AvgSignal",
                                                                                                   "H3K27ac_AvgSignal","H3K79me2_AvgSignal","chrm_acc_AvgSignal","methylome")]
negative_control_feature_matrix_numeric <- negative_control_feature_matrix[, sapply(negative_control_feature_matrix, is.numeric)] # Matrix for PCA by Sequence
negative_control_feature_matrix_numeric_t <- t(negative_control_feature_matrix_numeric) # Matrix for PCA by Feature without segregation by sequence type


# Perform PCA by Sequence
negative_control_pca_result <- prcomp(negative_control_feature_matrix_numeric, scale. = TRUE, center = TRUE, rank. = 5)
negative_control_pc_scores <- as.data.frame(negative_control_pca_result$x[, 1:2])
negative_control_pca_df <- cbind(negative_control_pc_scores, matrix_negative_control[, !names(matrix_negative_control) %in% c("Functional", "Dataset")])

# Add 'Functional' and 'Dataset' columns back to the PCA by Sequence result data frame
negative_control_pca_df$Functional <- negative_control_functional_column
negative_control_pca_df$Dataset <- negative_control_dataset_column

# Create a factor for sequence_type for shapes
negative_control_pca_df$sequence_type <- factor(negative_control_pca_df$Dataset, 
                                           levels = c("negative-control"))

########################
# 2D contour density plot by Sequence
negative_control_density <- kde2d(pca_df[pca_df$`Gene type` == "Negative controls",]$PC1, 
                                  pca_df[pca_df$`Gene type` == "Negative controls",]$PC2, 
                                  n = 100, 
                                  lims = lims)
filled.contour(negative_control_density, xlab = "PC1", ylab = "PC2", main = "Filled Contour Plot of PCA Density", nlevels = 25)

negative_control_df_density <- expand.grid(PC1 = negative_control_density$x, PC2 = negative_control_density$y)
negative_control_df_density$z <- as.vector(negative_control_density$z)

# Create a factor for sequence_type for shapes
negative_control_df_density$`Gene type` <- factor(c("Negative controls"), 
                                               levels = c("Negative controls"))

ggplot(negative_control_df_density, aes(x = PC1, y = PC2, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c() +
  geom_contour(aes(z = z), bins = 10) +
  labs(title = "Density Contour Plot of PCA density by Sequence", x = "PC1", y = "PC2")



# Plot PCA by Sequence
ggplot(pca_df[pca_df$Dataset == "negative-control",], aes(x = PC1, y = PC2, shape = sequence_type, color = sequence_type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(4)) + 
  scale_color_manual(values = c("navy")) +  # Specify colors
  labs(title = "Principal Component Analysis by Sequence", x = "PC1", y = "PC2") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10),  # Increase title size
    axis.title = element_text(size = 10),  # Increase axis title size
    axis.text = element_text(size = 10),    # Increase axis label size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 18),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  # Add semi-transparent density contours
  geom_contour(data = negative_control_df_density, aes(x = PC1, y = PC2, z = z), 
               bins = 10, color = "blue2", alpha = 1)


###############
# Join all PCAs
all_pca_df <- rbind(protein_pca_df,lncrna_pca_df,short_ncrna_pca_df,negative_control_pca_df)
all_density_df <- rbind(protein_df_density,lncrna_df_density,short_ncrna_df_density,negative_control_df_density)


# Plot PCA by Sequence
ggplot(pca_df, aes(x = PC1, y = PC2, shape = `Gene type`, color = `Gene type`)) +
  geom_point(size = 2, alpha = 0.4) +
  scale_shape_manual(values = c(15,17,19,4)) +  #c(19, 17, 15))
  scale_color_manual(values = c("firebrick","limegreen","darkorange","blue")) +  # Specify colors
  labs(title = "Principal Component Analysis by Sequence", x = "PC1 (16.15%)", y = "PC2 (9.80%)") +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 34),  # Increase title size
    axis.title = element_text(size = 28),  # Increase axis title size
    axis.text = element_text(size = 24),    # Increase axis label size
    legend.title = element_text(size = 28),  # Increase legend title size
    legend.text = element_text(size = 26),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
    
  ) +
  # Add semi-transparent density contours for lncRNA
  geom_contour(data = lncrna_df_density, aes(x = PC1, y = PC2, z = z), 
               bins = 80, color = "darkorange", alpha = 0.6, linewidth = 1) +
  # Add semi-transparent density contours for protein coding
  geom_contour(data = protein_df_density, aes(x = PC1, y = PC2, z = z), 
               bins = 80, color = "firebrick4", alpha = 0.6, linewidth = 1) +
  # Add semi-transparent density contours for sncRNA
  geom_contour(data = short_ncrna_df_density, aes(x = PC1, y = PC2, z = z), 
               bins = 80, color = "limegreen", alpha = 0.6, linewidth = 1) +
  
  # Add semi-transparent density contours
  geom_contour(data = negative_control_df_density, aes(x = PC1, y = PC2, z = z), 
               bins= 30, color = "blue2", alpha = 0.4, linewidth = 1)
help("geom_contour")
ggsave("2dcontour.png",path = "../results/latest1000all/violinPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


install.packages("rgl")  # For interactive 3D plots
install.packages("pca3d") # Helpful for PCA-specific 3D visualization
library(rgl)
# Extract PC scores for the first 3 components
pc_scores <- pca_result$x[, 1:3]
summary(pca_result)
# Example with 3 factor levels, adjust as needed
factor_levels <- levels(pca_df$`Gene type`)
group_colors <- c("firebrick4","darkorange","green", "blue") 

# Create a named color vector for easier mapping
color_mapping <- setNames(group_colors, factor_levels)
# Create the 3D plot
plot3d(pca_df,
       xlab = "PC1 (14.21%)", ylab = "PC2 (10.36%)", zlab = "PC3 (7.49%)",
       col = color_mapping[pca_df$`Gene type`], 
       size = 3) 
# Add legend
legend3d("topright", legend = factor_levels, pch = 16, col = group_colors)
# Specify the filename and desired resolution (e.g., 300 dots per inch)
help("rgl.snapshot")
rgl.snapshot("pca_plot_3d.png", fmt = "png", top = TRUE, webshot = FALSE, width = 1200, height = 900)
