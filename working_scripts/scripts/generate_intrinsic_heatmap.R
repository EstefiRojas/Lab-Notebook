# Generate a heatmap of intrinsic features correlation using a subset of intrinsic
# data. This script selects 100 positive cases and 1000 negative controls for each
# RNA type and their respective exons analyzed. To replicate the plots published in 
# the paper, all intrinsic data should be selected (i.e. no sub-setting).

# Load necessary libraries
library(Hmisc) # For rcorr() function
library(ggcorrplot)
library(dplyr)
library(reshape2)

library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
#install.packages("ggplotify")
library(ggplotify) ## to convert pheatmap to ggplot2
#install.packages("heatmaply")
library(heatmaply) ## for constructing interactive heatmap
library(scales)
library(plotly)   # For combining the final plots


#######################################
# Load KS statistics utility functions. (A custom version of ks to keep the sign in the D statistic)
source("scripts/config.R")
source("scripts/utils.R")

source("scripts/load_dinucleotide_features.R")
feature_matrix <- load_dinucleotide_features()
colnames(feature_matrix)
###############################
# Define features of interest #
SELECTED_INTRINSIC_FEATURES <- c("AA", "AC", "AG",
                               "AT", "CA", "CC",
                               "CpG", "CT", "GA", 
                               "GC",
                               "GG", "GT", "TA",
                               "TC", "TG", "TT",
                               "GC_percentage", "lowComplexity_density")

SELECTED_INTRINSIC_LABELS <- c("AA", "AC", "AG",
                             "AT", "CA", "CC",
                             "CpG", "CT", "GA", 
                             "GC",
                             "GG", "GT", "TA",
                             "TC", "TG", "TT",
                             "GC%", "Complexity")

# Select din features and convert to numeric
data_numeric <- feature_matrix %>% 
  dplyr::select(all_of(SELECTED_INTRINSIC_FEATURES)) %>%
  sapply(function(feature) as.numeric(as.character(feature))) %>%
  as.data.frame()
# Restore Dataset column
data_numeric$Dataset <- feature_matrix$Dataset

# Select a subset per RNA type. Choose 100 for positive cases, 1000 for negative controls.
set.seed(2025)  # Set seed for reproducibility
data_numeric_sample <- rbind(
  data_numeric %>%
    filter(Dataset %in% c("protein-coding-exon2")) %>%
    sample_n(100),
  data_numeric %>%
    filter(Dataset %in% c("protein-coding-exon3")) %>%
    sample_n(100),
  data_numeric %>%
    filter(Dataset %in% c("lncrna-exon1")) %>%
    sample_n(100),
  data_numeric %>%
    filter(Dataset %in% c("lncrna-exon2")) %>%
    sample_n(100),
  data_numeric %>%
    filter(Dataset %in% c("short-ncrna")) %>%
    sample_n(100),
  data_numeric %>%
    filter(Dataset %in% c("protein-exon2-negative-control")) %>%
    sample_n(1000),
  data_numeric %>%
    filter(Dataset %in% c("protein-exon3-negative-control")) %>%
    sample_n(1000),
  data_numeric %>%
    filter(Dataset %in% c("lncrna-exon1-negative-control")) %>%
    sample_n(1000),
  data_numeric %>%
    filter(Dataset %in% c("lncrna-exon2-negative-control")) %>%
    sample_n(1000),
  data_numeric %>%
    filter(Dataset %in% c("short-ncrna-negative-control")) %>%
    sample_n(1000)
)
# If the paper plot is desired, un-comment the following line, no other change needed.
#data_numeric_sample <- data_numeric

######################################
# Compute ks-stats for raw features. #
########
# mrna #
subset_m <- data_numeric_sample %>%
  filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")) %>%
  dplyr::select(-Dataset)
subset_nc_m <- data_numeric_sample %>%
  filter(Dataset %in% c("protein-exon2-negative-control", "protein-exon3-negative-control")) %>%
  dplyr::select(-Dataset)

# Run the K-S tests
ks_results_m <- run_ks_tests(subset_nc_m, subset_m)

#print(ks_results_m)
write.csv(ks_results_m, INTRINSIC_MRNA_KS_STAT_SAMPLE_FILE)

ks_d_values_m <- as.data.frame(t(ks_results_m))
rownames(ks_d_values_m) <- SELECTED_INTRINSIC_LABELS

##########
# lncrna #
subset_l <- data_numeric_sample %>%
  filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2")) %>%
  dplyr::select(-Dataset)
subset_nc_l <- data_numeric_sample %>%
  filter(Dataset %in% c("lncrna-exon1-negative-control","lncrna-exon2-negative-control")) %>%
  dplyr::select(-Dataset)

# Run the K-S tests
ks_results_l <- run_ks_tests(subset_nc_l, subset_l)

#print(ks_results_l)
write.csv(ks_results_l, INTRINSIC_LNCRNA_KS_STAT_SAMPLE_FILE)

ks_d_values_l <- as.data.frame(t(ks_results_l))
rownames(ks_d_values_l) <- SELECTED_INTRINSIC_LABELS

##########
# sncrna #
subset_s <- data_numeric_sample %>%
  filter(Dataset %in% c("short-ncrna")) %>%
  dplyr::select(-Dataset)
subset_nc_s <- data_numeric_sample %>%
  filter(Dataset %in% c("short-ncrna-negative-control")) %>%
  dplyr::select(-Dataset)

# Run the K-S tests
ks_results_s <- run_ks_tests(subset_nc_s, subset_s)

#print(ks_results_s)
write.csv(ks_results_s, INTRINSIC_SNCRNA_KS_STAT_SAMPLE_FILE)

ks_d_values_s <- as.data.frame(t(ks_results_s))
rownames(ks_d_values_s) <- SELECTED_INTRINSIC_LABELS

##################
# mean|signed_D| #
ks_results <- data.frame(
  meanAbsD = (abs(ks_d_values_m$signed_D) + abs(ks_d_values_s$signed_D) + abs(ks_d_values_l$signed_D))/3
)
rownames(ks_results) <- SELECTED_INTRINSIC_LABELS
#print(ks_results)
write.csv(ks_results, INTRINSIC_MEAN_KS_STAT_SAMPLE_FILE)


#########################################################
# Select data for analysis and run Spearman correlation #
allrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric_sample %>%
      #filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3",
      #                      "lncrna-exon1", "lncrna-exon2",
      #                      "short-ncrna")) %>%
      select(all_of(SELECTED_INTRINSIC_FEATURES))
  ),
  type = "spearman"
)

# Unpack correlation coefficients and p-values #
allrna_corr_matrix <- allrna_corr_obj[["r"]] # Don't compute the absolute value of the spearman coeff.
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

# Adjust p-values using Benjamini-Hochberg
#help("p.adjust")
allrna_p_values_adjusted <- p.adjust(allrna_pvalues, method="BH")

allrna_p_values_adjusted <- matrix(allrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(allrna_p_values_adjusted) <- diag

# Use short names
rownames(allrna_corr_matrix) <- SELECTED_INTRINSIC_LABELS
colnames(allrna_corr_matrix) <- SELECTED_INTRINSIC_LABELS

# Use short names
rownames(allrna_p_values_adjusted) <- SELECTED_INTRINSIC_LABELS
colnames(allrna_p_values_adjusted) <- SELECTED_INTRINSIC_LABELS


##############################
# Version 1 of heatmap plots, using ggcorrplot package #
ggcorrplot(allrna_corr_matrix, type = "lower", method = "square", lab = TRUE, 
           lab_col = "black", lab_size = 3, ggtheme = theme_void(), 
           title = "Spearman correlation Heatmap - All RNA types", show.diag = TRUE, 
           p.mat = allrna_p_values_adjusted, sig.level = corrected_alpha, insig = "pch"
) +
  theme(
    plot.title = element_text(size = 24),  # Increase title size
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    
  )

# Try in-built heatmap package from R #
heatmap(allrna_corr_matrix)


###################################
# Version 2, try pheatmap package #
pheatmap(allrna_corr_matrix, scale="row",
         
         color=colorRampPalette(c("navy", "white", "red"))(50),
         #cutree_cols=2, cutree_rows=2, 
         main="Spearman correlation Heatmap - All RNA types")

####################################
# Version 3, use heatmaply package #
# Calculate hierarchical clustering for rows and columns
hc_rows <- hclust(dist(abs(allrna_corr_matrix))) # Compute absolute values only for clustering.
hc_cols <- hclust(dist(t(abs(allrna_corr_matrix)))) # Compute absolute values only for clustering.

# Create dendrograms
row_dendro <- as.dendrogram(hc_rows)
col_dendro <- as.dendrogram(hc_cols)

# Get the column order from the clustering to align our plots later
col_order <- order.dendrogram(col_dendro)
row_order <- order.dendrogram(row_dendro)

# 3. Create the matrix for cell notes
cellnote_matrix <- matrix(NA_character_, nrow = nrow(allrna_corr_matrix), ncol = ncol(allrna_corr_matrix))
if (!is.null(dimnames(allrna_corr_matrix))) {
  dimnames(cellnote_matrix) <- dimnames(allrna_corr_matrix)
}

for (r_idx in 1:nrow(allrna_corr_matrix)) {
  for (c_idx in 1:ncol(allrna_corr_matrix)) {
    corr_val <- allrna_corr_matrix[r_idx, c_idx]
    p_val_adj <- allrna_p_values_adjusted[r_idx, c_idx]
    
    if (!is.na(corr_val)) { # This is the key condition
      text_val <- sprintf("%.1f", corr_val) 
      
      if (!is.na(p_val_adj) && p_val_adj >= corrected_alpha) {
        cellnote_matrix[r_idx, c_idx] <- " X"
      } else {
        cellnote_matrix[r_idx, c_idx] <- paste0("",text_val)
      }
    } else {
      # If corr_val is NA, cellnote_matrix entry remains NA_character_ (its initialized value)
      # cellnote_matrix[r_idx, c_idx] <- NA_character_ # Explicitly, though already NA
    }
  }
}

# Create the heatmap with dendrograms using heatmaply
heatmaply(allrna_corr_matrix, 
          Rowv = row_dendro, 
          Colv = col_dendro, 
          #col_side_colors = col_annotation_matrix,
          cellnote = cellnote_matrix, 
          cellnote_size = 8, 
          symm = TRUE,
          colors = colorRampPalette(c("navy", "white", "red"))(50), 
          limits = c(-1, 1), # MODIFIED: Explicitly set color scale limits
          grid_color = "grey80", 
          cellnote_textposition = "middle center", 
          grid_width = 0.001,    
          fontsize_row = 8,
          fontsize_col = 8,
          margins = c(60, 60, 40, 20) 
)


#####################################
#####################################
# VERSION 4, Add ks stats at the bottom 
# 3. Create the TWO individual plot components

# Plot A: Main Heatmap (Center)
main_heatmap <- heatmaply(
  allrna_corr_matrix,
  Rowv = row_dendro, Colv = col_dendro,
  show_dendrogram = c(FALSE, TRUE), 
  hide_colorbar = TRUE,
  
  xaxis_side = "top", 
  
  cellnote = cellnote_matrix,
  cellnote_size = 12, cellnote_textposition = "middle center",
  
  symm = TRUE,
  colors = colorRampPalette(c("navy", "white", "red"))(50),
  limits = c(-1, 1),
  grid_color = "grey80", grid_width = 0.001,
  
  fontsize_row = 10, fontsize_col = 10,
  margins = c(0, 0, 0, 0)
)


# Plot B: K-S D value annotation heatmap (Bottom)
ks_d_matrix <- data.frame(
  "mRNA KS" = ks_d_values_m$signed_D,
  "sncRNA KS" = ks_d_values_s$signed_D,
  "lncRNA KS" = ks_d_values_l$signed_D,
  "mean|KS|" = ks_results$meanAbsD,
  check.names = FALSE)

ks_d_matrix <- t(ks_d_matrix)
colnames(ks_d_matrix) <- SELECTED_INTRINSIC_LABELS

# Create cellnotes for the new multi-row matrix
ks_cellnote <- matrix(sprintf("%.2f", ks_d_matrix), nrow = nrow(ks_d_matrix), dimnames = dimnames(ks_d_matrix))

# We must reorder the columns to match the main heatmap's clustering
ks_d_matrix_ordered <- ks_d_matrix[, col_order, drop = FALSE]
ks_cellnote_ordered <- ks_cellnote[, col_order, drop = FALSE]

ks_heatmap <- heatmaply(
  ks_d_matrix_ordered,
  Rowv = FALSE, Colv = FALSE,
  cellnote = ks_cellnote_ordered,
  cellnote_color = "black", cellnote_size = 12, cellnote_textposition = "middle center",
  
  colors = colorRampPalette(c("navy", "white", "red"))(50),
  limits = c(-1, 1), # Use the same limits as the main heatmap for consistent color scale
  grid_color = "grey80", grid_width = 0.001,
  
  dendrogram = "none",
  # No longer using a single ylab, the row names of the matrix will be used
  
  fontsize_row = 10, fontsize_col = 10,
  margins = c(0, 0, 0, 0) # Margin on bottom for column labels
)


# 4. Combine the plots into a final layout
# Stack the two plots.
right_stack <- subplot(main_heatmap, ks_heatmap, 
                       nrows = 2, 
                       shareX = TRUE, 
                       margin = 0.005,
                       heights = c(0.9, 0.1), # Dendrogram, Main, Annotation
                       titleY = TRUE)


# Display the final plot
right_stack
