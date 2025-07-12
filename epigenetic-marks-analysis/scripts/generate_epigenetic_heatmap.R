# Load necessary libraries
library(corrr)
library(Hmisc) # For rcorr() function
library(ggcorrplot)
library(dplyr)
library(reshape2)

library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
install.packages("ggplotify")
library(ggplotify) ## to convert pheatmap to ggplot2
install.packages("heatmaply")
library(heatmaply) ## for constructing interactive heatmap

# Define features of interest
SELECTED_HISTONE_FEATURES <- c("H2AFZ_MaxScaledSignal", "H2AK5ac_MaxScaledSignal", "H2AK9ac_MaxScaledSignal",
                               "H2BK120ac_MaxScaledSignal", "H2BK12ac_MaxScaledSignal", "H2BK15ac_MaxScaledSignal",
                               "H2BK20ac_MaxScaledSignal", "H2BK5ac_MaxScaledSignal", "H3F3A_MaxScaledSignal", 
                               "H3K14ac_MaxScaledSignal",
                               "H3K18ac_MaxScaledSignal", "H3K23ac_MaxScaledSignal", "H3K23me2_MaxScaledSignal",
                               "H3K27ac_MaxScaledSignal", "H3K27me3_MaxScaledSignal", "H3K36me3_MaxScaledSignal",
                               "H3K4ac_MaxScaledSignal", "H3K4me1_MaxScaledSignal", "H3K4me2_MaxScaledSignal",
                               "H3K4me3_MaxScaledSignal", "H3K56ac_MaxScaledSignal", "H3K79me1_MaxScaledSignal",
                               "H3K79me2_MaxScaledSignal", "H3K9ac_MaxScaledSignal", "H3K9me1_MaxScaledSignal",
                               "H3K9me2_MaxScaledSignal", "H3K9me3_MaxScaledSignal", #"H3T11ph_MaxScaledSignal",
                               "H4K12ac_MaxScaledSignal", "H4K20me1_MaxScaledSignal", "H4K5ac_MaxScaledSignal",
                               "H4K8ac_MaxScaledSignal", "H4K91ac_MaxScaledSignal", "chrm_acc_MaxScaledSignal",
                               "methylome")

SELECTED_HISTONE_LABELS <- c("H2AFZ", "H2AK5ac", "H2AK9ac",
                             "H2BK120ac", "H2BK12ac", "H2BK15ac",
                             "H2BK20ac", "H2BK5ac", "H3F3A",
                             "H3K14ac",
                             "H3K18ac", "H3K23ac", "H3K23me2",
                             "H3K27ac", "H3K27me3", "H3K36me3",
                             "H3K4ac", "H3K4me1", "H3K4me2",
                             "H3K4me3", "H3K56ac", "H3K79me1",
                             "H3K79me2", "H3K9ac", "H3K9me1",
                             "H3K9me2", "H3K9me3", #"H3T11ph",
                             "H4K12ac", "H4K20me1", "H4K5ac",
                             "H4K8ac", "H4K91ac", "chrm_acc",
                             "methylome")

# Load data
feature_matrix_histones <- read.csv("...", header = TRUE)
# Convert to numeric
data_numeric_histones <- feature_matrix_histones %>% 
  #dplyr::select(-Dataset) %>%
  sapply(function(feature) as.numeric(as.character(feature))) %>%
  as.data.frame()
# Restore Dataset column
data_numeric_histones$Dataset <- feature_matrix$Dataset

# Select data for analysis and run spearman correlation
allrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric_histones %>%
      filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3",
                            "lncrna-exon1", "lncrna-exon2",
                            "short-ncrna")) %>%
      select(all_of(SELECTED_HISTONE_FEATURES))
  ),
  type = "spearman"
)

# Unpack correlation coefficients and p-values #
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

# Adjust p-values using Benjamini-Hochberg
help("p.adjust")
allrna_p_values_adjusted <- p.adjust(allrna_pvalues, method="BH")

allrna_p_values_adjusted <- matrix(allrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(allrna_p_values_adjusted) <- diag

# Use short names
rownames(allrna_corr_matrix) <- SELECTED_HISTONE_LABELS
colnames(allrna_corr_matrix) <- SELECTED_HISTONE_LABELS

# Use short names
rownames(allrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS
colnames(allrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS

# Plot the heatmap
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
heatmap(allrna_corr_matrix)

pheatmap(allrna_corr_matrix, scale="row",
         
         color=colorRampPalette(c("navy", "white", "red"))(50),
         #cutree_cols=2, cutree_rows=2, 
         main="Spearman correlation Heatmap - All RNA types")


# Calculate hierarchical clustering for rows and columns
help("dist")
hc_rows <- hclust(dist(allrna_corr_matrix))
hc_cols <- hclust(dist(t(allrna_corr_matrix)))

# Create dendrograms
row_dendro <- as.dendrogram(hc_rows)
col_dendro <- as.dendrogram(hc_cols)

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

# --- DIAGNOSTICS FOR cellnote_matrix ---
message("--- DIAGNOSTICS: cellnote_matrix ---")
message("Head of cellnote_matrix:")
print(head(cellnote_matrix, 5))
num_actual_strings <- sum(!is.na(as.vector(cellnote_matrix)))
message(paste("Number of non-NA_character_ entries in cellnote_matrix:", num_actual_strings))
message(paste("Number of NA_character_ entries in cellnote_matrix:", sum(is.na(as.vector(cellnote_matrix)))))
if (num_actual_strings == 0 && nrow(allrna_corr_matrix) > 0) {
  warning("cellnote_matrix appears to be all NAs! This means allrna_corr_matrix might be all NAs or the loop logic isn't assigning strings.")
}
message("------------------------------------")
# --- END DIAGNOSTICS ---


# Create the heatmap with dendrograms using heatmaply
help("heatmaply")
heatmaply(allrna_corr_matrix, Rowv = row_dendro, Colv = col_dendro, 
          cellnote = cellnote_matrix, cellnote_size = 8, symm = TRUE,
          colors = colorRampPalette(c("navy", "white", "red"))(50), 
          limits = c(-1, 1), # MODIFIED: Explicitly set color scale limits
          grid_color = "grey80", cellnote_textposition = "middle center", 
          grid_width = 0.001,    
          fontsize_row = 8,
          fontsize_col = 8,
          margins = c(60, 60, 40, 20) 
          )
