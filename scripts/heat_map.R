# 0. Load necessary libraries
library(randomForest)
library(ggplot2)
library(reshape2)
library('corrr')
library(ggcorrplot)

# 1. Load data, choose one of the following:
file <- "lncrna_epigenetic_features_matrix.csv"
file <- "short_ncrna_epigenetic_features_matrix.csv"
file <- "protein_epigenetic_features_matrix.csv"

data <- read.csv(file, header = TRUE, sep = "\t")

# 2. Convert Functional column to numeric
data$Functionality <- factor(data$Functionality, levels=c("No","Yes"))
functional_numeric <- as.numeric(factor(data$Functionality)) - 1
data$Functionality <- functional_numeric

# 3. Remove NA entries
data <- na.roughfix(data)

# 4. Compute Spearman Correlation
cor_matrix <- cor(data, method = "spearman")

# 5. Create triangular shaped Heatmap
# Keep only upper values from correlation matrix
cor_matrix[lower.tri(cor_matrix)] <- NA

# Melt the correlation matrix into long format
# Reshape to suit ggplot, remove NAs, and sort the labels
cor_melted <- na.omit(melt(cor_matrix))

# Triangular heatmap without p-values
ggplot(cor_melted, aes(Var1, Var2)) +
  ggtitle('Features Spearman Correlation - lincRNA') +
  theme_bw() +
  geom_tile(aes(fill = value), 
            color='white') +
  geom_text(aes(Var2, Var1, l
                abel = round(value, 2), 
                color = value), 
            nudge_x = 0, 
            nudge_y = 0, 
            fontface = "bold") + 
  #scale_color_gradient2(low = 'dodgerblue', mid = "bisque", high = 'darkorchid4', space = 'Lab') +
  scale_color_gradient2(low = "dodgerblue", 
                        mid = "bisque", 
                        high = "darkorchid4", 
                        midpoint = 0, 
                        limits = c(-1, 1)) +
  scale_fill_gradient2(low = "white", 
                       high = "white", 
                       midpoint = 0, 
                       guide = "none") +  # For background
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='white'),
        axis.text.y = element_text(hjust = 1) ) +
  labs(title = "Features Spearman Correlation - lincRNA", 
       x = "", 
       y = "")


# 6. Compute p-values and Bonferroni correction
# Compute correlation using rcorr library
corr_obj <- rcorr(as.matrix(data), type = "spearman")

# Unpack correlation and p-pvalues
corr_matrix <- corr_obj[["r"]]
pvalues <- corr_obj[["P"]]

# Number of features
n <- ncol(data_normalized_sncrna)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
p_values_adjusted <- p.adjust(pvalues, method="bonferroni")
p_values_adjusted <- matrix(p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(p_values_adjusted) <- diag

# Triangular heatmap with significance marking
ggcorrplot(corr_matrix, 
           type = "lower", 
           method = "square", 
           lab = TRUE, 
           lab_col = "black", 
           lab_size = 3, 
           ggtheme = theme_void, 
           title = "Spearman correlation Heatmap - lincRNA", 
           show.diag = TRUE, 
           p.mat = p_values_adjusted, 
           sig.level = 0.05, 
           insig = "pch",
           ) +
  theme(
    plot.title = element_text(size = 24),  # Increase title size
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    
  )
ggsave("spearmanCorrlincRNA.png",path = "../results/SpearmanCorAnalysis/", scale = 3, width = 1920, height = 1080, units = "px", bg = "white")

