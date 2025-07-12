file <- "lncrna_epigenetic_features_matrix.csv"
file <- "short_ncrna_epigenetic_features_matrix.csv"
file <- "protein_epigenetic_features_matrix.csv"

data <- read.csv(file, header = TRUE, sep = "\t")
data$Functionality <- factor(data$Functionality, levels=c("No","Yes"))
functional_numeric <- as.numeric(factor(data$Functionality)) - 1
data$Functionality <- functional_numeric
data <- na.roughfix(data)
#summary(data)
cor_matrix <- cor(data, method = "spearman")

# Function to compute p-value matrix
calculate_pval_matrix <- function(data, method = "spearman") {
  n <- ncol(data)
  pval_matrix <- matrix(NA, nrow = n, ncol = n)
  diag(pval_matrix) <- 1 # Diagonals will be 1 (correlation of a variable with itself)
  
  col_names <- colnames(data) # Store column names
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      result <- cor.test(data[, i], data[, j], method = method)
      pval_matrix[i, j] <- pval_matrix[j, i] <- result$p.value
    }
  }
  
  # Set column and row names
  colnames(pval_matrix) <- col_names
  rownames(pval_matrix) <- col_names
  
  pval_matrix
}

# Calculate the p-value matrix (using Spearman correlation here)
pval_matrix <- calculate_pval_matrix(data)
#print(pval_matrix)


# Number of features
n <- ncol(cor_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Extract p-values
pvals <- pval_matrix[lower.tri(pval_matrix, diag = TRUE)]

# Adjust p-values using Bonferroni
p_values_adjusted <- p.adjust(pvals, method="bonferroni")

# Find significant correlations
#significant <- pvals < corrected_alpha
#significant

cor_matrix[lower.tri(cor_matrix)] <- NA
pval_matrix[lower.tri(pval_matrix)] <- NA

# Heatmap
library(ggplot2)
library(reshape2)
# Melt the correlation matrix into long format
# Reshape to suit ggplot, remove NAs, and sort the labels
cor_melted <- na.omit(melt(cor_matrix))
pval_melted <- na.omit(melt(pval_matrix))

combined_matrix <- merge(cor_melted, pval_melted, by = c("Var1", "Var2")) 

heatmap(cor_matrix)

# Inspiration from: http://pseudofish.com/triangle-heatmaps-in-r-using-ggplot.html
# Triangle heatmap
ggplot(combined_matrix, aes(Var1, Var2)) +
  ggtitle('Features Spearman Correlation - protein coding') +
  theme_bw() +
  geom_tile(aes(fill = value.x), color='white') +
  geom_text(aes(Var2, Var1, label = ifelse(value.y < corrected_alpha, 
                                           round(value.x, 2), "X"), 
                                           color = value.x), 
            nudge_x = 0, nudge_y = 0, fontface = "bold") + 
  scale_color_gradient2(low = "dodgerblue", mid = "bisque", high = "darkorchid4", midpoint = 0, limits = c(-1, 1)) +
  scale_fill_gradient2(low = "white", high = "white", midpoint = 0, guide = "none") +  # For background
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='white'),
        axis.text.y = element_text(hjust = 1) ) +
  labs(title = "Features Spearman Correlation - protein coding", x = "", y = "")


