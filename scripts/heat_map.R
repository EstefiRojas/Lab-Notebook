file <- "lncrna_epigenetic_features_matrix.csv"
file <- "short_ncrna_epigenetic_features_matrix.csv"
file <- "protein_epigenetic_features_matrix.csv"

data <- read.csv(file, header = TRUE, sep = "\t")
data$Functionality <- factor(data$Functionality, levels=c("No","Yes"))
functional_numeric <- as.numeric(factor(data$Functionality)) - 1
data$Functionality <- functional_numeric
data <- na.roughfix(data)
cor_matrix <- cor(data, method = "spearman")
cor_matrix[lower.tri(cor_matrix)] <- NA
# Heatmap
library(ggplot2)
library(reshape2)
# Melt the correlation matrix into long format
# Reshape to suit ggplot, remove NAs, and sort the labels
cor_melted <- na.omit(melt(cor_matrix))

#cor_melted <- melt(cor_matrix)

heatmap(cor_matrix)

# Triangle heatmap to compare cohorts
ggplot(cor_melted, aes(Var1, Var2)) +
  ggtitle('Features Spearman Correlation - protein coding') +
  theme_bw() +
  geom_tile(aes(fill = value), color='white') +
  geom_text(aes(Var2, Var1, label = round(value, 2), color = value), nudge_x = 0, nudge_y = 0, fontface = "bold") + 
  #scale_color_gradient2(low = 'dodgerblue', mid = "bisque", high = 'darkorchid4', space = 'Lab') +
  scale_color_gradient2(low = "dodgerblue", mid = "bisque", high = "darkorchid4", midpoint = 0, limits = c(-1, 1)) +
  scale_fill_gradient2(low = "white", high = "white", midpoint = 0, guide = "none") +  # For background
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='white'),
        axis.text.y = element_text(hjust = 1) ) +
  labs(title = "Features Spearman Correlation - protein coding", x = "", y = "")
