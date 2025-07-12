# install.packages("ggpmisc") # Uncomment this line if you need to install ggpmisc

# Load libraries
library(dplyr)
library(ggplot2)
library(ggpmisc) # For adding regression equation and R-squared

# Load positive csv files
mrna_bigwig_data <- read.csv("../data/datasets/histone_feature/bigWig/H3K79me2/H3K79me2_protein-exon2-histone-feature.csv", header = TRUE)
mrna_narrowpeak_data <- read.csv("../data/datasets/histone_feature/H3K79me2/H3K79me2_protein-exon2-histone-feature.csv", header = TRUE)

# Load negative csv files
mrna_bigwig_nc_data <- read.csv("../data/datasets/histone_feature/bigWig/H3K79me2/H3K79me2_protein-exon2-NC-histone-feature.csv", header = TRUE)
mrna_narrowpeak_nc_data <- read.csv("../data/datasets/histone_feature/H3K79me2/H3K79me2_protein-exon2-NC-histone-feature.csv", header = TRUE)

# Join datasets
mrna_signal_data <- data.frame("bigWig_signal" = mrna_bigwig_data$H3K79me2_AvgSignal, 
                               "narrowPeak_signal" = mrna_narrowpeak_data$H3K79me2_AvgSignal)
mrna_nc_signal_data <- data.frame("bigWig_signal" = mrna_bigwig_nc_data$H3K79me2_AvgSignal, 
                                  "narrowPeak_signal" = mrna_narrowpeak_nc_data$H3K79me2_AvgSignal)

# Add a 'source' column to each data frame to identify them
mrna_signal_data$source <- "mRNA (+)"
mrna_nc_signal_data$source <- "mRNA (-)"

# Combine the two data frames into one
combined_data <- rbind(mrna_signal_data, mrna_nc_signal_data)


# Define the formula for the regression model
formula <- y ~ x

# Create the plot using the combined data
# Map the 'source' column to the color aesthetic
ggplot(combined_data, aes(y = narrowPeak_signal, x = bigWig_signal)) + # Base aesthetics
  geom_point(aes(color = source), alpha = 0.7, size = 2) + # Color points by source
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = formula) + # Add overall linear regression line (black)
  # Updated stat_poly_eq call to address warnings:
  stat_poly_eq(formula = formula,
               # Use after_stat() instead of ..notation..
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               parse = TRUE, # Parse the label for mathematical notation
               # Use label.x/label.y with hjust/vjust for positioning
               label.x = Inf, # Position horizontally (Inf = right edge)
               label.y = -Inf, # Position vertically (-Inf = bottom edge)
               hjust = 1.05,  # Horizontal justification ( > 1 moves left from edge)
               vjust = -6.2, # Vertical justification ( < 0 moves down from edge)
               size = 10) + # Adjust text size
  scale_color_manual(values = c("mRNA (+)" = "blue", "mRNA (-)" = "red")) + # Define specific colors for points
  labs(
    title = "Comparison of H3K79me2 Signal Data",
    x = "bigWig Signal",
    y = "narrowPeak Signal",
    color = "Gene Type" # Customize legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    panel.grid.major = element_line(color = "gray90"),  # Lighter grid lines
    panel.grid.minor = element_line(color = "gray95")   # Lighter grid lines
  )

