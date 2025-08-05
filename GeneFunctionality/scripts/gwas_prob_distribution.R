# Load libraries
library(dplyr)
library(tidyr)
library(stringr) # For str_detect
library(ggplot2)
library(ggsignif)

# Load csv file
gwas_data <- read.csv("../data/gwas/lncrna-gwas-ranked-5-aug-2025.csv", header = TRUE)

filter_data <- gwas_data %>%
  separate(Probability_Functional, 
           into = c("tl_prob", "tr_prob"), 
           sep = "\\|",
           remove = FALSE) %>% 
  mutate(tl_prob = as.numeric(tl_prob),
         tr_prob = as.numeric(tr_prob),
         tl_abs_beta_for_min_p = abs(tl_beta_for_min_p),
         tr_abs_beta_for_min_p = abs(tr_beta_for_min_p),
         # First determine which beta to use based on probability
         selected_abs_beta = ifelse(tl_prob >= tr_prob, tl_abs_beta_for_min_p, tr_abs_beta_for_min_p),
         # Then classify based on the selected beta
         abs_beta_ratio_category = case_when(
           selected_abs_beta < 0.025 ~ "<0.025",
           selected_abs_beta >= 0.025 & selected_abs_beta < 0.05 ~ "0.025-0.05",
           selected_abs_beta >= 0.05 & selected_abs_beta < 0.5 ~ "0.05-0.5",
           selected_abs_beta >= 0.5 ~ "≥0.5",
           TRUE ~ "NA"
         )
  )

filter_data$abs_beta_ratio_category <- factor(filter_data$abs_beta_ratio_category,
                                              levels = c("NA","<0.025","0.025-0.05","0.05-0.5","≥0.5"))

# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filter_data %>%
  count(abs_beta_ratio_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$abs_beta_ratio_category, " (", legend_data$n, ")"),
  legend_data$abs_beta_ratio_category
)
help("scale_x_log10")
ggplot(filter_data, aes(x = selected_abs_beta)) +
  geom_histogram() +
  scale_x_log10(breaks=c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100)) +
  labs(
    title = "GWAS Beta",
    x = "Absolute Beta"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )


# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("NA","<0.025"),
                       c("NA","0.025-0.05"),
                       c("NA","0.05-0.5"),
                       c("NA","≥0.5"))

# Custom function to perform K-S test and format the D-statistic and p-value stars
ks_test_custom <- function(x, y) {
  test <- ks.test(x, y)
  
  # Convert p-value to significance stars
  # 'cut' is a clean way to assign stars based on p-value ranges
  stars <- cut(test$p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
               labels = c("***", "**", "*", "ns"))
  
  # Format the label as: D-value (stars)
  label <- paste0(format(test$statistic, digits = 2), " (", format(test$p.value, digits = 2), ")")
  
  # Return a list with the custom label for ggsignif
  return(list(p.value = label))
}


plot_modified <- ggplot(filter_data, aes(x = abs_beta_ratio_category, y = highest_prob, fill = abs_beta_ratio_category)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.3, outliers=TRUE, na.rm = TRUE, position = position_dodge(width = 0.9), width=0.2) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS Beta",
    x = "Absolute Beta",
    y = "lncRNA Probability"
  ) +
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(), # Hiding minor grid lines for a cleaner look
    legend.position = "none"
  )
  
# Add the statistical comparison layer
plot_with_stats <- plot_modified +
  geom_signif(
    comparisons = my_comparisons,
    test = "ks_test_custom",
    step_increase = 0.2,
    textsize = 9.5,
    tip_length = 0.01,
    y_position = 1.2
  ) +
  # KEY CHANGE: Set explicit breaks for the y-axis and use coord_cartesian to set the visual range.
  # This prevents nonsensical axis ticks (e.g., > 1) while keeping room for annotations.
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.7), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations


# Display the final plot
print(plot_with_stats)
  
  
  
