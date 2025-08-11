# Script to analyze Montero et.al 2024 essentiality groups against our model predictions.
# Load libraries
library(dplyr)
library(stringr) # For str_detect
library(ggplot2)
library(ggsignif)

# Load csv file
essentials_exon1_data <- rbind(read.csv("../data/Montero/processed/annotated_table9_grna1_vs_exon1_prob.csv", header = TRUE),
                               read.csv("../data/Montero/processed/annotated_table9_grna2_vs_exon1_prob.csv", header = TRUE))
essentials_exon2_data <- rbind(read.csv("../data/Montero/processed/annotated_table9_grna1_vs_exon2_prob.csv", header = TRUE),
                               read.csv("../data/Montero/processed/annotated_table9_grna2_vs_exon2_prob.csv", header = TRUE))

ess_e1_d <- essentials_exon1_data %>% select(Q_name_suffix,GeneID,TranscriptID,Probability_Functional,Essential_Status)
ess_e2_d <- essentials_exon2_data %>% select(Q_name_suffix,GeneID,TranscriptID,Probability_Functional,Essential_Status)

# Join the two datasets by Q_name and GeneID
ess_d_join <- full_join(ess_e1_d, ess_e2_d, by = c("Q_name_suffix", "GeneID", "Essential_Status"))

# Get the max probability.
ess_d_join <- ess_d_join %>% mutate(max_prob = pmax(Probability_Functional.x, Probability_Functional.y, na.rm = TRUE))

# Define the desired order for the categorical variable
essentiality_labels <- c("Non-essential", "Essential")

# Reorder the 'Essential_Status' column by converting it to a factor with specified levels
ess_d_join <- ess_d_join %>%
  mutate(Essential_Status = factor(Essential_Status, levels = c("NOT ESSENTIAL","ESSENTIAL"), labels = essentiality_labels))

# Keep only the max probability per GeneID
filtered_df <- ess_d_join %>%
  group_by(GeneID) %>%
  # For each GeneID, keep only the row with the maximum Probability_Functional.
  # `with_ties = FALSE` ensures that if there's a tie, only one row is returned.
  slice_max(order_by = max_prob, n = 1, with_ties = FALSE) %>%
  ungroup() # It's good practice to ungroup after the operation.


# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filtered_df %>%
  count(Essential_Status)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$Essential_Status, " (", legend_data$n, ")"),
  legend_data$Essential_Status
)



# Generate the Violin Plot
plot_violin <- ggplot(data = filtered_df,
                      aes(x = Essential_Status, y = max_prob, fill = Essential_Status)) +
  
  # Add the violin layer
  # trim = FALSE ensures the violin tails are drawn to the full range of the data
  geom_violin(scale = "width", na.rm = TRUE, adjust = 1, trim = FALSE) +
  
  
  # Add a box plot inside the violins
  geom_boxplot(alpha=0.0, outliers=TRUE, na.rm = TRUE, position = position_dodge(width = 0.9), width=0.2) +
  
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_manual(values = c("Essential" = "#1a53ff", "Non-essential" = "#beb9db"),
                    labels = new_labels) +
  
  # Update the labels for the new plot layout
  labs(
    title = "Distribution of Functional Probability by Essentiality Group, Montero et.al 2024",
    x = "Essentiality Group",
    y = "Functional Probability",
    fill = "Essentiality Group" # Update legend title to match 'fill'
  ) +
  
  theme_minimal() +
  theme(
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(), # Hiding minor grid lines for a cleaner look
    legend.position = "right"
  )


# Define the pairwise comparisons to be performed
my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)

# Custom function to perform K-S test and format the D-statistic and p-value stars
ks_test_custom <- function(x, y) {
  test <- ks.test(x, y)
  
  # Convert p-value to significance stars
  # 'cut' is a clean way to assign stars based on p-value ranges
  stars <- cut(test$p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
               labels = c("***", "**", "*", "ns"))
  
  # Format the label as: D-value (stars)
  label <- paste0(format(test$statistic, digits = 2), " (", stars, ")")
  
  # Return a list with the custom label for ggsignif
  return(list(p.value = label))
}

# Add the statistical comparison layer using the custom function
plot_with_stats <- plot_violin +
  geom_signif(
    comparisons = my_comparisons,
    test = "ks_test_custom",  # Use our new custom function
    step_increase = 0.1,      # Spacing can be reduced for single-line labels
    textsize = 6.5,
    tip_length = 0.01,
    y_position = 1.5      # Manually set the y-position for the first bracket
  ) +
  # Expand the y-axis slightly to make room for the highest comparison brackets
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))


# To display the final plot
print(plot_with_stats)

