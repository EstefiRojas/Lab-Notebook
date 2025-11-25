# Load libraries
library(dplyr)
library(stringr) # For str_detect
library(ggplot2)
library(ggsignif)

# Load csv file
essentials1_lncrna_data <- read.csv("../results/unified_genome_alignments.csv", header = TRUE)
table(essentials1_lncrna_data$Study)

ess_d <- essentials1_lncrna_data %>% filter(Protein_Off_Target == "NO" & !is.na(Probability_Functional) & Study == "Liang")
table(ess_d$Essentiality)
# Define the desired order for the categorical variable
essentiality_values <- c("Non-essential", "Rare", "Common", "Core")
essentiality_labels <- c("Non", "Rare", "Common", "Core")

# Reorder the 'Essentiality' column by converting it to a factor with specified levels
ess_d <- ess_d %>%
  mutate(Essentiality = factor(Essentiality, levels = essentiality_values, labels = essentiality_labels))


# Keep only the max probability per GeneID
filtered_df <- ess_d %>%
  group_by(ENSG_ID, Essentiality) %>%
  # For each GeneID, keep only the row with the maximum Probability_Functional.
  # `with_ties = FALSE` ensures that if there's a tie, only one row is returned.
  slice_max(order_by = Probability_Functional, n = 1, with_ties = FALSE) %>%
  ungroup() # It's good practice to ungroup after the operation.


table(filtered_df$Essentiality)

# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filtered_df %>%
  count(Essentiality)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$Essentiality, "\n(n=", legend_data$n, ")"),
  legend_data$Essentiality
)


# --- Generate the Violin Plot ---

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("Non","Rare"),
                       c("Non","Common"),
                       c("Non","Core"))

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


# Compute KS stats

# --- Step 1: Pre-calculate statistics for labels ---

# Define your reference group and comparison groups
reference_group <- "Non"
comparison_groups <- c("Rare", "Common", "Core")

# Extract the data for the reference group
reference_data <- filtered_df %>%
  filter(Essentiality == reference_group) %>%
  pull(Probability_Functional)

# Calculate KS statistic for each comparison group
stats_list <- lapply(comparison_groups, function(group) {
  # Extract data for the current comparison group
  comparison_data <- filtered_df %>%
    filter(Essentiality == group) %>%
    pull(Probability_Functional)
  
  # Perform the KS test
  ks_result <- ks.test(reference_data, comparison_data)
  
  # Return a data frame with the necessary info for plotting
  data.frame(
    Essentiality = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

# Combine the list of data frames into a single data frame
stats_labels <- do.call(rbind, stats_list)

# Define the y-position for the labels (adjust as needed)
stats_labels$y_position <- 1.1 


#help("geom_boxplot")
# 14 Jul 2025 #
# --- Start of the plot code ---
plot_modified <- ggplot(data = filtered_df,
                        aes(x = Essentiality, y = Probability_Functional, 
                            fill = Essentiality, color = Essentiality)) +
  
  # 1. Add box plot
  geom_boxplot(data = ~ subset(., Essentiality != "Core"), linewidth = 0.9,
               na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  
  # 2. Add a jitter plot for the 'Shared' group
  geom_jitter(data = ~ subset(., Essentiality == "Core"),
             size = 6, 
             shape = 21, 
             color = "black",
             stroke = 1) +
  
  # 3. Manually set fill colors
  scale_fill_brewer(palette = "Set2") +
  
  # 4. Manually set border colors
  scale_colour_brewer(palette = "Set2") +
  
  scale_x_discrete(labels = new_labels) +
  
  # 5. Override the legend glyph to be a filled circle
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5, color = "black"))) +
  
  # 6. Update labels (No changes here)
  labs(
    title = "Gene Essentiality",
    y = "lncRNA Probability",
    fill = "Essentiality Group"
  ) +
  
  # 7. Apply custom theme (No changes here)
  theme_minimal() +
  theme(
    text = element_text(size = 36),
    plot.title = element_text(size = 46, hjust = 0.5),
    axis.text.x = element_text(hjust = 0.5, size = 30),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 44),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3)
  )

# Add the statistical comparison layer
plot_with_text_stats <- plot_modified +
  
  # Add the pre-calculated stats as text
  geom_text(
    data = stats_labels,
    aes(x = Essentiality, y = y_position, label = label),
    inherit.aes = FALSE, 
    size = 9.5,
    color = "black",
    fontface = "bold"
  ) +
  
  # Adjust y-axis limits to ensure text is visible
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.15), clip = "off") # Increased ylim slightly

# Display the final plot
print(plot_with_text_stats)

