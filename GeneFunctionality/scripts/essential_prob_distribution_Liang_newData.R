# Load libraries
library(dplyr)
library(stringr) # For str_detect
library(ggplot2)
library(ggsignif)

# Load csv file
essentials_exon1_data <- read.csv("../data/Liang/processed/annotated_tableS3_gRNAs_vs_exon1_prob.csv", header = TRUE)
essentials_exon2_data <- read.csv("../data/Liang/processed/annotated_tableS3_gRNAs_vs_exon2_prob.csv", header = TRUE)


ess_e1_d <- essentials_exon1_data %>% select(Q_name_suffix,GeneID,TranscriptID,Probability_Functional,Essential_Status)
ess_e2_d <- essentials_exon2_data %>% select(Q_name_suffix,GeneID,TranscriptID,Probability_Functional,Essential_Status)

# Join the two datasets by Q_name and GeneID
ess_d_join <- full_join(ess_e1_d, ess_e2_d, by = c("Q_name_suffix", "GeneID", "Essential_Status"))

# Get the max probability.
ess_d_join <- ess_d_join %>% mutate(max_prob = pmax(Probability_Functional.x, Probability_Functional.y, na.rm = TRUE))


table(ess_d_join$Essential_Status)

# Define the desired order for the categorical variable
essentiality_values <- c("Non-essential", "Cell-type specific", "Partially shared", "Shared")
essentiality_labels <- c("Non", "Specific", "Partial", "Shared")

# Reorder the 'Essential_Status' column by converting it to a factor with specified levels
ess_d_join <- ess_d_join %>%
  mutate(Essential_Status = factor(Essential_Status, levels = essentiality_values, labels = essentiality_labels))

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
  paste0(legend_data$Essential_Status, "\n(n=", legend_data$n, ")"),
  legend_data$Essential_Status
)


# --- Generate the Violin Plot ---

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("Non","Specific"),
                       c("Non","Partial"),
                       c("Non","Shared"))

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
comparison_groups <- c("Specific", "Partial", "Shared")

# Extract the data for the reference group
reference_data <- filtered_df %>%
  filter(Essential_Status == reference_group) %>%
  pull(max_prob)

# Calculate KS statistic for each comparison group
stats_list <- lapply(comparison_groups, function(group) {
  # Extract data for the current comparison group
  comparison_data <- filtered_df %>%
    filter(Essential_Status == group) %>%
    pull(max_prob)
  
  # Perform the KS test
  ks_result <- ks.test(reference_data, comparison_data)
  
  # Return a data frame with the necessary info for plotting
  data.frame(
    Essential_Status = group,
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
                        aes(x = Essential_Status, y = max_prob, 
                            fill = Essential_Status, color = Essential_Status)) +
  
  # 1. Add box plot
  geom_boxplot(data = ~ subset(., Essential_Status != "Shared"), linewidth = 0.9,
               na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  
  # 2. Add a jitter plot for the 'Shared' group
  geom_point(data = ~ subset(., Essential_Status == "Shared"),
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
    aes(x = Essential_Status, y = y_position, label = label),
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



# Generate a Jitter Plot
plot_jitter <- ggplot(data = filtered_df,
                      aes(x = Essential_Status, y = max_prob, color = Essential_Status)) +
  
  # Add the jitter layer
  # 'width' controls the horizontal spread of points
  # 'alpha' makes points semi-transparent to show density
  geom_jitter(na.rm = TRUE, width = 0.25, alpha = 0.6, size = 2) +
  
  # Use scale_color_manual since we mapped the 'color' aesthetic
  scale_color_manual(values = c("Cell-type specific" = "#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db"),
                     labels = new_labels) +
  
  # Update the labels for the new plot layout
  labs(
    title = "Distribution of Functional Probability by Essentiality Group, Liang et.al 2024",
    x = "Essentiality Group",
    y = "Functional Probability",
    color = "Essentiality Group" # Update legend title to match 'color'
  ) +
  
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    text = element_text(size = 20), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.title = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold"),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# Add the statistical comparison layer using the custom function
plot_with_stats <- plot_jitter +
  geom_signif(
    comparisons = my_comparisons,
    test = "ks_test_custom",
    step_increase = 0.2,
    textsize = 6.5,
    tip_length = 0.01,
    y_position = 1.2
  ) +
  # KEY CHANGE: Set explicit breaks for the y-axis and use coord_cartesian to set the visual range.
  # This prevents nonsensical axis ticks (e.g., > 1) while keeping room for annotations.
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.7), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations


# To display the final plot
print(plot_with_stats)


plot_scatter <- ggplot(data = ess_d_join,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -selected, y = Probability_Functional, color = group_day14)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Optional: Apply a specific color scale (e.g., ColorBrewer Set1)
  # You might want scale_color_manual() to assign specific colors to "Essential" and "Non-essential"
  # Add smooth trend lines (linear model 'lm' for each group)
  # se = TRUE adds the confidence interval ribbon (default)
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) + # Added trend line layer
  scale_color_manual(values = c("Cell-type specific" ="#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Functional Probability (exon1) vs gRNA Depletion Liang 2024", # Your title
    x = "Sum depletion value day 14 and day7 (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )

print(plot_scatter)






##################################
# OLD LIANG DATA
essentials_data_w <- essentials_data_w %>%
  mutate(
    day7_sum_all = rowSums(
      select(., contains(".Day7.")),
      na.rm = TRUE
    ),
    day14_sum_all = rowSums(
      select(., contains(".Day14.")),
      na.rm = TRUE
    ),
    sum_all = rowSums(
      select(., contains(".Day")),
      na.rm = TRUE
    ),
    delta14_7 = day14_sum_all - day7_sum_all
  )

essentials_exon2_data_w <- essentials_exon2_data_w %>%
  mutate(
    day7_sum_all = rowSums(
      select(., contains(".Day7.")),
      na.rm = TRUE
    ),
    day14_sum_all = rowSums(
      select(., contains(".Day14.")),
      na.rm = TRUE
    ),
    sum_all = rowSums(
      select(., contains(".Day")),
      na.rm = TRUE
    ),
    delta14_7 = day14_sum_all - day7_sum_all
  )

# --- Format group_day14 into a new Factor ---
# Create a new column 'Essentiality_Group' based on 'group_day14'
essentials_data_w <- essentials_data_w %>%
  mutate(
    Essentiality_Group = case_when(
      # Map the first three groups to "Essential"
      group_day14 %in% c("Cell-type specific", "Partially shared", "Shared") ~ "Essential",
      # Map "Non-essential" to itself
      group_day14 == "Non-essential" ~ "Non-essential",
      # Optional: handle any unexpected values (assign NA or a default category)
      TRUE ~ NA_character_
    ),
    # Convert the new column to a factor with specified levels
    # This ensures the desired order in legends/plots
    Essentiality_Group = factor(Essentiality_Group, levels = c("Essential", "Non-essential"))
  )

essentials_exon2_data_w <- essentials_exon2_data_w %>%
  mutate(
    Essentiality_Group = case_when(
      # Map the first three groups to "Essential"
      group_day14 %in% c("Cell-type specific", "Partially shared", "Shared") ~ "Essential",
      # Map "Non-essential" to itself
      group_day14 == "Non-essential" ~ "Non-essential",
      # Optional: handle any unexpected values (assign NA or a default category)
      TRUE ~ NA_character_
    ),
    # Convert the new column to a factor with specified levels
    # This ensures the desired order in legends/plots
    Essentiality_Group = factor(Essentiality_Group, levels = c("Essential", "Non-essential"))
  )

# --- Plotting the Distribution ---

# 1. Histogram
# geom_histogram bins the data and shows the count in each bin.
# Adjust 'binwidth' or 'bins' to change the appearance.
plot_histogram <- ggplot(data = essentials_data_w, aes(x = p_results_min_p)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "white", alpha = 0.8) +
  labs(
    title = "Distribution of Sum Difference (Day7 - Day14)",
    x = "Sum Difference",
    y = "Frequency (Count)"
  ) +
  theme_minimal() # A clean theme

print(plot_histogram)


# 2. Density Plot
# geom_density shows a smoothed estimate of the distribution.
plot_density <- ggplot(data = essentials_data_w, aes(x = delta14_7)) +
  geom_density(fill = "lightcoral", alpha = 0.7) +
  labs(
    title = "Density Plot of Sum Difference (Day7 - Day14)",
    x = "Sum Difference",
    y = "Density"
  ) +
  theme_minimal()

print(plot_density)

# 3. Histogram with Density Overlay (Optional)
# Using aes(y = ..density..) tells the histogram to plot density instead of count.
plot_hist_density <- ggplot(data = essentials_data_w, aes(x = delta14_7)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 5, fill = "skyblue", color = "white", alpha = 0.7) + # Use density scale
  geom_density(color = "red", linewidth = 1) + # Overlay density line
  labs(
    title = "Distribution of Sum Difference with Density Overlay",
    x = "Sum Difference",
    y = "Density"
  ) +
  theme_minimal()

print(plot_hist_density)


# --- Plotting the Relationship ---

# 4. Scatter Plot (-log2 Min Difference vs. Functional Probability)
# geom_point creates a scatter plot to show potential correlation/relationship.
# Plotting -sum_difference on x-axis as requested, assuming sum_difference is already log2 scale.
plot_scatter <- ggplot(data = essentials_data_w, aes(x = - delta14_7, y = Probability_Functional)) +
  geom_point(alpha = 0.6, color = "darkgreen", size = 1.5) + # Add points with some transparency and size
  labs(
    title = "Functional Probability vs gRNA Depletion Liang 2024", # Updated title
    x = "delta(day14 - day7) (-log2 Scale)", # Updated x-axis label to reflect the transformation
    y = "Functional Probability"
  ) +
  theme_minimal() + 
  theme(
    text = element_text(size = 30),
    panel.grid.major = element_line(color = "gray90"),  # Lighter grid lines
    panel.grid.minor = element_line(color = "gray95")   # Lighter grid lines
  )

print(plot_scatter)


# 5. Scatter Plot (-log2 Sum Difference vs. Probability Functional)
# geom_point creates a scatter plot to show potential correlation/relationship.
# Plotting -sum_difference on x-axis as requested, assuming sum_difference is already log2 scale.
plot_scatter <- ggplot(data = essentials_data_w, aes(x = - sum_all, y = Probability_Functional)) + # Use -sum_difference for x-axis
  geom_point(alpha = 0.6, color = "darkgreen", size = 1.5) + # Add points with some transparency and size
  labs(
    title = "Functional Probability vs gRNA Depletion Liang 2024", # Updated title
    x = "sum(day14 + day7) (-log2 Scale)", # Updated x-axis label to reflect the transformation
    y = "Functional Probability"
  ) +
  theme_minimal() + 
  theme(
    text = element_text(size = 30),
    panel.grid.major = element_line(color = "gray90"),  # Lighter grid lines
    panel.grid.minor = element_line(color = "gray95")   # Lighter grid lines
  )

print(plot_scatter)



# --- Plotting the Relationship ---
# 4. Scatter Plot (-log2 Sum All vs. Probability Functional, colored by Essentiality_Group)
# geom_point creates a scatter plot. Color is mapped to the new Essentiality_Group factor.
plot_scatter <- ggplot(data = essentials_data_w,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -day14_min_calculated, y = Probability_Functional, color = group_day14)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Add smooth trend lines (linear model 'lm' for each group)
  # se = TRUE adds the confidence interval ribbon (default)
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) + # Added trend line layer
  # You might want scale_color_manual() to assign specific colors to "Essential" and "Non-essential"
  scale_color_manual(values = c("Cell-type specific" ="#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Functional Probability (exon 1) vs gRNA Depletion Liang 2024", # Your title
    x = "Minimum depletion value at day 14 (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )

print(plot_scatter)

#Sum
#essentials_data_w_filtered <- essentials_data_w %>%
#  filter(sum_all_calculated < 0)
# 4. Scatter Plot (-log2 Sum All vs. Probability Functional, colored by Essentiality_Group)
# geom_point creates a scatter plot. Color is mapped to the new Essentiality_Group factor.
plot_scatter <- ggplot(data = essentials_data_w,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -sum_all_calculated, y = Probability_Functional, color = group_day14)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Optional: Apply a specific color scale (e.g., ColorBrewer Set1)
  # You might want scale_color_manual() to assign specific colors to "Essential" and "Non-essential"
  # Add smooth trend lines (linear model 'lm' for each group)
  # se = TRUE adds the confidence interval ribbon (default)
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) + # Added trend line layer
  scale_color_manual(values = c("Cell-type specific" ="#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Functional Probability (exon1) vs gRNA Depletion Liang 2024", # Your title
    x = "Sum depletion value day 14 and day7 (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )

print(plot_scatter)

# 4. Scatter Plot (-log2 Sum All vs. Probability Functional, colored by Essentiality_Group)
# geom_point creates a scatter plot. Color is mapped to the new Essentiality_Group factor.
plot_scatter <- ggplot(data = essentials_exon2_data_w,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -day14_min_calculated, y = Probability_Functional, color = group_day14)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Add smooth trend lines (linear model 'lm' for each group)
  # se = TRUE adds the confidence interval ribbon (default)
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) + # Added trend line layer
  # Optional: Apply a specific color scale (e.g., ColorBrewer Set1)
  # You might want scale_color_manual() to assign specific colors to "Essential" and "Non-essential"
  scale_color_manual(values = c("Cell-type specific" ="#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Functional Probability (exon2) vs gRNA Depletion Liang 2024", # Your title
    x = "Minimum depletion value at day 14 (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )

print(plot_scatter)

#Sum
#essentials_data_w_filtered <- essentials_data_w %>%
#  filter(sum_all_calculated < 0)
# 4. Scatter Plot (-log2 Sum All vs. Probability Functional, colored by Essentiality_Group)
# geom_point creates a scatter plot. Color is mapped to the new Essentiality_Group factor.
plot_scatter <- ggplot(data = essentials_exon2_data_w,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -sum_all_calculated, y = Probability_Functional, color = group_day14)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Add smooth trend lines (linear model 'lm' for each group)
  # se = TRUE adds the confidence interval ribbon (default)
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) + # Added trend line layer
  # Optional: Apply a specific color scale (e.g., ColorBrewer Set1)
  # You might want scale_color_manual() to assign specific colors to "Essential" and "Non-essential"
  scale_color_manual(values = c("Cell-type specific" ="#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Functional Probability (exon2) vs gRNA Depletion Liang 2024", # Your title
    x = "Sum depletion value day 14 and day7 (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )

print(plot_scatter)



######################
# Join both datasets #
essentials_data_all <- rbind(essentials_data_w, essentials_exon2_data_w)

# Plot -log10(min p_value)
plot_scatter <- ggplot(data = essentials_data_all,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -log10(p_results_min_p), y = Probability_Functional, color = group_day14)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Add smooth trend lines (linear model 'lm' for each group)
  # se = TRUE adds the confidence interval ribbon (default)
  geom_smooth(method = "loess", span = 2, se = FALSE, linewidth = 0.8, alpha = 0.15) + # Added trend line layer
  # Optional: Apply a specific color scale (e.g., ColorBrewer Set1)
  # You might want scale_color_manual() to assign specific colors to "Essential" and "Non-essential"
  scale_color_manual(values = c("Cell-type specific" ="#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Functional Probability vs min gRNA Depletion p_value Liang 2024", # Your title
    x = "Minimum p value(-log10 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )

print(plot_scatter)


# geom_point creates a scatter plot. Color is mapped to the new Essentiality_Group factor.
plot_scatter <- ggplot(data = essentials_data_all,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -day14_min_calculated, y = Probability_Functional, color = group_day14)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Add smooth trend lines (linear model 'lm' for each group)
  # se = TRUE adds the confidence interval ribbon (default)
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) + # Added trend line layer
  # Optional: Apply a specific color scale (e.g., ColorBrewer Set1)
  # You might want scale_color_manual() to assign specific colors to "Essential" and "Non-essential"
  scale_color_manual(values = c("Cell-type specific" ="#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Functional Probability vs gRNA Depletion Liang 2024", # Your title
    x = "Minimum depletion value at day 14 (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )

print(plot_scatter)

#Sum
#essentials_data_w_filtered <- essentials_data_w %>%
#  filter(sum_all_calculated < 0)
# 4. Scatter Plot (-log2 Sum All vs. Probability Functional, colored by Essentiality_Group)
# geom_point creates a scatter plot. Color is mapped to the new Essentiality_Group factor.
plot_scatter <- ggplot(data = essentials_data_all,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -sum_all_calculated, y = Probability_Functional, color = group_day14)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Add smooth trend lines (linear model 'lm' for each group)
  # se = TRUE adds the confidence interval ribbon (default)
  geom_smooth(method = "loess", span = 0.4, se = FALSE, linewidth = 0.8) +
  # Add RIBBON for DOTTED confidence interval BOUNDS only
  # stat = "smooth" calculates the smooth fit and CI
  # geom = "ribbon" is used implicitly by geom_ribbon
  # aes(ymin/ymax = after_stat(...)) maps the calculated CI bounds
  # fill = NA prevents filling the ribbon area, leaving only the lines
  # linetype = "dotted" makes the boundary lines dotted
  # show.legend = FALSE prevents adding duplicate legend entries
  geom_ribbon(stat = "smooth", method = "loess", span = 0.4, se = TRUE,
              aes(ymin = after_stat(ymin), ymax = after_stat(ymax), color = group_day14), # Map ymin/ymax, inherit color
              fill = NA, # Don't fill the ribbon area
              linetype = "dotted", # Make the boundary lines dotted
              linewidth = 0.6, # Control line thickness
              alpha = 0.5, # Control transparency
              show.legend = FALSE) + # Hide from legend
  # Apply a specific color scale
  scale_color_manual(values = c("Cell-type specific" ="#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Functional Probability vs gRNA Depletion Liang 2024", # Your title
    x = "Sum depletion value day 14 and day7 (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )

print(plot_scatter)



##################################
# Filter Non-essential observations
essential_data_positives <- essentials_data_all %>%
  filter(Essentiality_Group != "Non-essential")

plot_scatter <- ggplot(data = essential_data_positives,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -sum_all_calculated, y = Probability_Functional, color = group_day14)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Optional: Apply a specific color scale (e.g., ColorBrewer Set1)
  # You might want scale_color_manual() to assign specific colors to "Essential" and "Non-essential"
  scale_color_manual(values = c("Cell-type specific" ="#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Functional Probability vs gRNA Depletion Liang 2024", # Your title
    x = "Sum depletion value day 14 and day7 (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )

print(plot_scatter)


# geom_point creates a scatter plot. Color is mapped to the new Essentiality_Group factor.
plot_scatter <- ggplot(data = essentials_data_w,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -sum_all_calculated, y = Probability_Functional, color = Essentiality_Group)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Optional: Apply a specific color scale (e.g., ColorBrewer Set1)
  # You might want scale_color_manual() to assign specific colors to "Essential" and "Non-essential"
  scale_color_manual(values = c("Essential" = "#b30000", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Functional Probability vs gRNA Depletion Liang 2024", # Your title
    x = "Depletion sum(day14 + day7) (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )

print(plot_scatter)

# geom_point creates a scatter plot. Color is mapped to the new Essentiality_Group factor.
plot_scatter <- ggplot(data = essentials_data_w,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -HAP1.Day14.Replicate.1, y = Probability_Functional, color = Essentiality_Group)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Optional: Apply a specific color scale (e.g., ColorBrewer Set1)
  # You might want scale_color_manual() to assign specific colors to "Essential" and "Non-essential"
  scale_color_manual(values = c("Essential" = "#b30000", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Functional Probability vs gRNA Depletion Liang 2024", # Your title
    x = "HAP1 day14 fold change (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )

print(plot_scatter)






# 20 - Jun - 2025
library(ggsignif)
# Define the desired order for the categorical variable
essentiality_order <- c("Non-essential", "Cell-type specific", "Partially shared", "Shared")

# Reorder the 'group_day14' column by converting it to a factor with specified levels
essentials_data_all <- essentials_data_all %>%
  mutate(group_day14 = factor(group_day14, levels = essentiality_order))

# Generate the Violin Plot
plot_violin <- ggplot(data = essentials_data_all,
                      # For a violin plot, we map the categorical variable to x
                      # and the continuous variable to y.
                      # We use 'fill' to color the violins based on the group.
                      aes(x = group_day14, y = Probability_Functional, fill = group_day14)) +
  
  # Add the violin layer
  # trim = FALSE ensures the violin tails are drawn to the full range of the data
  geom_violin(scale = "width", na.rm = TRUE, adjust = 1, trim = FALSE) +
  
  
  # (Optional but recommended) Add a box plot inside the violins
  # This shows the median, interquartile range, and outliers.
  geom_boxplot(alpha=0.0, outliers=TRUE, na.rm = TRUE, position = position_dodge(width = 0.9), width=0.2) +
  
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  # The color values are the same as your original plot.
  scale_fill_manual(values = c("Cell-type specific" = "#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  
  # Update the labels for the new plot layout
  labs(
    title = "Distribution of Functional Probability by Essentiality Group",
    x = "Essentiality Group",
    y = "Functional Probability",
    fill = "Essentiality Group" # Update legend title to match 'fill'
  ) +
  
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    text = element_text(size = 16), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(), # Hiding minor grid lines for a cleaner look
    legend.position = "right"
  )


# Define the pairwise comparisons to be performed
# This creates all unique pairs from your ordered groups.
my_comparisons <- combn(essentiality_order, 2, simplify = FALSE)

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
    textsize = 3.5,
    tip_length = 0.01,
    y_position = 1.2      # Manually set the y-position for the first bracket
  ) +
  # Expand the y-axis slightly to make room for the highest comparison brackets
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))


# To display the final plot
print(plot_with_stats)





# 12 - Jun - 2025
dplyr_left_join_df <- left_join(essentials_data_w, essentials_exon2_data_w, by = "GeneID")

final_df <- dplyr_left_join_df %>% 
  dplyr::select(day14_min_calculated.x, day14_min_calculated.y, sum_all_calculated.x,sum_all_calculated.y,Probability_Functional.x,Probability_Functional.y,group_day14.x,group_day14.y) %>%
  mutate(max_prob_functional = pmax(Probability_Functional.x,Probability_Functional.y, na.rm = TRUE),
         # Create the new column `sum_max_prob`
         # `case_when` is great for defining logic with multiple conditions.
         sum_max_prob = case_when(
           is.na(Probability_Functional.y) | Probability_Functional.x >= Probability_Functional.y ~ sum_all_calculated.x,
           Probability_Functional.x < Probability_Functional.y ~ sum_all_calculated.y,
           TRUE ~ NA_real_ 
         ),
         day14_min = case_when(
           is.na(Probability_Functional.y) | Probability_Functional.x >= Probability_Functional.y ~ day14_min_calculated.x,
           Probability_Functional.x < Probability_Functional.y ~ day14_min_calculated.y,
           TRUE ~ NA_real_
         ),
         day14_group = case_when(
           is.na(Probability_Functional.y) | Probability_Functional.x >= Probability_Functional.y ~ group_day14.x,
           Probability_Functional.x < Probability_Functional.y ~ group_day14.y,
           TRUE ~ NA_character_
         ),
         
         )

# 4. Scatter Plot (-log2 Sum All vs. Probability Functional, colored by Essentiality_Group)
# geom_point creates a scatter plot. Color is mapped to the new Essentiality_Group factor.
plot_scatter <- ggplot(data = final_df,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -day14_min, y = max_prob_functional, color = day14_group)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Add smooth trend lines (linear model 'lm' for each group)
  # se = TRUE adds the confidence interval ribbon (default)
  geom_smooth(method = "loess", span = 0.4, se = FALSE, linewidth = 0.8) +
  # Apply a specific color scale
  scale_color_manual(values = c("Cell-type specific" ="#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Max Functional Probability vs gRNA Depletion Liang 2024", # Your title
    x = "Minimum depletion value at day 14 (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )
#help(geom_smooth)
print(plot_scatter)


plot_scatter <- ggplot(data = final_df,
                       # Map x, y, and color aesthetics to the new factor
                       aes(x = -sum_max_prob, y = max_prob_functional, color = day14_group)) +
  # Add points
  geom_point(alpha = 1, size = 1.5) +
  # Add smooth trend lines (linear model 'lm' for each group)
  # se = TRUE adds the confidence interval ribbon (default)
  geom_smooth(method = "loess", span = 0.4, se = FALSE, linewidth = 0.8) +
  # Apply a specific color scale
  scale_color_manual(values = c("Cell-type specific" ="#1a53ff", "Shared" = "#b30000", "Partially shared" = "#87bc45", "Non-essential" = "#beb9db")) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    title = "Max Functional Probability vs gRNA Depletion Liang 2024", # Your title
    x = "Sum depletion value (-log2 Scale)", # Your x-axis label
    y = "Functional Probability", # Your y-axis label
    color = "Essentiality Group" # Updated legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), # Your text size
    panel.grid.major = element_line(color = "gray90"), # Your grid lines
    panel.grid.minor = element_line(color = "gray95"), # Your grid lines
    legend.position = "right" # Optional: Move legend to bottom
  )
#help(geom_smooth)
print(plot_scatter)
