# Load libraries
library(dplyr)
library(stringr) 
library(ggplot2)
library(ggsignif)

# Load csv file
essentials1_lncrna_data <- read.csv("../results/unified_genome_alignments.csv", header = TRUE)

ess_d <- essentials1_lncrna_data %>% filter(Protein_Off_Target == "NO")

# Define the desired order for the categorical variable
essentiality_values <- c("Non-essential", "Specific", "Common", "Core")
essentiality_labels <- c("Non", "Specific", "Common", "Core")

# Reorder the 'Essentiality' column
ess_d <- ess_d %>%
  mutate(Essentiality = factor(Essentiality, levels = essentiality_values, labels = essentiality_labels))

# Group by Study AND GeneID
filtered_df <- ess_d %>%
  group_by(Study, ENSG_ID) %>%
  slice_max(order_by = Probability_Functional, n = 1, with_ties = FALSE) %>%
  ungroup()

# --- MODIFICATION: Create Labels with Counts per Study ---
# 1. Calculate counts per Study + Essentiality
counts_df <- filtered_df %>%
  group_by(Study, Essentiality) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Label_With_Count = paste0(Essentiality, "\n(n=", n, ")"))

# 2. Join back to filtered_df
filtered_df <- filtered_df %>%
  left_join(counts_df, by = c("Study", "Essentiality"))

# 3. Ensure the new labels are ordered correctly (Non -> Specific -> Common -> Core)
# We sort by the original Essentiality factor to define the factor levels for the new labels
label_levels <- filtered_df %>%
  arrange(Study, Essentiality) %>%
  pull(Label_With_Count) %>%
  unique()

filtered_df$Label_With_Count <- factor(filtered_df$Label_With_Count, levels = label_levels)


# --- Calculate KS stats PER STUDY ---

# Define groups
reference_group <- "Non"
comparison_groups <- c("Specific", "Common", "Core")

# Create a function to compute stats for a specific study dataframe
compute_study_stats <- function(df_subset, study_name) {
  
  # Get reference data for this study
  ref_data <- df_subset %>%
    filter(Essentiality == reference_group) %>%
    pull(Probability_Functional)
  
  # Loop through comparison groups
  stats_list <- lapply(comparison_groups, function(group) {
    comp_data <- df_subset %>%
      filter(Essentiality == group) %>%
      pull(Probability_Functional)
    
    # Only run test if we have data points in both groups
    if(length(comp_data) > 0 && length(ref_data) > 0) {
      ks_result <- ks.test(ref_data, comp_data)
      
      return(data.frame(
        Study = study_name, 
        Essentiality = group,
        label = paste0("KS=", round(ks_result$statistic, 2)),
        y_position = 1.1 
      ))
    } else {
      return(NULL)
    }
  })
  
  do.call(rbind, stats_list)
}

# Apply the function to each study and combine results
all_stats_labels <- filtered_df %>%
  group_split(Study) %>%
  lapply(function(subset_df) {
    compute_study_stats(subset_df, unique(subset_df$Study))
  }) %>%
  do.call(rbind, .) %>%
  # Join with counts_df to get the correct Label_With_Count for placing the text
  left_join(counts_df, by = c("Study", "Essentiality"))


# --- Generate the Faceted Plot ---

plot_modified <- ggplot(data = filtered_df,
                        # Update x to use the new label with counts
                        aes(x = Label_With_Count, y = Probability_Functional, 
                            fill = Essentiality, color = Essentiality)) +
  
  # 1. Add box plot (excluding Core)
  geom_boxplot(data = ~ subset(., Essentiality != "Core"), linewidth = 0.9,
               na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  
  # 2. Add points (for Core)
  geom_point(data = ~ subset(., Essentiality == "Core"),
             size = 4,
             shape = 21, 
             color = "black",
             stroke = 1) +
  
  # 3. Colors
  scale_fill_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Set2") +
  
  # 4. Facet by Study with scales="free_x"
  # This allows each facet to have its own unique x-axis labels (with different counts)
  facet_wrap(~ Study, scales = "free_x") + 
  
  # 5. Legend guides
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5, color = "black"))) +
  
  # 6. Labels
  labs(
    title = "Gene Essentiality by Study",
    y = "lncRNA Probability",
    fill = "Essentiality Group"
  ) +
  
  # 7. Theme
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 30, hjust = 0.5),
    strip.text = element_text(size = 24, face = "bold"), 
    # Angle labels slightly if they get too long, or keep flat if space allows
    axis.text.x = element_text(hjust = 0.5, size = 14), 
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 24),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3)
  )

# Add the statistical comparison layer
plot_final <- plot_modified +
  geom_text(
    data = all_stats_labels,
    # Update x to use the new label with counts for text placement
    aes(x = Label_With_Count, y = y_position, label = label),
    inherit.aes = FALSE, 
    size = 5,
    color = "black",
    fontface = "bold"
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.2), clip = "off") 

# Display the final plot
print(plot_final)
