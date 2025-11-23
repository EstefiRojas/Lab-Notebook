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

# --- MODIFICATION: Group by Study AND GeneID ---
# We must include 'Study' in the grouping to ensure we keep the best candidate 
# per gene *within each study*, rather than globally.
filtered_df <- ess_d %>%
  group_by(Study, ENSG_ID) %>%
  slice_max(order_by = Probability_Functional, n = 1, with_ties = FALSE) %>%
  ungroup()

# --- MODIFICATION: Calculate KS stats PER STUDY ---

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
    
    # Only run test if we have data points
    if(length(comp_data) > 0 && length(ref_data) > 0) {
      ks_result <- ks.test(ref_data, comp_data)
      
      return(data.frame(
        Study = study_name, # Important: Include Study column for faceting
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
  do.call(rbind, .)


# --- Generate the Faceted Plot ---

plot_modified <- ggplot(data = filtered_df,
                        aes(x = Essentiality, y = Probability_Functional, 
                            fill = Essentiality, color = Essentiality)) +
  
  # 1. Add box plot (excluding Core)
  geom_boxplot(data = ~ subset(., Essentiality != "Core"), linewidth = 0.9,
               na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  
  # 2. Add points (for Core)
  geom_point(data = ~ subset(., Essentiality == "Core"),
             size = 4, # Reduced size slightly to fit facets
             shape = 21, 
             color = "black",
             stroke = 1) +
  
  # 3. Colors
  scale_fill_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Set2") +
  
  # 4. Facet by Study
  facet_wrap(~ Study, scales = "fixed") + # Use scales="free_x" if x-groups vary wildly
  
  # 5. Legend guides
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5, color = "black"))) +
  
  # 6. Labels
  labs(
    title = "Gene Essentiality by Study",
    y = "lncRNA Probability",
    fill = "Essentiality Group"
  ) +
  
  # 7. Theme (Adjusted text sizes for faceted view)
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 30, hjust = 0.5),
    strip.text = element_text(size = 24, face = "bold"), # Style for Study headers
    axis.text.x = element_text(hjust = 0.5, size = 16),
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
    aes(x = Essentiality, y = y_position, label = label),
    inherit.aes = FALSE, 
    size = 5, # Reduced text size to fit panels
    color = "black",
    fontface = "bold"
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.2), clip = "off") 

# Display the final plot
print(plot_final)
