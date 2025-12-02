# Load libraries
library(dplyr)
library(stringr) 
library(ggplot2)
library(ggsignif)
#install.packages("ggforce")
library(ggforce)
#install.packages("ggbeeswarm")
library(ggbeeswarm)

# Load csv file
essentials1_lncrna_data <- read.csv("../results/annotated_unified_genome_alignments.csv", header = TRUE)

ess_d <- essentials1_lncrna_data %>% filter(Protein_Off_Target == "NO" & !is.na(Probability_Functional) & Antisense_to_CDS == "NO")

# Define the desired order for the categorical variable
essentiality_values <- c("Non-essential", "Rare", "Common", "Core")
essentiality_labels <- c("Non", "Rare", "Common", "Core")

# Reorder the 'Essentiality' column
ess_d <- ess_d %>%
  mutate(Essentiality = factor(Essentiality, levels = essentiality_values, labels = essentiality_labels))

# Group by Study AND GeneID
filtered_df <- ess_d %>%
  group_by(ENSG_ID, Study, Essentiality) %>%
  slice_max(order_by = Probability_Functional, n = 1, with_ties = FALSE) %>%
  ungroup()




# --- MODIFICATION: Create Unique Labels with Counts per Study ---
# 1. Calculate counts per Study + Essentiality
counts_df <- filtered_df %>%
  group_by(Study, Essentiality) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    # The visible label
    Display_Label = paste0(Essentiality, "\n(n=", n, ")"),
    # A unique key for sorting (Study + Label) to prevent "n=44" in Liu merging with "n=44" in Montero
    Label_Key = paste(Study, Display_Label, sep = "___") 
  )

# 2. Join back to filtered_df
filtered_df <- filtered_df %>%
  left_join(counts_df, by = c("Study", "Essentiality"))

# 3. Ensure the factor levels follow the strict Study -> Essentiality order
ordered_levels <- counts_df %>%
  arrange(Study, Essentiality) %>%
  pull(Label_Key)

filtered_df$Label_Key <- factor(filtered_df$Label_Key, levels = ordered_levels)


# --- Calculate KS stats PER STUDY ---

# Define groups
reference_group <- "Non"
comparison_groups <- c("Rare", "Common", "Core")

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
  # Join with counts_df to get the correct Label_Key for placing the text
  left_join(counts_df, by = c("Study", "Essentiality"))


# --- Generate the Faceted Plot ---
help("geom_beeswarm")
plot_modified <- ggplot(data = filtered_df,
                        # Use the Unique 'Label_Key' for x-axis to maintain sort order
                        aes(x = Label_Key, y = Probability_Functional, 
                            fill = Essentiality, color = Essentiality)) +
  
  # 1. Add box plot (excluding Core)
  geom_boxplot(data = ~ subset(., Essentiality != "Core"), linewidth = 0.9,
               na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  
  # 2. Add points (for Core)
  geom_beeswarm(data = ~ subset(., Essentiality == "Core"),
              size = 4,
              cex = 3,
              shape = 21, 
              color = "black",
              priority = "ascending",
              stroke = 1) +
  
  # 3. Colors
  scale_fill_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Set2") +
  
  # 4. Facet by Study with scales="free_x"
  facet_wrap(~ Study, scales = "free_x") + 
  
  # 5. Clean up the X-axis labels to remove the "Study___" prefix
  scale_x_discrete(labels = function(x) sub(".*___", "", x)) +
  
  # 6. Legend guides
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5, color = "black"))) +
  
  # 7. Labels
  labs(
    title = "Gene Essentiality by Study",
    y = "lncRNA Probability",
    fill = "Essentiality Group"
  ) +
  
  # 8. Theme
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 30, hjust = 0.5),
    strip.text = element_text(size = 24, face = "bold"), 
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
    # Use Label_Key for placement here as well
    aes(x = Label_Key, y = y_position, label = label),
    inherit.aes = FALSE, 
    size = 5,
    color = "black",
    fontface = "bold"
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.2), clip = "off") 

# Display the final plot
print(plot_final)

