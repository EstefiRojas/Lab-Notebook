# Generate jitter plots of intrinsic features using a subset of intrinsic
# data. This script selects 100 positive cases and 1000 negative controls for each
# RNA type and their respective exons analyzed. To replicate the plots published in 
# the paper, all intrinsic data should be selected (i.e. no sub-setting).

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load configuration and data
source("scripts/config.R")
source("scripts/load_dinucleotide_features.R")
feature_matrix <- load_dinucleotide_features()

# Define features of interest
SELECTED_INTRINSIC_FEATURES <- c("AA", "AC", "AG",
                               "AT", "CA", "CC",
                               "CpG", "CT", "GA", 
                               "GC",
                               "GG", "GT", "TA",
                               "TC", "TG", "TT",
                               "GC_percentage", "lowComplexity_density")

# Select din features and convert to numeric
data_numeric <- feature_matrix %>%
  dplyr::select(all_of(c("Dataset", SELECTED_INTRINSIC_FEATURES))) %>%
  mutate(across(all_of(SELECTED_INTRINSIC_FEATURES), ~as.numeric(as.character(.))))

# Select a subset per RNA type. Choose 100 for positive cases, 1000 for negative controls.
set.seed(2025)  # Set seed for reproducibility
data_numeric_sample <- rbind(
  data_numeric %>%
    filter(Dataset %in% c("protein-coding-exon2")) %>%
    sample_n(100),
  data_numeric %>%
    filter(Dataset %in% c("protein-coding-exon3")) %>%
    sample_n(100),
  data_numeric %>%
    filter(Dataset %in% c("lncrna-exon1")) %>%
    sample_n(100),
  data_numeric %>%
    filter(Dataset %in% c("lncrna-exon2")) %>%
    sample_n(100),
  data_numeric %>%
    filter(Dataset %in% c("short-ncrna")) %>%
    sample_n(100),
  data_numeric %>%
    filter(Dataset %in% c("protein-exon2-negative-control")) %>%
    sample_n(1000),
  data_numeric %>%
    filter(Dataset %in% c("protein-exon3-negative-control")) %>%
    sample_n(1000),
  data_numeric %>%
    filter(Dataset %in% c("lncrna-exon1-negative-control")) %>%
    sample_n(1000),
  data_numeric %>%
    filter(Dataset %in% c("lncrna-exon2-negative-control")) %>%
    sample_n(1000),
  data_numeric %>%
    filter(Dataset %in% c("short-ncrna-negative-control")) %>%
    sample_n(1000)
)
# If the paper plot is desired, un-comment the following line, no other change needed.
#data_numeric_sample <- data_numeric

# Rename datasets for plotting
data_numeric_sample <- data_numeric_sample %>%
  mutate(Dataset = case_when(
    Dataset == "protein-coding-exon2" ~ "mRNA-exon2(+)",
    Dataset == "protein-coding-exon3" ~ "mRNA-exon3(+)",
    Dataset == "lncrna-exon1" ~ "lncRNA-exon1(+)",
    Dataset == "lncrna-exon2" ~ "lncRNA-exon2(+)",
    Dataset == "short-ncrna" ~ "sncRNA(+)",
    Dataset == "protein-exon2-negative-control" ~ "mRNA-exon2(-)",
    Dataset == "protein-exon3-negative-control" ~ "mRNA-exon3(-)",
    Dataset == "lncrna-exon1-negative-control" ~ "lncRNA-exon1(-)",
    Dataset == "lncrna-exon2-negative-control" ~ "lncRNA-exon2(-)",
    Dataset == "short-ncrna-negative-control" ~ "sncRNA(-)",
    TRUE ~ as.character(Dataset)
  ))

# Reshape data for plotting
data_long <- data_numeric_sample %>%
  pivot_longer(cols = all_of(SELECTED_INTRINSIC_FEATURES), names_to = "Feature", values_to = "Value")

# Create directory for plots if it doesn't exist
dir.create(INTRINSIC_JITTER_PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Define a color map for consistent coloring
color_map <- c(
  "mRNA-exon2(+)" = "#8B0000",      # Dark Red (warm)
  "mRNA-exon2(-)" = "#00008B",      # Dark Blue (cold)
  "mRNA-exon3(+)" = "#FF8C00",      # Dark Orange (warm)
  "mRNA-exon3(-)" = "#006400",      # Dark Green (cold)
  "lncRNA-exon1(+)" = "#DC143C",    # Crimson (warm)
  "lncRNA-exon1(-)" = "#4B0082",    # Indigo (cold)
  "lncRNA-exon2(+)" = "#FF4500",    # OrangeRed (warm)
  "lncRNA-exon2(-)" = "#4682B4",    # SteelBlue (cold)
  "sncRNA(+)" = "#FFD700",        # Gold (warm)
  "sncRNA(-)" = "#800080"         # Purple (cold)
)

# Generate and save jitter plots for each feature
for (feature_name in SELECTED_INTRINSIC_FEATURES) {
  p <- ggplot(data_numeric_sample, aes(x = Dataset, y = .data[[feature_name]], color = Dataset)) +
    geom_jitter(width = 0.25, alpha = 0.6) +
    scale_color_manual(values = color_map) +
    theme_minimal() +
    labs(title = paste("Jitter plot for", feature_name),
         x = "Dataset",
         y = "Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(filename = file.path(INTRINSIC_JITTER_PLOTS_DIR, paste0(feature_name, "_jitter.png")),
         plot = p,
         bg = "white",
         width = 10,
         height = 6)
}

# Generate and save a faceted plot of all features
p_faceted <- ggplot(data_long, aes(x = Dataset, y = Value, color = Dataset)) +
  geom_jitter(width = 0.25, alpha = 0.6) +
  facet_wrap(~ Feature, scales = "free_y") +
  scale_color_manual(values = color_map) +
  theme_minimal() +
  labs(title = "Jitter plots for Intrinsic Features",
       x = "Dataset",
       y = "Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(filename = file.path(INTRINSIC_JITTER_PLOTS_DIR, "all_features_jitter.png"),
       plot = p_faceted,
       bg = "white",
       width = 16,
       height = 12)

cat("Jitter plots saved to", INTRINSIC_JITTER_PLOTS_DIR, "\n")
