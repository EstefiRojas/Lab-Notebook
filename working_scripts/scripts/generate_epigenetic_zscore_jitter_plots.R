# Create jitter plots from subset of epigenetic z-score data.

# Load necessary libraries.
library(ggplot2)
library(tidyr)

source("scripts/config.R")
source("scripts/load_epigenetic_zscores.R")
source("scripts/utils.R")

zscores_all <- load_epigenetic_zscores()
zscores_subset <- load_epigenetic_zscores(subset = TRUE)

#If the plots need to be generated with all data instead of the random subset, comment the following line:
zscores_all <- zscores_subset

# Define labels for all epigenetic features of interest.
SELECTED_HISTONE_LABELS <- c(
  "H2AFZ", "H2AK5ac", "H2AK9ac", "H2BK120ac", "H2BK12ac", "H2BK15ac", "H2BK20ac", "H2BK5ac", "H3F3A", "H3K14ac",
  "H3K18ac", "H3K23ac", "H3K23me2", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4ac", "H3K4me1", "H3K4me2", "H3K4me3",
  "H3K56ac", "H3K79me1", "H3K79me2", "H3K9ac", "H3K9me1", "H3K9me2", "H3K9me3", #"H3T11ph",
  "H4K12ac", "H4K20me1", "H4K5ac", "H4K8ac", "H4K91ac", "chrm_acc", "methylome"
)

# Rename column names of zscores_all dataframe.
colnames(zscores_all) <- c(SELECTED_HISTONE_LABELS, "Dataset")

# Melt data for plotting.
ALL_FEATURES <- colnames(zscores_all %>% dplyr::select(-Dataset))

df_long_prot <- melt_zscore_data(
  zscores_all %>% filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")) %>% dplyr::select(-Dataset),
  ALL_FEATURES
)

df_long_lncrna <- melt_zscore_data(
  zscores_all %>% filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2")) %>% dplyr::select(-Dataset),
  ALL_FEATURES
)

df_long_sncrna <- melt_zscore_data(
  zscores_all %>% filter(Dataset == "short-ncrna") %>% dplyr::select(-Dataset),
  ALL_FEATURES
)

df_long_prot_nc <- melt_zscore_data(
  zscores_all %>% filter(Dataset %in% c("protein-exon2-negative-control", "protein-exon3-negative-control")) %>% dplyr::select(-Dataset),
  ALL_FEATURES
)

df_long_lncrna_nc <- melt_zscore_data(
  zscores_all %>% filter(Dataset %in% c("lncrna-exon1-negative-control", "lncrna-exon2-negative-control")) %>% dplyr::select(-Dataset),
  ALL_FEATURES
)

df_long_sncrna_nc <- melt_zscore_data(
  zscores_all %>% filter(Dataset == "short-ncrna-negative-control") %>% dplyr::select(-Dataset),
  ALL_FEATURES
)


# Create a list of dataframes
list_df <- list(
  df_long_prot, df_long_lncrna, df_long_sncrna,
  df_long_prot_nc, df_long_lncrna_nc, df_long_sncrna_nc
)

# Create a list of names for the dataframes
list_names <- c(
  "mRNA (+)", "lncRNA (+)", "sncRNA (+)",
  "mRNA (-)", "lncRNA (-)", "sncRNA (-)"
)

# Add names to the dataframes in the list
names(list_df) <- list_names

# Combine all dataframes into one
combined_df <- dplyr::bind_rows(list_df, .id = "Dataset")
unique(zscores_all$Dataset)
unique(combined_df$Dataset)

# Convert Dataset to a factor with specified order
combined_df$Dataset <- factor(combined_df$Dataset,
  levels = c(
    "mRNA (+)", "lncRNA (+)", "sncRNA (+)",
    "mRNA (-)", "lncRNA (-)", "sncRNA (-)"
  )
)

# Create the plot
create_jitter_plot(combined_df, EPIGENETIC_ZSCORES_JITTER_PLOT_FILE)
