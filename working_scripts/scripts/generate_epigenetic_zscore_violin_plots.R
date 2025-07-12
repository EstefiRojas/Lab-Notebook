# Create violin plots from subset of epigenetic z-score data.
# Load necessary libraries.
library(ggplot2)
library(tidyr)

source("scripts/config.R")
source("scripts/load_epigenetic_zscores.R")
source("scripts/utils.R")

zscores_all <- load_epigenetic_zscores()
zscores_subset <- load_epigenetic_zscores(subset = TRUE)

# If the plots need to be generated with all data instead of the random subset, comment the following line:
zscores_all <- zscores_subset

# Define labels for all epigenetic features of interest.
SELECTED_HISTONE_LABELS <- c("H2AFZ", "H2AK5ac", "H2AK9ac",
                             "H2BK120ac", "H2BK12ac", "H2BK15ac",
                             "H2BK20ac", "H2BK5ac", "H3F3A",
                             "H3K14ac",
                             "H3K18ac", "H3K23ac", "H3K23me2",
                             "H3K27ac", "H3K27me3", "H3K36me3",
                             "H3K4ac", "H3K4me1", "H3K4me2",
                             "H3K4me3", "H3K56ac", "H3K79me1",
                             "H3K79me2", "H3K9ac", "H3K9me1",
                             "H3K9me2", "H3K9me3", #"H3T11ph",
                             "H4K12ac", "H4K20me1", "H4K5ac",
                             "H4K8ac", "H4K91ac", "chrm_acc",
                             "methylome")

# Rename column names of zscores_all dataframe.
colnames(zscores_all) <- c(SELECTED_HISTONE_LABELS,"Dataset")

# Melt data for plotting.
ALL_FEATURES <- colnames(zscores_all %>% dplyr::select(-Dataset))
df_long_prot <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")) %>%
                                   dplyr::select(-Dataset), ALL_FEATURES)
df_long_lncrna <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2")) %>%
                                     dplyr::select(-Dataset), ALL_FEATURES)
df_long_sncrna <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("short-ncrna")) %>%
                                     dplyr::select(-Dataset), ALL_FEATURES)

df_long_prot_neg <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("protein-exon2-negative-control", "protein-exon3-negative-control")) %>%
                                       dplyr::select(-Dataset), ALL_FEATURES)
df_long_lncrna_neg <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("lncrna-exon1-negative-control", "lncrna-exon2-negative-control")) %>%
                                         dplyr::select(-Dataset), ALL_FEATURES)
df_long_sncrna_neg <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("short-ncrna-negative-control")) %>%
                                         dplyr::select(-Dataset), ALL_FEATURES)

# Create factor for grouping.
df_long_prot$`Gene type` <- "mRNA\n(+)"
df_long_sncrna$`Gene type` <- "sncRNA\n(+)"
df_long_lncrna$`Gene type` <- "lncRNA\n(+)"
df_long_prot_neg$`Gene type` <- "mRNA\n(-)"
df_long_sncrna_neg$`Gene type` <- "sncRNA\n(-)"
df_long_lncrna_neg$`Gene type` <- "lncRNA\n(-)"

df_long_all <- rbind(df_long_prot, df_long_prot_neg,
                     df_long_sncrna, df_long_sncrna_neg,
                     df_long_lncrna, df_long_lncrna_neg)

# OR Convert 'feature' to a factor and specify the desired order
df_long_all$feature <- factor(df_long_all$feature, 
                              levels = ALL_FEATURES, 
                              labels = ALL_FEATURES)

df_long_all$`Gene type` <- factor(df_long_all$`Gene type`, levels = c("mRNA\n(+)",
                                                                      "mRNA\n(-)",
                                                                      "sncRNA\n(+)",
                                                                      "sncRNA\n(-)",
                                                                      "lncRNA\n(+)",
                                                                      "lncRNA\n(-)"))

# Generate violin plot for the first 6 most informative features (i.e. highest ks D value).
create_violin_plot(df_long_all %>% dplyr::filter(feature %in% c("H3K36me3",
                                                                "H3K4me3",
                                                                "H3K79me1",
                                                                "H3K79me2",
                                                                "H3K9ac",
                                                                "H4K20me1")), c(-5,30), "Epigenetic_p1_1", plot_path = EPIGENETIC_P1_1_PLOT_FILE)
create_violin_plot(df_long_all %>% dplyr::filter(feature %in% c("H3K36me3",
                                                                "H3K4me3",
                                                                "H3K79me1",
                                                                "H3K79me2",
                                                                "H3K9ac",
                                                                "H4K20me1")), c(-50,250), "Epigenetic_p1_2", plot_path = EPIGENETIC_P1_2_PLOT_FILE)

# Generate violin plot for the second 6 most informative features.
create_violin_plot(df_long_all %>% dplyr::filter(feature %in% c("H3K4me1",
                                                                "chrm_acc",
                                                                "H3K27ac",
                                                                "H3K4me2",
                                                                "methylome",
                                                                "H3F3A")), c(-5,10), "Epigenetic_p2_1", plot_path = EPIGENETIC_P2_1_PLOT_FILE)
create_violin_plot(df_long_all %>% dplyr::filter(feature %in% c("H3K4me1",
                                                                "chrm_acc",
                                                                "H3K27ac",
                                                                "H3K4me2",
                                                                "methylome",
                                                                "H3F3A")), c(-10,50), "Epigenetic_p2_2", plot_path = EPIGENETIC_P2_2_PLOT_FILE)

# Generate violin plots for the rest of features in two groups of 9.
create_violin_plot(df_long_all %>% dplyr::filter(feature %in% c("H3K9me3",
                                                                "H3K18ac",
                                                                "H3K14ac",
                                                                "H3K23ac",
                                                                "H2BK5ac",
                                                                "H2AFZ",
                                                                "H3K56ac",
                                                                "H3K27me3",
                                                                "H4K8ac")), c(-10,15), "Epigenetic_p3_1", plot_path = EPIGENETIC_P3_1_PLOT_FILE)

create_violin_plot(df_long_all %>% dplyr::filter(feature %in% c("H3K4ac",
                                                                "H2AK5ac",
                                                                "H4K91ac",
                                                                "H2AK9ac",
                                                                "H2BK120ac",
                                                                "H2BK12ac",
                                                                "H2BK15ac",
                                                                "H4K5ac",
                                                                "H3K23me2")), c(-10,15), "Epigenetic_p4_1", plot_path = EPIGENETIC_P4_1_PLOT_FILE)

# Generate violin plot for the 4 least informative features.
create_violin_plot(df_long_all %>% dplyr::filter(feature %in% c("H4K12ac",
                                                                "H3K9me1",
                                                                "H2BK20ac",
                                                                "H3K9me2")), c(-10,30), "Epigenetic_p5_1", plot_path = EPIGENETIC_P5_1_PLOT_FILE)