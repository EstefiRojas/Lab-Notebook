# Create violin plots for z-score distribution of features.
# Load necessary libraries
library(ggplot2)
library(tidyr)

# Load z-scores
file <- "pasted_histone_outputs/gene-functionality-zscores.csv"
zscores_data <- read.csv(file, header = TRUE, check.names = FALSE)
# OR
zscores_data <- rbind(functional_protein_exon2_dataset_zscores,functional_protein_exon3_dataset_zscores,
                      functional_lncrna_exon1_dataset_zscores, functional_lncrna_exon2_dataset_zscores,
                      functional_short_ncrna_dataset_zscores,
                      protein_exon2_negative_control_dataset_zscores, protein_exon3_negative_control_dataset_zscores,
                      lncrna_exon1_negative_control_dataset_zscores, lncrna_exon2_negative_control_dataset_zscores,
                      short_ncrna_negative_control_dataset_zscores)
# OR
zscores_data <- rbind(protein_functional_z_scores,
                      lncrna_functional_z_scores,
                      sncrna_functional_z_scores,
                      protein_negative_z_scores,
                      lncrna_negative_z_scores,
                      sncrna_negative_z_scores)

# OR
zscores_data <- rbind(read.csv("z_scores/functional_protein_z_scores.csv",header = TRUE,check.names = FALSE),
                      read.csv("z_scores/functional_sncrna_z_scores.csv",header = TRUE,check.names = FALSE))

# Divide data in corresponding gene type category
list_with_zscores_prot <- zscores_data %>% filter(Dataset=="protein-coding-exon2" | Dataset=="protein-coding-exon3")
list_with_zscores_lncrna <- zscores_data %>% filter(Dataset=="lncrna-exon1" | Dataset=="lncrna-exon2")
list_with_zscores_sncrna<- zscores_data %>% filter(Dataset=="short-ncrna")
list_with_zscores_protein_neg <- zscores_data %>% filter(Dataset=="protein-exon2-negative-control" | Dataset=="protein-exon3-negative-control")
list_with_zscores_lncrna_neg <- zscores_data %>% filter(Dataset=="lncrna-exon1-negative-control" | Dataset=="lncrna-exon2-negative-control")
list_with_zscores_sncrna_neg <- zscores_data %>% filter(Dataset=="short-ncrna-negative-control")
# OR
list_with_zscores_prot <- data.frame(Dataset="protein", read.csv("z_scores/functional_protein_z_scores.csv",header = TRUE,check.names = FALSE))
list_with_zscores_lncrna <- data.frame(Dataset="lncrna", read.csv("z_scores/functional_lncrna_z_scores.csv",header = TRUE,check.names = FALSE))
list_with_zscores_sncrna<- data.frame(Dataset="sncrna", read.csv("z_scores/functional_sncrna_z_scores.csv",header = TRUE,check.names = FALSE))
list_with_zscores_protein_neg <- data.frame(Dataset="protein_negative_control", read.csv("z_scores/protein_negative_control_z_scores.csv",header = TRUE,check.names = FALSE))
list_with_zscores_lncrna_neg <- data.frame(Dataset="lncrna_negative_control", read.csv("z_scores/lncrna_negative_control_z_scores.csv",header = TRUE,check.names = FALSE))
list_with_zscores_sncrna_neg <- data.frame(Dataset="sncrna_negative_control", read.csv("z_scores/sncrna_negative_control_z_scores.csv",header = TRUE,check.names = FALSE))
# OR
list_with_zscores_protein <- protein_functional_z_scores
list_with_zscores_protein$Dataset <- "protein"
list_with_zscores_lncrna <- lncrna_functional_z_scores
list_with_zscores_lncrna$Dataset <- "lncrna"
list_with_zscores_sncrna <- sncrna_functional_z_scores
list_with_zscores_sncrna$ Dataset <- "sncrna"

list_with_zscores_protein_neg <- protein_negative_z_scores
list_with_zscores_protein_neg$Dataset <- "protein_negative_control"
list_with_zscores_lncrna_neg <- lncrna_negative_z_scores
list_with_zscores_lncrna_neg$Dataset <- "lncrna_negative_control"
list_with_zscores_sncrna_neg <- sncrna_negative_z_scores
list_with_zscores_sncrna_neg$ Dataset <- "sncrna_negative_control"


colnames(list_with_zscores_protein)
VIOLIN_PLOT_SELECTED_FEATURES <- c("Random.number","GC.content","AA","AC",
                                   "AG","AT","CA","CC",
                                   "CpG","CT","GA","GC",
                                   "GG","GT","TA","TC",
                                   "TG","TT","low_complexity_density","phyloP.max_241w",
                                   "phyloP.max_100w","GERP_91_mammals_max","GERP_63_amniotes_max","RPKM_tissue",
                                   "RPKM_primary.cell","Copy.number","Repeat.free","Fickett_score",
                                   "RNAcode","Max.covariance","MFE","accessibility",
                                   "RNAalifold","Interaction_ave","gnomAD_SNP_density","gnomAD_MAF",
                                   "H3K27ac","H3K36me3","H3K79me2","chromatin_acc",
                                   "methylome")
#OR
VIOLIN_PLOT_SELECTED_FEATURES <- c("H3K9ac_MaxScaledSignal","H3K79me2_MaxScaledSignal",
                                   "chrm_acc_MaxScaledSignal","methylome")

VIOLIN_PLOT_SELECTED_FEATURES_LABELS <- c("Random number",
                                          "GC%",
                                          "AA dinucleotide content","AC dinucleotide content","AG dinucleotide content",
                                          "AT dinucleotide content","CA dinucleotide content","CC dinucleotide content",
                                          "CpG dinucleotide content",
                                          "CT dinucleotide content","GA dinucleotide content","GC dinucleotide content",
                                          "GG dinucleotide content",
                                          "GT dinucleotide content","TA dinucleotide content","TC dinucleotide content",
                                          "TG dinucleotide content",
                                          "TT dinucleotide content",
                                          "Low complexity density","PhyloP max (241 mammals)","PhyloP max (100 vertebrates)",
                                          "GERP max (91 eutherian mammals)","GERP max (63 amniota vertebrates)",
                                          "Tissue RPKM","Primary cell RPKM","Genomic copy number","Repeat free region",
                                          "Fickett score","RNAcode coding potential","Covariance max",
                                          "Secondary structure MFE","RNA accessibility","RNAalifold","RNA-RNA interaction",
                                          "SNP density","Minor allele frequency average",
                                          "H3K27ac","H3K36me3","H3K79me2","Chromatin accessibility","Methylome")
#OR
VIOLIN_PLOT_SELECTED_FEATURES_LABELS <- c("H3K9ac","H3K79me2", "Chromatin Acc","Methylome")

VIOLIN_PLOT_SELECTED_FEATURES_LABELS_SHORT <- c("Random",
                                                "GC%",
                                                "AA","AC","AG","AT","CA","CC",
                                                "CpG","CT","GA","GC",
                                                "GG","GT","TA","TC","TG","TT",
                                                "Complexity","PhyloP-mammals","PhyloP-vertebrates",
                                                "GERP-mammals","GERP-vertebrates",
                                                "Tissue RPKM","Primary cell RPKM","Copies","Repeat free",
                                                "Fickett","RNAcode","Covariance",
                                                "MFE","Accessibility","RNAalifold","Interactions",
                                                "SNPs","MAF",
                                                "H3K27ac","H3K36me3","H3K79me2","Chromatin", "Methylome")
# OR
VIOLIN_PLOT_SELECTED_FEATURES_LABELS_SHORT <- c("H3K9ac","H3K79me2","Chromatin Accessibility","Methylome")

# Function to Melt z-scores for plotting
melt_zscore_data <- function(dataframe, selected_features) {
  plot_list <- list()
  # Create an empty list to store data frames
  long_data_frames <- list()
  
  for (col in selected_features) {
    # Melt z-scores for violin plot
    df_with_zscores <- dataframe %>%
      pivot_longer(cols = col, 
                   names_to = "feature", 
                   values_to = "z_score") %>%
      dplyr::select(feature, z_score)
    
    # Store the result in the list
    long_data_frames[[col]] <- df_with_zscores
  }
  
  # Combine all melted data frames using bind_rows
  df_long <- dplyr::bind_rows(long_data_frames)
  
  return(df_long)
}

# Melt z-scores in a long vector just of selected features
df_long_prot <- melt_zscore_data(list_with_zscores_protein, VIOLIN_PLOT_SELECTED_FEATURES)
df_long_lncrna <- melt_zscore_data(list_with_zscores_lncrna, VIOLIN_PLOT_SELECTED_FEATURES)
df_long_sncrna <- melt_zscore_data(list_with_zscores_sncrna, VIOLIN_PLOT_SELECTED_FEATURES)

df_long_prot_neg <- melt_zscore_data(list_with_zscores_protein_neg, VIOLIN_PLOT_SELECTED_FEATURES)
df_long_lncrna_neg <- melt_zscore_data(list_with_zscores_lncrna_neg, VIOLIN_PLOT_SELECTED_FEATURES)
df_long_sncrna_neg <- melt_zscore_data(list_with_zscores_sncrna_neg, VIOLIN_PLOT_SELECTED_FEATURES)

# OR Melt z-scores in a long vector of all features
ALL_FEATURES <- colnames(zscores_data %>% dplyr::select(-Dataset))
df_long_prot <- melt_zscore_data(list_with_zscores_prot %>% dplyr::select(-Dataset), ALL_FEATURES)
df_long_lncrna <- melt_zscore_data(list_with_zscores_lncrna %>% dplyr::select(-Dataset), ALL_FEATURES)
df_long_sncrna <- melt_zscore_data(list_with_zscores_sncrna %>% dplyr::select(-Dataset), ALL_FEATURES)

df_long_prot_neg <- melt_zscore_data(list_with_zscores_protein_neg %>% dplyr::select(-Dataset), ALL_FEATURES)
df_long_lncrna_neg <- melt_zscore_data(list_with_zscores_lncrna_neg %>% dplyr::select(-Dataset), ALL_FEATURES)
df_long_sncrna_neg <- melt_zscore_data(list_with_zscores_sncrna_neg %>% dplyr::select(-Dataset), ALL_FEATURES)

# Create factor for grouping
df_long_prot$`Gene type` <- "mRNA\n(+)"
df_long_sncrna$`Gene type` <- "sncRNA\n(+)"
df_long_lncrna$`Gene type` <- "lncRNA\n(+)"
df_long_prot_neg$`Gene type` <- "mRNA\n(-)"
df_long_sncrna_neg$`Gene type` <- "sncRNA\n(-)"
df_long_lncrna_neg$`Gene type` <- "lncRNA\n(-)"

df_long_all <- rbind(df_long_prot, df_long_prot_neg,
                     df_long_sncrna, df_long_sncrna_neg,
                     df_long_lncrna, df_long_lncrna_neg)

# Convert 'feature' to a factor and specify the desired order
df_long_all$feature <- factor(df_long_all$feature, 
                              levels = VIOLIN_PLOT_SELECTED_FEATURES, 
                              labels = VIOLIN_PLOT_SELECTED_FEATURES_LABELS_SHORT)
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

create_violin_plot <- function(df_long_format, ylimit, title, violin_adjust = 1, violin_scale = "width") {
  p <- ggplot(df_long_format, aes(x = feature, y = z_score, fill = `Gene type`)) +
    #geom_jitter() + 
    geom_violin(scale = violin_scale, na.rm = TRUE, adjust = violin_adjust, trim = FALSE) +
    geom_boxplot(alpha=0.0, outliers=TRUE, na.rm = TRUE, position = position_dodge(width = 0.9), width=0.2) +
    facet_wrap(~ feature, scales = "free") + 
    labs(title = "Conservation", x = "Feature", y = "Robust Z-score") +
    theme_minimal(base_size = 34) +
    theme(
      axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
      axis.text.y = element_text(size = 20),  # Increase y-axis text size
      axis.title.x = element_text(size = 0),  # Increase x-axis title size
      axis.title.y = element_text(size = 26),  # Increase y-axis title size
      legend.position = "none",
      legend.title = element_text(size = 20),  # Increase legend title size
      legend.text = element_text(size = 18),  # Increase legend text size
      plot.title = element_text(size = 0, hjust = 0.5), # Increase plot title size and center it
      plot.subtitle = element_blank(),
      axis.line.x = element_blank(),          # Remove x-axis line
      axis.ticks.x = element_blank(),         # Remove x-axis ticks
      panel.grid.major.x = element_blank(),   # Remove major grid lines along x-axis
      panel.grid.minor.x = element_blank(),    # Remove minor grid lines along x-axis
      strip.text = element_text(margin = ggplot2::margin(0,0,40,0))
    ) +
    scale_fill_manual(values = c("#F4A582FF","#c9e3f6FF","#D6604DFF","#56bdfcFF","#e37b88FF","#53a4f5FF")) + # colors
    coord_cartesian(ylim = ylimit) # Zoom into an specific area without removing data points
  
  ggsave(paste(title,".png",sep = ""),path = "../results/latest1000all/violinPlots/Subset/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)
  return(p)
}
create_violin_plot(df_long_all %>% dplyr::filter(feature %in% c("H2AFZ_MaxScaledSignal",
                                                                "H3F3A_MaxScaledSignal",
                                                                "H2BK120ac_MaxScaledSignal",
                                                                "H2BK12ac_MaxScaledSignal",
                                                                "H2BK15ac_MaxScaledSignal",
                                                                "H2BK20ac_MaxScaledSignal",
                                                                "H2BK5ac_MaxScaledSignal",
                                                                "H3K27me3_MaxScaledSignal",
                                                                "H3K9me3_MaxScaledSignal",
                                                                "H4K91ac_MaxScaledSignal",
                                                                "H3K4me1_MaxScaledSignal",
                                                                "H3K4me2_MaxScaledSignal")), c(-5,10), "Epigenetic_p1")

create_violin_plot(df_long_all %>% dplyr::filter(feature %in% c("H2AK5ac_MaxScaledSignal",
                                                                "H3K14ac_MaxScaledSignal",
                                                                "H3K9ac_MaxScaledSignal",
                                                                "H3K9me1_MaxScaledSignal",
                                                                "H3K9me2_MaxScaledSignal",
                                                                "H3K18ac_MaxScaledSignal",
                                                                "H3K23ac_MaxScaledSignal",
                                                                "H4K5ac_MaxScaledSignal",
                                                                "H4K8ac_MaxScaledSignal",
                                                                "H3K4ac_MaxScaledSignal")), c(-10,15), "Epigenetic_p2")

create_violin_plot(df_long_all %>% dplyr::filter(feature %in% c("H2AK9ac_MaxScaledSignal",
                                                                "H3K23me2_MaxScaledSignal",
                                                                "H3K27ac_MaxScaledSignal",
                                                                "H3K36me3_MaxScaledSignal",
                                                                "H3K56ac_MaxScaledSignal",
                                                                "H3K79me1_MaxScaledSignal" )), c(-10,25), "Epigenetic_p3")

create_violin_plot(df_long_all %>% dplyr::filter(feature %in% c("H3T11ph_MaxScaledSignal",
                                                                "H4K12ac_MaxScaledSignal",
                                                                "H4K20me1_MaxScaledSignal",
                                                                "H3K4me3_MaxScaledSignal",
                                                                "H3K79me2_MaxScaledSignal",
                                                                "chrm_acc_MaxScaledSignal",
                                                                "methylome")), c(NA,NA), "Epigenetic_p4")
# Create plots
# Conservation:
df_long_conservation_paper <- df_long_all[df_long_all$feature=="PhyloP-mammals" 
                                          | df_long_all$feature=="PhyloP-vertebrates" 
                                          | df_long_all$feature=="GERP-mammals" 
                                          | df_long_all$feature=="GERP-vertebrates",]
create_violin_plot(df_long_conservation_paper, c(-2, 6), "conservation_review")

# Expression:
df_long_expression_paper <- df_long_all[df_long_all$feature=="Tissue RPKM" 
                                        | df_long_all$feature=="Primary cell RPKM", ]
create_violin_plot(df_long_expression_paper, c(-1, 500), "expression_review_test")

# Dinucleotides:
df_long_intrinsic1 <- df_long_all[df_long_all$feature=="GC%" 
                                  | df_long_all$feature=="Complexity"
                                  | df_long_all$feature=="CpG"
                                  | df_long_all$feature=="GG"
                                  | df_long_all$feature=="TA"
                                  | df_long_all$feature=="GA"
                                  | df_long_all$feature=="GT"
                                  | df_long_all$feature=="AC"
                                  | df_long_all$feature=="CC", ]
custom_order <- c("GC%","GA", "CpG", 
                  "TA","GG","GT",
                  "AC","CC", "Complexity")
df_long_intrinsic1$feature <- factor(df_long_intrinsic1$feature, levels = custom_order)
df_long_intrinsic1_ordered <- df_long_intrinsic1 %>%
  arrange(feature)
create_violin_plot(df_long_intrinsic1_ordered, c(0,0.5), "intrinsic2")
summary(list_with_zscores_prot)
summary(list_with_zscores_lncrna)
summary(list_with_zscores_sncrna)

summary(list_with_zscores_protein_neg)
summary(list_with_zscores_lncrna_neg)
summary(list_with_zscores_sncrna_neg)

# Epigenetic:
df_long_epigen_1_paper <- df_long_all[df_long_all$feature=="H3K9ac" 
                                      | df_long_all$feature=="Chromatin Accessibility"
                                      | df_long_all$feature=="H3K79me2"
                                      | df_long_all$feature=="Methylome", ]
create_violin_plot(df_long_epigen_1_paper, c(-5,15), "epigenetic_review")
create_violin_plot(df_long_epigen_1_paper, c(-10,40), "epigenetic_review2")
create_violin_plot(df_long_epigen_1_paper, c(NA,NA), "epigenetic_review3")

#OR
df_long_epigen_1_paper <- df_long_all[df_long_all$feature=="H3K27ac" 
                                      | df_long_all$feature=="H3K36me3" 
                                      | df_long_all$feature=="H3K79me2", ]
create_violin_plot(df_long_epigen_1_paper, c(NA,NA), "epigenetic_maxScaledSignal")
create_violin_plot(df_long_epigen_1_paper, c(0,40), "epigenetic_maxScaledSignal2")
create_violin_plot(df_long_epigen_1_paper, c(-3,15), "epigenetic_maxSacaledSignal3")


# Association:
df_long_repeat_assoc_coding <- df_long_all[df_long_all$feature=="Repeat free" 
                                           | df_long_all$feature=="Copies", ]
custom_order <- c("Repeat free", "Copies")
df_long_repeat_assoc_coding$feature <- factor(df_long_repeat_assoc_coding$feature, levels = custom_order)
df_long_repeat_assoc_coding_ordered <- df_long_repeat_assoc_coding %>%
  arrange(feature)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,5), "association_review_trimtrue", violin_scale = "width")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-1,1000), "association_review2")

# Structure:
df_long_struct_paper <- df_long_all[df_long_all$feature=="Covariance" 
                                    | df_long_all$feature=="MFE" 
                                    | df_long_all$feature=="Interactions"
                                    | df_long_all$feature=="RNAcode"
                                    | df_long_all$feature=="Fickett", ]
custom_order <- c("Fickett","RNAcode","Interactions","Covariance","MFE")
df_long_struct_paper$feature <- factor(df_long_struct_paper$feature, levels = custom_order)
df_long_struct_pape_ordered <- df_long_struct_paper %>%
  arrange(feature)
create_violin_plot(df_long_struct_pape_ordered, c(NA,NA), "structure_review")
create_violin_plot(df_long_struct_pape_ordered, c(-3,17), "structure_review2")
create_violin_plot(df_long_struct_pape_ordered, c(-7,2), "structure_review3")

# Population Variation:
df_long_pop_var <- df_long_all[df_long_all$feature=="MAF"
                               | df_long_all$feature=="SNPs", ]
create_violin_plot(df_long_pop_var, c(-2,10), "popvar_review")
