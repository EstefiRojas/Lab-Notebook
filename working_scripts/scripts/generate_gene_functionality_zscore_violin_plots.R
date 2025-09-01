# Create violin plots from gene functionality z-score data.
# Load necessary libraries.
library(ggplot2)
library(tidyr)

source("scripts/config.R")
source("scripts/load_gene_functionality_zscores.R")
source("scripts/utils.R")
zscores_all <- load_gene_functionality_zscores()

colnames(zscores_all)
# Define names of selected features.
VIOLIN_PLOT_SELECTED_FEATURES <- c("Random","GC_percentage",
                                   "AA","AC","AG","AT","CA","CC",
                                   "CpG","CT","GA","GC",
                                   "GG","GT","TA","TC","TG","TT",
                                   "lowComplexity_density","phyloP_max_241w","phyloP_max_100w",
                                   "GERP_91_mammals_max","GERP_63_amniotes_max",
                                   "RPKM_tissue","RPKM_primary.cell",
                                   "copy_number","repeat_distance",
                                   "fickett","coding_potential","Max_covariance",
                                   "MFE","accessibility","RNAalifold_score","Interaction_ave",
                                   "SNP_density","MAF_avg",
                                   "H3K9ac_MaxScaledSignal","H3K79me2_MaxScaledSignal","H3K79me1_MaxScaledSignal",
                                   "chrm_acc_MaxScaledSignal",
                                   "methylome")
length(VIOLIN_PLOT_SELECTED_FEATURES)

VIOLIN_PLOT_SELECTED_FEATURES_LABELS <- c("Random number","GC%",
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
                                          "H3K9ac","H3K79me2","H3K79me1","Chromatin Accessibility","Methylome")

VIOLIN_PLOT_SELECTED_FEATURES_LABELS_SHORT <- c("Random","GC%",
                                                "AA","AC","AG","AT","CA","CC",
                                                "CpG","CT","GA","GC",
                                                "GG","GT","TA","TC","TG","TT",
                                                "Low Complexity","PhyloP-mammals","PhyloP-vertebrates",
                                                "GERP-mammals","GERP-vertebrates",
                                                "Tissue RPKM","Primary cell RPKM",
                                                "Copies","Repeat free",
                                                "Fickett","RNAcode","Covariance",
                                                "MFE","Accessibility","RNAalifold","Interactions",
                                                "SNPs","MAF",
                                                "H3K9ac","H3K79me2","H3K79me1",
                                                "Chromatin Accessibility",
                                                "Methylome")
#length(VIOLIN_PLOT_SELECTED_FEATURES_LABELS_SHORT)
# OR
# VIOLIN_PLOT_SELECTED_FEATURES_LABELS_SHORT <- c("H3K9ac","H3K79me2","Chromatin Accessibility","Methylome")

# Melt data for plotting.
df_long_prot <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")) %>%
                                   dplyr::select(-Dataset), VIOLIN_PLOT_SELECTED_FEATURES)
df_long_lncrna <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2")) %>%
                                     dplyr::select(-Dataset), VIOLIN_PLOT_SELECTED_FEATURES)
df_long_sncrna <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("short-ncrna")) %>%
                                     dplyr::select(-Dataset), VIOLIN_PLOT_SELECTED_FEATURES)

df_long_prot_neg <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("protein-exon2-negative-control", "protein-exon3-negative-control")) %>%
                                       dplyr::select(-Dataset), VIOLIN_PLOT_SELECTED_FEATURES)
df_long_lncrna_neg <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("lncrna-exon1-negative-control", "lncrna-exon2-negative-control")) %>%
                                         dplyr::select(-Dataset), VIOLIN_PLOT_SELECTED_FEATURES)
df_long_sncrna_neg <- melt_zscore_data(zscores_all %>% filter(Dataset %in% c("short-ncrna-negative-control")) %>%
                                         dplyr::select(-Dataset), VIOLIN_PLOT_SELECTED_FEATURES)


# Create factor for grouping
df_long_prot$`Gene type` <- "mRNA (+)"
df_long_sncrna$`Gene type` <- "sncRNA (+)"
df_long_lncrna$`Gene type` <- "lncRNA (+)"
df_long_prot_neg$`Gene type` <- "mRNA (-)"
df_long_sncrna_neg$`Gene type` <- "sncRNA (-)"
df_long_lncrna_neg$`Gene type` <- "lncRNA (-)"

df_long_all <- rbind(df_long_prot, df_long_prot_neg,
                     df_long_sncrna, df_long_sncrna_neg,
                     df_long_lncrna, df_long_lncrna_neg)

# Convert 'feature' to a factor and specify the desired order
df_long_all$feature <- factor(df_long_all$feature, 
                              levels = VIOLIN_PLOT_SELECTED_FEATURES, 
                              labels = VIOLIN_PLOT_SELECTED_FEATURES_LABELS_SHORT)

# Convert Gene type to factor.
df_long_all$`Gene type` <- factor(df_long_all$`Gene type`, levels = c("mRNA (+)",
                                                                      "mRNA (-)",
                                                                      "sncRNA (+)",
                                                                      "sncRNA (-)",
                                                                      "lncRNA (+)",
                                                                      "lncRNA (-)"))

# Function to generate the violin plots.
create_violin_plot <- function(df_long_format, ylimit, title, violin_adjust = 1, 
                               violin_scale = "width", kernel = "gaussian", bw = "nrd0",
                               width=3840, height=2160) {
  p <- ggplot(df_long_format, aes(x = feature, y = z_score, fill = `Gene type`), color = `Gene type`) +
    #geom_jitter() + 
    geom_violin(scale = violin_scale, na.rm = TRUE, adjust = violin_adjust, trim = FALSE, kernel = kernel, bw = bw, bounds = ylimit) +
    geom_boxplot(alpha=0.0, outliers=TRUE, na.rm = TRUE, position = position_dodge(width = 0.9), width=0.2) +
    # Jitter layer - will now inherit the 'color' mapping and be colored.
    #geom_jitter(
    #  position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2),
    #  alpha = 0.5, # Slightly increased alpha for visibility
    #  size = 0.5,
    #  na.rm = TRUE
    #) +
    facet_wrap(~ feature, scales = "free") + 
    labs(title = "Conservation", x = "Feature", y = "Robust Z-score") +
    theme_minimal(base_size = 34) +
    theme(
      axis.text.x = element_text(size = 0, angle = 45, hjust = 1),
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
    scale_fill_manual(values = c("mRNA (+)"="#F4A582FF","mRNA (-)"="#c9e3f6FF",
                                 "sncRNA (+)"="#D6604DFF","sncRNA (-)"="#56bdfcFF",
                                 "lncRNA (+)"="#e37b88FF","lncRNA (-)"="#53a4f5FF")) + # colors
    scale_color_manual(values = c("mRNA (+)"="#F4A582FF","mRNA (-)"="#c9e3f6FF",
                                  "sncRNA (+)"="#D6604DFF","sncRNA (-)"="#56bdfcFF",
                                  "lncRNA (+)"="#e37b88FF","lncRNA (-)"="#53a4f5FF")) + # colors
    coord_cartesian(ylim = ylimit) # Zoom into an specific area without removing data points
  
  ggsave(paste(title,".png",sep = ""),path = VIOLIN_PLOTS_DIR, scale = 3, width = width, height = height, units = "px", bg = "white", dpi = 600)
  return(p)
}

# Create plots for:
#################
# Conservation. #
df_long_conservation_paper <- df_long_all[df_long_all$feature=="PhyloP-mammals" 
                                          | df_long_all$feature=="PhyloP-vertebrates" 
                                          | df_long_all$feature=="GERP-mammals" 
                                          | df_long_all$feature=="GERP-vertebrates",]
create_violin_plot(df_long_conservation_paper, c(-2, 6), "conservation_review")
create_violin_plot(df_long_conservation_paper, c(-2, 8), "conservation_review2")


###############
# Expression. #
df_long_expression_paper <- df_long_all[df_long_all$feature=="Tissue RPKM" 
                                        | df_long_all$feature=="Primary cell RPKM", ]
create_violin_plot(df_long_expression_paper, c(-1, 500), "expression_review_test")

create_violin_plot(df_long_expression_paper, c(-1,500), "expression_scale_width_adjust_1000_bw_nrd0", violin_scale = "width", violin_adjust = 1000, bw = "nrd0")
create_violin_plot(df_long_expression_paper, c(-1,500), "expression_scale_width_adjust_1000_bw_nrd", violin_scale = "width", violin_adjust = 1000, bw = "nrd")
create_violin_plot(df_long_expression_paper, c(-1,500), "expression_scale_width_adjust_1000_bw_uvc", violin_scale = "width", violin_adjust = 1000, bw = "ucv")
create_violin_plot(df_long_expression_paper, c(-1,500), "expression_scale_width_adjust_1000_bw_bvc", violin_scale = "width", violin_adjust = 1000, bw = "bcv")
create_violin_plot(df_long_expression_paper, c(-1,50), "expression_scale_width_adjust_1000_bw_SJ", violin_scale = "width", violin_adjust = 500, bw = "SJ", bounds = c(-1,500))

create_violin_plot(df_long_expression_paper, c(-1,50), "expression_scale_width_adjust_1_bw_nrd0", violin_scale = "width", violin_adjust = 1/2, bw = "nrd0")
create_violin_plot(df_long_expression_paper, c(-1,50), "expression_scale_width_adjust_0.5_bw_nrd", violin_scale = "width", violin_adjust = 1/2, bw = "nrd")
create_violin_plot(df_long_expression_paper, c(-1,500), "expression_scale_width_adjust_0.5_bw_uvc", violin_scale = "width", violin_adjust = 1/2, bw = "ucv")
create_violin_plot(df_long_expression_paper, c(-1,500), "expression_scale_width_adjust_0.5_bw_bvc", violin_scale = "width", violin_adjust = 1/2, bw = "bcv")
create_violin_plot(df_long_expression_paper, c(-1,500), "expression_scale_width_adjust_0.25_bw_SJ", violin_scale = "width", violin_adjust = 1/4, bw = "SJ", bounds = c(-1,700))

create_violin_plot(df_long_expression_paper, c(-1,500), "expression_scale_width_adjust_0.5_bw_nrd0_jitter", violin_scale = "width", violin_adjust = 1/2, bw = "nrd0")

# Plot separating by gene type:
df_long_expression_mRNA_paper <- df_long_all[(df_long_all$feature=="Tissue RPKM" 
                                        | df_long_all$feature=="Primary cell RPKM")
                                        & 
                                        (df_long_all$`Gene type`=="mRNA (+)"
                                         | df_long_all$`Gene type`=="mRNA (-)"), ]

create_violin_plot(df_long_expression_mRNA_paper, c(-1,175), "exp_mRNA_sc_width_adj_0.5_bw_nrd", violin_scale = "width", violin_adjust = 0.5, bw = "nrd")

df_long_expression_sncRNA_paper <- df_long_all[(df_long_all$feature=="Tissue RPKM" 
                                              | df_long_all$feature=="Primary cell RPKM")
                                             & 
                                               (df_long_all$`Gene type`=="sncRNA (+)"
                                                | df_long_all$`Gene type`=="sncRNA (-)"), ]

create_violin_plot(df_long_expression_sncRNA_paper, c(-1,1200), "exp_sncRNA_sc_width_adj_0.5_bw_nrd", violin_scale = "width", violin_adjust = 0.5, bw = "nrd")

df_long_expression_lncRNA_paper <- df_long_all[(df_long_all$feature=="Tissue RPKM" 
                                                | df_long_all$feature=="Primary cell RPKM")
                                               & 
                                                 (df_long_all$`Gene type`=="lncRNA (+)"
                                                  | df_long_all$`Gene type`=="lncRNA (-)"), ]

create_violin_plot(df_long_expression_lncRNA_paper, c(-1,25), "exp_lncRNA_sc_width_adj_0.5_bw_nrd", violin_scale = "width", violin_adjust = 0.5, bw = "nrd")


##############
# Intrinsic. #
df_long_intrinsic1 <- df_long_all[df_long_all$feature=="GT"
                                  | df_long_all$feature=="GA"
                                  | df_long_all$feature=="CpG"
                                  | df_long_all$feature=="TA"
                                  | df_long_all$feature=="GC%"
                                  | df_long_all$feature=="Low Complexity", ]
custom_order <- c("GC%", "GA", "CpG", "TA", "GT", "Low Complexity")
df_long_intrinsic1$feature <- factor(df_long_intrinsic1$feature, levels = custom_order)
df_long_intrinsic1_ordered <- df_long_intrinsic1 %>%
  arrange(feature)

create_violin_plot(df_long_intrinsic1_ordered, c(-4,3), "intrinsic1", height=1440)
create_violin_plot(df_long_intrinsic1_ordered, c(0,1), "intrinsic2_LowComplexity", height = 1440)



###############
# Epigenetic. #
df_long_epigen_1_paper <- df_long_all[df_long_all$feature=="H3K9ac" 
                                      | df_long_all$feature=="Chromatin Accessibility"
                                      | df_long_all$feature=="H3K79me2"
                                      | df_long_all$feature=="H3K79me1"
                                      | df_long_all$feature=="Methylome", ]
custom_order <- c("H3K79me2", "H3K79me1", "H3K9ac", "Chromatin Accessibility", "Methylome")
df_long_epigen_1_paper$feature <- factor(df_long_epigen_1_paper$feature, levels = custom_order)
create_violin_plot(df_long_epigen_1_paper, c(-5,15), "epigenetic_review")
create_violin_plot(df_long_epigen_1_paper, c(-10,40), "epigenetic_review2")
create_violin_plot(df_long_epigen_1_paper, c(-3,5), "epigenetic_review3")

#OR
df_long_epigen_1_paper <- df_long_all[df_long_all$feature=="H3K27ac" 
                                      | df_long_all$feature=="H3K36me3" 
                                      | df_long_all$feature=="H3K79me2", ]
create_violin_plot(df_long_epigen_1_paper, c(NA,NA), "epigenetic_maxScaledSignal")
create_violin_plot(df_long_epigen_1_paper, c(0,40), "epigenetic_maxScaledSignal2")
create_violin_plot(df_long_epigen_1_paper, c(-3,15), "epigenetic_maxSacaledSignal3")


################
# Association. #
df_long_repeat_assoc_coding <- df_long_all[df_long_all$feature=="Repeat free" 
                                           | df_long_all$feature=="Copies", ]
custom_order <- c("Repeat free", "Copies")
df_long_repeat_assoc_coding$feature <- factor(df_long_repeat_assoc_coding$feature, levels = custom_order)
df_long_repeat_assoc_coding_ordered <- df_long_repeat_assoc_coding %>%
  arrange(feature)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2,8), "association_scale_width_adjust_1", violin_scale = "width", violin_adjust = 1)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_1", violin_scale = "area", violin_adjust = 1)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_count_adjust_1", violin_scale = "count", violin_adjust = 1)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust-0.5", violin_scale = "area", violin_adjust = 1/2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_2", violin_scale = "area", violin_adjust = 2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_0.5", violin_scale = "width", violin_adjust = 1/2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_2", violin_scale = "width", violin_adjust = 2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_count_adjust_0.5", violin_scale = "count", violin_adjust = 1/2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_count_adjust_2", violin_scale = "count", violin_adjust = 2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(NA,NA), "association_scale_area_adjust_1_bounds_null", violin_scale = "area", violin_adjust = 1)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_1_kernel_gauss", violin_scale = "area", violin_adjust = 1, kernel = "gaussian")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_1_kernel_rectangular", violin_scale = "area", violin_adjust = 1, kernel = "rectangular")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_1_kernel_triangular", violin_scale = "area", violin_adjust = 1, kernel = "triangular")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_1_kernel_epanechnikov", violin_scale = "area", violin_adjust = 1, kernel = "epanechnikov")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_1_kernel_biweight", violin_scale = "area", violin_adjust = 1, kernel = "biweight")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_1_kernel_cosine", violin_scale = "area", violin_adjust = 1, kernel = "cosine")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2,8), "association_scale_area_adjust_1_bw_nrd0", violin_scale = "area", violin_adjust = 1, bw = "nrd0") # Good one
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_1_bw_nrd", violin_scale = "area", violin_adjust = 1, bw = "nrd")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_1_bw_ucv", violin_scale = "area", violin_adjust = 1, bw = "ucv")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_1_bw_bcv", violin_scale = "area", violin_adjust = 1, bw = "bcv")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_area_adjust_1_bw_SJ", violin_scale = "area", violin_adjust = 1, bw = "SJ")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_1_bw_nrd0", violin_scale = "width", violin_adjust = 1, bw = "nrd0")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_1_bw_nrd", violin_scale = "width", violin_adjust = 1, bw = "nrd")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_1_bw_uvc", violin_scale = "width", violin_adjust = 1, bw = "ucv")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_1_bw_bvc", violin_scale = "width", violin_adjust = 1, bw = "bcv")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_1_bw_SJ", violin_scale = "width", violin_adjust = 1, bw = "SJ")

create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_1000_bw_nrd0", violin_scale = "width", violin_adjust = 1000, bw = "nrd0")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_1000_bw_nrd", violin_scale = "width", violin_adjust = 1000, bw = "nrd")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_1000_bw_uvc", violin_scale = "width", violin_adjust = 1000, bw = "ucv")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_1000_bw_bvc", violin_scale = "width", violin_adjust = 1000, bw = "bcv")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-2.5,8), "association_scale_width_adjust_1000_bw_SJ", violin_scale = "width", violin_adjust = 1000, bw = "SJ")

create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-100,5000), "association_scale_width_adjust_1000_bw_nrd0", violin_scale = "width", violin_adjust = 1000, bw = "nrd0")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-100,5000), "association_scale_width_adjust_1000_bw_nrd", violin_scale = "width", violin_adjust = 1000, bw = "nrd")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-100,5000), "association_scale_width_adjust_1000_bw_uvc", violin_scale = "width", violin_adjust = 1000, bw = "ucv")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-100,5000), "association_scale_width_adjust_1000_bw_bvc", violin_scale = "width", violin_adjust = 1000, bw = "bcv")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-100,4000), "association_scale_width_adjust_1000_bw_SJ", violin_scale = "width", violin_adjust = 1000, bw = "SJ")
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(-100,5000), "association_scale_width_adjust_500_bw_SJ", violin_scale = "width", violin_adjust = 500, bw = "SJ")



create_violin_plot(df_long_repeat_assoc_coding_ordered, c(0,1000), "association_review2", violin_adjust = 1/2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(0,1000), "association_review2", violin_adjust = 1/2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(0,1000), "association_review2", violin_adjust = 1/2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(0,1000), "association_review2", violin_adjust = 1/2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(0,1000), "association_review2", violin_adjust = 1/2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(0,1000), "association_review2", violin_adjust = 1/2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(0,1000), "association_review2", violin_adjust = 1/2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(0,1000), "association_review2", violin_adjust = 1/2)
create_violin_plot(df_long_repeat_assoc_coding_ordered, c(0,1000), "association_review2", violin_adjust = 1/2)


help("geom_violin")
# Structure.
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

# Population Variation.
df_long_pop_var <- df_long_all[df_long_all$feature=="MAF"
                               | df_long_all$feature=="SNPs", ]
create_violin_plot(df_long_pop_var, c(-2,10), "popvar_review")
