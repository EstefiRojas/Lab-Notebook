# Create violin plots from gene functionality z-score data.
# Load necessary libraries.
library(ggplot2)
library(tidyr)
library(patchwork)
library(scales)

source("scripts/config.R")
source("scripts/load_gene_functionality_zscores.R")
source("scripts/utils.R")
zscores_all <- load_gene_functionality_zscores()

PLOT_SELECTED_FEATURES <- c("Random","GC_percentage",
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
                            "H3K9ac_MaxScaledSignal","H3K79me2_MaxScaledSignal","chrm_acc_MaxScaledSignal",
                            "methylome")

# Function to create distribution plots
create_distribution_plots <- function(data, gene_type, color, output_dir = "plots", axis_label = "Robust Z-score", type = "Functional", title, tag) {
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Prepare data for plotting
  long_data <- data %>%
    tidyr::pivot_longer(everything(), names_to = "Feature", values_to = "Value") %>%
    mutate(Type = type)
  
  # List of features to plot
  features <- unique(long_data$Feature)
  plots_list <- list()
  
  # Create plots for each feature
  for (feature in features) {
    feature_data <- long_data %>% filter(Feature == feature)
    # Remove NAs
    feature_data <- feature_data[!is.na(feature_data$Value), ]
    # Ensure column is numeric
    feature_data$Value <- as.numeric(feature_data$Value)
    
    # 1. Histogram
    p1 <- ggplot(feature_data, aes(x = Value, fill = Type, colour = Type)) +
      geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
      labs(x = axis_label, y = "Count") +
      scale_fill_manual(values = c(color)) +
      scale_color_manual(values = c(color)) +
      #scale_x_log10(breaks = scales::log_breaks(n = 10)) +
      #scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
      #                   breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
      #                   labels = scales::label_number()
      #) +
      theme_minimal(base_size = 84) +
      theme(axis.text.x = element_text(size = 28, angle = 45),
            axis.text.y = element_text(size = 28),  # Increase y-axis text size
            axis.title.x = element_text(size = 44),  # Increase x-axis title size
            axis.title.y = element_text(size = 44),  # Increase y-axis title size
            legend.position = "none"
      )
    
    #ggsave(file.path(output_dir, paste0(gene_type, "_", feature, "_histogram.png")), p1, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
    
    # 2. Density plot
    p2 <- ggplot(feature_data, aes(x = Value, fill = Type, colour = Type)) +
      geom_density(alpha = 0.6) +
      theme_minimal(base_size = 84) +
      labs(x = axis_label, y = "Density") +
      scale_fill_manual(values = c(color)) +
      scale_color_manual(values = c(color)) +
      #scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
      #                   breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
      #                   labels = scales::label_number()
      #) +
      theme(axis.text.x = element_text(size = 28, angle = 45),
            axis.text.y = element_text(size = 28),  # Increase y-axis text size
            axis.title.x = element_text(size = 44),  # Increase x-axis title size
            axis.title.y = element_text(size = 44),  # Increase y-axis title size
            legend.position = "none"#,
            #plot.title = element_text(size = 64, face = "italic", hjust = 0.5)  # center the title
      )
    
    #ggsave(file.path(output_dir, paste0(gene_type, "_", feature, "_density.png")), p2, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
    
    # 3. Boxplot
    p3 <- ggplot(feature_data, aes(x = feature, y = Value, fill = Type)) +
      geom_boxplot(outlier.colour = color, outlier.size = 5) +
      theme_minimal(base_size = 84) +
      labs(y = axis_label) +
      #scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
      #                   breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
      #                   labels = scales::label_number()
      #) +
      theme(axis.text.x = element_text(size = 0, angle = 45),
            axis.text.y = element_text(size = 28),  # Increase y-axis text size
            axis.title.x = element_text(size = 0),  # Increase x-axis title size
            axis.title.y = element_text(size = 44),  # Increase y-axis title size
            legend.position = "none"
      ) +
      scale_fill_manual(values = c(color))
    
    #ggsave(file.path(output_dir, paste0(gene_type, "_", feature, "_boxplot.png")), p3, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
    
    # 4. Join all plots in one
    combined_plot <- (p1 + p2 + p3)
    
    combined_plot_annotated <- wrap_elements(combined_plot + 
                                               plot_annotation(title = title,
                                                               tag_levels = list(tag)) &
                                               theme(plot.title = element_text(size = 64,
                                                                               face = "italic",
                                                                               hjust = 0.5),
                                                     strip.text = element_text(margin = margin(0,0,30,0)),
                                                     plot.tag.position = c(0, 1),
                                                     plot.tag = element_text(size = 38, face = "bold", hjust = 0, vjust = 0, margin = margin(40,40,40,40))))
    
    ggsave(file.path(output_dir, paste0(gene_type, "_", feature, "_combined.png")), combined_plot_annotated, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
    plots_list[[feature]] <- combined_plot_annotated
  }
  return(plots_list)
}


# Generate plots to visualize distributions:
# Methylome
unique(zscores_all$Dataset)
p1 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_FEATURES)),
                                "zscores_mrna", color = "#F4A582FF", 
                                output_dir = "results/distributionPlots/", title = "mRNA(+)", tag = c('A','B','C'))
p2 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("protein-exon2-negative-control", "protein-exon3-negative-control")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_FEATURES)),
                                "zscores_mrna_negative", "#c9e3f6FF", 
                                output_dir = "results/distributionPlots/",
                                type = "Negative Control", title = "mRNA(-)", tag = c('D','E','F'))

p3 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("short-ncrna")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_FEATURES)),
                                "zscores_sncRNA", color = "#D6604DFF", 
                                output_dir = "results/distributionPlots/", title = "sncRNA(+)", tag = c('A','B','C'))
p4 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("short-ncrna-negative-control")) %>% 
                                  dplyr::select(all_of(c("methylome"))),
                                "zscores_sncRNA_negative", color = "#56bdfcFF", 
                                output_dir = "results/distributionPlots/",
                                type = "Negative Control", title = "sncRNA(-)", tag = c('D','E','F'))

p5 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("lncrna-exon2", "lncrna-exon3")) %>% 
                                  dplyr::select(all_of(c("methylome"))),
                                "zscores_lncRNA", color = "#e37b88FF", 
                                output_dir = "results/distributionPlots/", title = "lncRNA(+)", tag = c('A','B','C'))
p6 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("lncrna-exon2-negative-control", "lncrna-exon3-negative-control")) %>% 
                                  dplyr::select(all_of(c("methylome"))),
                                "zscores_lncRNA_negative", color = "#53a4f5FF", 
                                output_dir = "results/distributionPlots/",
                                type = "Negative Control", title = "lncRNA(-)", tag = c('D','E','F'))

pp1 <- p1[["methylome"]]
pp2 <- p2[["methylome"]]

final_plot_methylome_mrna <- (pp1 / pp2) +
  plot_annotation(
    title = "Methylome"
  ) &
  theme(
    plot.title = element_text(size = 84,
                              hjust = 0.5)  # center the title
  )

final_plot_methylome_sncrna <- (p3[["methylome"]] / p4[["methylome"]]) +
  plot_annotation(
    title = "Methylome"
  ) &
  theme(
    plot.title = element_text(size = 84,
                              hjust = 0.5)  # center the title
  )

final_plot_methylome_lncrna <- (p5[["methylome"]] / p6[["methylome"]]) +
  plot_annotation(
    title = "Methylome"
  ) &
  theme(
    plot.title = element_text(size = 84,
                              hjust = 0.5)  # center the title
  )

ggsave(file.path("results/distributionPlots/", paste0("mRNA_methylome_combined_log.png")), final_plot_methylome_mrna, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
ggsave(file.path("results/distributionPlots/", paste0("sncRNA_methylome_combined_log.png")), final_plot_methylome_sncrna, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
ggsave(file.path("results/distributionPlots/", paste0("lncRNA_methylome_combined_log.png")), final_plot_methylome_lncrna, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)



