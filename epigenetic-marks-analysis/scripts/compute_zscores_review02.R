###############
# Script to Compute z-scores for all gene functionality features separated by Gene Type.
# The aim is to explore the resulting z-scores and determine the reason for some 
# extremely big values obtained for some features.
# Can be executed from the command line using:
#
#  Rscript compute_zscores_review02.R ../data/features/gene_functionality_features_latest1000all_extra.csv --subset=100
###############
options(repos = c(CRAN = "https://cran.r-project.org"))
library(dplyr)
library(ggplot2)
library(tidyr)
install.packages("patchwork", dependencies = TRUE)
library(patchwork)
library(scales)

# Function to compute robust z-scores
compute_robust_zscores_optimized <- function(positive_matrix, negative_matrix, verbose = FALSE) {
  # Use NEGATIVE CONTROLS median absolute deviations (MAD) and medians to compute robust z-scores
  mad_values <- apply(negative_matrix, 2, mad, na.rm = TRUE)  # Calculate MAD for each column in the negative control set
  median_values <- apply(negative_matrix, 2, median, na.rm = TRUE)  # Calculate median for each column in the negative control set
  meanAD_values <- apply(negative_matrix, 2, function(col) mean(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE))
  
  # Handle case where MAD is 0 or very small by using meanAD instead
  mad_values <- apply(negative_matrix, 2, function(col) {
    mad_value <- mad(col, na.rm = TRUE)
    if (mad_value <= 1e-6) {  # If MAD is 0 or very small, calculate mean absolute deviation (meanAD)
      meanAD <- mean(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE)
      warning(sprintf("Column %s with very small MAD of %f encountered: Using meanAD value of %f.\n", col, mad_value, meanAD))
      return(1.2533 * meanAD)  # Use meanAD if it is greater than the threshold
    } else {
      return(mad_value)  # Use MAD if it is sufficient
    }
  })
  
  # Calculate z-scores: subtract the median and divide by MAD for each feature
  zscores <- sweep(positive_matrix, 2, median_values, "-")  # Center each column by subtracting the median
  zscores <- sweep(zscores, 2, mad_values, "/")  # Scale each column by dividing by MAD
  
  print("\n\nMedian values:\n")
  print(median_values)
  print("MAD values:\n")
  print(mad_values)
  print("meanAD values:\n")
  print(meanAD_values)
  
  zscores
}

# Main function
main <- function() {
  # Capture command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Please provide the path to the data file as an argument.")  # Ensure a file path is provided
  }
  names(feature_matrix)[1] <- "row"
  # Load and preprocess data
  data_file <- args[1]  # Get the file path from the command line arguments
  data <- read.csv(data_file, header = TRUE, check.names = TRUE)# Load the data from the specified CSV file
  data <- feature_matrix %>% dplyr::select(-Dataset,-ID,-Functional,-Chromosome,-Start,-End,-Sequence)
  
  # Check if subset option is provided
  #subset_value <- as.numeric(gsub("--subset=", "", args[grep("--subset=", args)]))
  #if (is.na(subset_value)) {
    subset_value <- nrow(data)  # Use nrow to indicate no subsetting by default
    subset_value <- 100
  #}
  
  # Convert to numeric
  data_numeric <- data %>% 
    #dplyr::select(-Dataset) %>%
    sapply(function(feature) as.numeric(as.character(feature))) %>%
    as.data.frame()
  summary(data_numeric)
  # Restore Dataset column
  data_numeric$Dataset <- feature_matrix$Dataset
  unique(data_numeric$Dataset)
  table(data_numeric$Dataset)
  #Protein coding
  # Separate the data into positive and negative feature matrices
  protein_positive_feature_matrix <- data_numeric %>% 
    filter(Dataset == "protein-coding-exon2" | Dataset == "protein-coding-exon3") #%>% 
    #dplyr::select(-Dataset)
  protein_negative_feature_matrix <- data_numeric %>% 
    filter(Dataset == "protein-exon2-negative-control" | Dataset == "protein-exon3-negative-control")# %>% 
    #dplyr::select(-Dataset)
  
  # Select random samples from positive and negative datasets
  #set.seed(197)  # Set seed for reproducibility
  protein_positive_feature_matrix <- protein_positive_feature_matrix %>% sample_n(min(2 * subset_value, nrow(protein_positive_feature_matrix)))  # Subset rows from positive data
  protein_negative_feature_matrix <- protein_negative_feature_matrix %>% sample_n(min(2 * 10 * subset_value, nrow(protein_negative_feature_matrix)))  # Subset rows from negative data

  print(summary(protein_positive_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  print(summary(protein_negative_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))  
  #print(apply(protein_positive_feature_matrix, 2, sd, na.rm = TRUE))
  #apply(protein_negative_feature_matrix, 2, sd, na.rm = TRUE)
  #print(sd(protein_negative_feature_matrix$RPKM_tissue, na.rm = TRUE))
  # Compute z-scores
  protein_functional_z_scores <- compute_robust_zscores_optimized(protein_positive_feature_matrix %>% dplyr::select(-Dataset), 
                                                                  protein_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                                  verbose)
  protein_negative_z_scores <- compute_robust_zscores_optimized(protein_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                                protein_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                                verbose)
  
  print(summary(protein_functional_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  print(count(protein_positive_feature_matrix))
  print(summary(protein_negative_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  print(summary(protein_negative_feature_matrix))
  #print(sd(protein_functional_z_scores$RPKM_tissue, na.rm = TRUE))
  #print(sd(protein_negative_z_scores$RPKM_tissue, na.rm = TRUE))
  
  #Restore Dataset
  protein_functional_z_scores$Dataset <- (protein_positive_feature_matrix %>% 
    filter(Dataset == "protein-coding-exon2" | Dataset == "protein-coding-exon3") %>% 
      dplyr::select(Dataset))$Dataset
  protein_negative_z_scores$Dataset <- (protein_negative_feature_matrix %>% 
    filter(Dataset == "protein-exon2-negative-control" | Dataset == "protein-exon3-negative-control") %>% 
      dplyr::select(Dataset))$Dataset
  #Join in one dataframe
  protein_zscores_all <- rbind(protein_functional_z_scores,protein_negative_z_scores)
  # Save z-scores to CSV
  z_scores_file <- "pasted_histone_outputs/mrna_z_scores.csv"
  write.csv(protein_functional_z_scores, z_scores_file, row.names = FALSE)
  
  z_scores_file <- "pasted_histone_outputs/mrna_negative_control_z_scores.csv"
  write.csv(protein_negative_z_scores, z_scores_file, row.names = FALSE)
  
  # Save raw data subset
  subset_file <- "pasted_histone_outputs/mrna_features.csv"
  write.csv(protein_positive_feature_matrix, subset_file, row.names = FALSE)
  
  subset_file <- "pasted_histone_outputs/mrna_negative_control_features.csv"
  write.csv(protein_negative_feature_matrix, subset_file, row.names = FALSE)
  
  #lncRNA
  # Separate the data into positive and negative feature matrices
  lncrna_positive_feature_matrix <- data_numeric %>% 
    filter(Dataset == "lncrna-exon1" | Dataset == "lncrna-exon2") #%>% 
    #dplyr::select(-Dataset)
  lncrna_negative_feature_matrix <- data_numeric %>% 
    filter(Dataset == "lncrna-exon1-negative-control" | Dataset == "lncrna-exon2-negative-control") #%>% 
    #dplyr::select(-Dataset)
  
  # Select random samples from positive and negative datasets
  #set.seed(197)  # Set seed for reproducibility
  lncrna_positive_feature_matrix <- lncrna_positive_feature_matrix %>% sample_n(min(2 * subset_value, nrow(lncrna_positive_feature_matrix)))  # Subset rows from positive data
  lncrna_negative_feature_matrix <- lncrna_negative_feature_matrix %>% sample_n(min(2 * 10 * subset_value, nrow(lncrna_negative_feature_matrix)))  # Subset rows from negative data
  
  # Compute z-scores
  lncrna_functional_z_scores <- compute_robust_zscores_optimized(lncrna_positive_feature_matrix %>% dplyr::select(-Dataset), 
                                                                 lncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                                 verbose)
  lncrna_negative_z_scores <- compute_robust_zscores_optimized(lncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                               lncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                               verbose)
  
  print(summary(lncrna_positive_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  count(lncrna_positive_feature_matrix)
  print(summary(lncrna_negative_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  count(lncrna_negative_feature_matrix)
  
  print(summary(lncrna_functional_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  print(summary(lncrna_negative_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  
  #Restore Dataset
  lncrna_functional_z_scores$Dataset <- (lncrna_positive_feature_matrix %>% 
    filter(Dataset == "lncrna-exon1" | Dataset == "lncrna-exon2") %>% 
    dplyr::select(Dataset))$Dataset
  lncrna_negative_z_scores$Dataset <- (lncrna_negative_feature_matrix %>% 
    filter(Dataset == "lncrna-exon1-negative-control" | Dataset == "lncrna-exon2-negative-control") %>% 
    dplyr::select(Dataset))$Dataset
  
  # Save z-scores to CSV
  z_scores_file <- "pasted_histone_outputs/lncrna_z_scores.csv"
  write.csv(lncrna_functional_z_scores, z_scores_file, row.names = FALSE)
  
  z_scores_file <- "pasted_histone_outputs/lncrna_negative_control_z_scores.csv"
  write.csv(lncrna_negative_z_scores, z_scores_file, row.names = FALSE)
  
  # Save raw subset feature to CSV
  subset_file <- "pasted_histone_outputs/lncrna_features.csv"
  write.csv(lncrna_positive_feature_matrix, subset_file, row.names = FALSE)
  
  subset_file <- "pasted_histone_outputs/lncrna_negative_control_features.csv"
  write.csv(lncrna_negative_feature_matrix, subset_file, row.names = FALSE)
  
  #sncRNA
  # Separate the data into positive and negative feature matrices
  sncrna_positive_feature_matrix <- data_numeric %>% 
    filter(Dataset == "short-ncrna" | Dataset == "short-ncrna") #%>% 
    #dplyr::select(-Dataset)
  sncrna_negative_feature_matrix <- data_numeric %>% 
    filter(Dataset == "short-ncrna-negative-control" | Dataset == "short-ncrna-negative-control") #%>% 
    #dplyr::select(-Dataset)
  
  # Select random samples from positive and negative datasets
  #set.seed(197)  # Set seed for reproducibility
  sncrna_positive_feature_matrix <- sncrna_positive_feature_matrix %>% sample_n(min(subset_value, nrow(sncrna_positive_feature_matrix)))  # Subset rows from positive data
  sncrna_negative_feature_matrix <- sncrna_negative_feature_matrix %>% sample_n(min(10 * subset_value, nrow(sncrna_negative_feature_matrix)))  # Subset rows from negative data
  
  # Compute z-scores
  sncrna_functional_z_scores <- compute_robust_zscores_optimized(sncrna_positive_feature_matrix %>% dplyr::select(-Dataset), 
                                                                 sncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                                 verbose)
  sncrna_negative_z_scores <- compute_robust_zscores_optimized(sncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                               sncrna_negative_feature_matrix %>% dplyr::select(-Dataset), 
                                                               verbose)
  
  print(summary(sncrna_positive_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  count(sncrna_positive_feature_matrix)
  print(summary(sncrna_negative_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  count(sncrna_negative_feature_matrix)
  
  print(summary(sncrna_functional_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  print(summary(sncrna_negative_z_scores %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  
  
  #Restore Dataset
  sncrna_functional_z_scores$Dataset <- (sncrna_positive_feature_matrix %>% 
    filter(Dataset ==  "short-ncrna" | Dataset == "short-ncrna") %>% 
    dplyr::select(Dataset))$Dataset
  sncrna_negative_z_scores$Dataset <- (sncrna_negative_feature_matrix %>% 
    filter(Dataset == "short-ncrna-negative-control" | Dataset == "short-ncrna-negative-control") %>% 
    dplyr::select(Dataset))$Dataset
  # Save z-scores to CSV
  z_scores_file <- "pasted_histone_outputs/sncrna_z_scores.csv"
  write.csv(sncrna_functional_z_scores, z_scores_file, row.names = FALSE)
  
  z_scores_file <- "pasted_histone_outputs/sncrna_negative_control_z_scores.csv"
  write.csv(sncrna_negative_z_scores, z_scores_file, row.names = FALSE)
  
  # Save raw subset to CSV
  subset_file <- "pasted_histone_outputs/sncrna_features.csv"
  write.csv(sncrna_positive_feature_matrix, subset_file, row.names = FALSE)
  
  subset_file <- "pasted_histone_outputs/sncrna_negative_control_features.csv"
  write.csv(sncrna_negative_feature_matrix, subset_file, row.names = FALSE)
  
  # Join all z-scores in one df:
  zscores_all <- rbind(protein_functional_z_scores,lncrna_functional_z_scores,sncrna_functional_z_scores,
                       protein_negative_z_scores,lncrna_negative_z_scores,sncrna_negative_z_scores)
  
  #############
  # Compute ks-stats for z-scores
  #mrna
  subset_m <- protein_functional_z_scores %>% dplyr::select(-Dataset)
  subset_nc_m <- protein_negative_z_scores %>% dplyr::select(-Dataset)
  
  # Run the K-S tests
  ks_results_m <- run_ks_tests(subset_nc_m, subset_m)
  
  #lncrna
  subset_l <- lncrna_functional_z_scores %>% dplyr::select(-Dataset)
  subset_nc_l <- lncrna_negative_z_scores %>% dplyr::select(-Dataset)
  
  # Run the K-S tests
  ks_results_l <- run_ks_tests(subset_nc_l, subset_l)
  
  #sncrna
  subset_s <- sncrna_functional_z_scores %>% dplyr::select(-Dataset)
  subset_nc_s <- sncrna_negative_z_scores %>% dplyr::select(-Dataset)
  
  # Run the K-S tests
  ks_results_s <- run_ks_tests(subset_nc_s, subset_s)
  
  ###############
  # Compute ks-stats for raw features #
  #mrna
  subset_m <- protein_positive_feature_matrix %>% dplyr::select(-Dataset)
  subset_nc_m <- protein_negative_feature_matrix %>% dplyr::select(-Dataset)
  
  # Run the K-S tests
  ks_results_m <- run_ks_tests(subset_nc_m, subset_m)
  
  print(ks_results_m)
  write.csv(ks_results_m,"mrna_ks_stat.csv")
  
  #lncrna
  subset_l <- lncrna_positive_feature_matrix %>% dplyr::select(-Dataset)
  subset_nc_l <- lncrna_negative_feature_matrix %>% dplyr::select(-Dataset)
  
  # Run the K-S tests
  ks_results_l <- run_ks_tests(subset_nc_l, subset_l)
  write.csv(ks_results_m,"lncrna_ks_stat.csv")
  
  #sncrna
  subset_s <- sncrna_positive_feature_matrix %>% dplyr::select(-Dataset)
  subset_nc_s <- sncrna_negative_feature_matrix %>% dplyr::select(-Dataset)
  
  # Run the K-S tests
  ks_results_s <- run_ks_tests(subset_nc_s, subset_s)
  write.csv(ks_results_m,"sncrna_ks_stat.csv")
  
  
  # Generate plots to visualize distributions:
  # Z-scores
  p1 <- create_distribution_plots(protein_functional_z_scores %>% dplyr::select(-Dataset),
                            "zscores_mrna", color = "#F4A582FF", 
                            output_dir = "plots/zscores/new_subset", title = "mRNA(+)", tag = c('A','B','C'))
  p2 <- create_distribution_plots(protein_negative_z_scores %>% dplyr::select(-Dataset), 
                            "zscores_mrna_negative", "#c9e3f6FF", output_dir = "plots/zscores/new_subset",
                            type = "Negative Control", title = "mRNA(-)", tag = c('D','E','F'))
  
  p3 <- create_distribution_plots(sncrna_functional_z_scores, 
                            "zscores_sncRNA", color = "#D6604DFF", 
                            output_dir = "plots/zscores/new", title = "sncRNA(+)", tag = c('A','B','C'))
  p4 <- create_distribution_plots(sncrna_negative_z_scores, 
                            "zscores_sncRNA_negative", color = "#56bdfcFF", 
                            output_dir = "plots/zscores/new",
                            type = "Negative Control", title = "sncRNA(-)", tag = c('D','E','F'))
  
  pp1 <- p1[["H3K27ac_MaxScaledSignal"]]
  pp2 <- p2[["H3K27ac_MaxScaledSignal"]]
  
  final_plot_pyhloP <- (pp1 / pp2) +
    plot_annotation(
      title = "PhyloP mammals"
    ) &
    theme(
      plot.title = element_text(size = 84,
                                hjust = 0.5)  # center the title
    )
  
  final_plot_copies <- (p1[["Copy.number"]] / p2[["Copy.number"]]) +
    plot_annotation(
      title = "Copies"
    ) &
    theme(
      plot.title = element_text(size = 84,
                                hjust = 0.5)  # center the title
    )
  
  final_plot_RPKM <- (p3[["RPKM_primary.cell"]] / p4[["RPKM_primary.cell"]]) +
    plot_annotation(
      title = "Primary cell RPKM"
    ) &
    theme(
      plot.title = element_text(size = 84,
                                hjust = 0.5)  # center the title
    )
  
  ggsave(file.path("plots/zscores/new", paste0("mRNA_phyloP.max241w_combined_log.png")), final_plot_pyhloP, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
  ggsave(file.path("plots/zscores/new", paste0("mRNA_Copy.number_combined_log.png")), final_plot_copies, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
  ggsave(file.path("plots/zscores/new", paste0("sncRNA_RPKM_primary.cell_combined_log.png")), final_plot_RPKM, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
}

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
      scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
                         labels = scales::label_number()
                         ) +
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
      scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
                         labels = scales::label_number()
      ) +
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
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
                         labels = scales::label_number()
      ) +
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

# Run the main function
main()
