# Effect on Distance Analysis, delta z-scores based on distance to gene
# Load libraries:
library(dplyr)
library(tidyr)
library(ggplot2)
library(robustbase) # For the Sn scale estimator
install.packages(c("ggpubr", "broom", "AICcmodavg"))
library(ggpubr)
library(broom)
library(AICcmodavg)

##########################
# Define helper functions:
##########################
# Rename features to new paper names
rename_features <- function(dataframe, current_names, new_names) {
  for (i in 1:length(current_names)) {  # Iterate from 1 to 5
    names(dataframe)[ names(dataframe) == current_names[i] ] <- new_names[i]
  }
  return(dataframe)
}

# Subset negative datasets by Distance to Gene:
subset_by_Distance_to_Gene <- function(dataset) {
  list_of_subsets <- list()
  list_of_subsets[['1-5k']] <- dataset[dataset$DistanceGene > 0 & dataset$DistanceGene <= 5000,]
  list_of_subsets[['5k-50k']] <- dataset[dataset$DistanceGene > 5000 & dataset$DistanceGene <= 50000,]
  list_of_subsets[['50k-500k']] <- dataset[dataset$DistanceGene > 50000 & dataset$DistanceGene <= 500000,]
  list_of_subsets[['500k-5M']] <- dataset[dataset$DistanceGene > 500000 & dataset$DistanceGene <= 5000000,]
  return(list_of_subsets)
}

# Function to compute gene coords upstream or downstream
get_gene_coords <- function(dataset) {
  dataset$GeneUpstreamEnd <- pmax(1, dataset$Start - dataset$DistanceGene - 1)
  dataset$GeneDownstreamStart <- dataset$End + dataset$DistanceGene + 1
  return(dataset)
}

# Function to add sequence length as a feature
add_seq_length <- function(dataset, is_positive=TRUE) {
  if(!is_positive){
    dataset_1 <- dataset
    dataset_2 <- dataset
    dataset_exact <- dataset
    dataset_1$seqLength <- nchar(dataset$Sequence) - 1
    dataset_2$seqLength <- nchar(dataset$Sequence) - 2
    dataset_exact$seqLength <- nchar(dataset$Sequence)
    dataset <- rbind(dataset_exact, dataset_1, dataset_2)
  } else
    dataset$seqLength <- nchar(dataset$Sequence)
  return(dataset)
}

# Function to add sequence length as a feature
add_seq_length1 <- function(dataset) {
  dataset$seqLength <- nchar(dataset$Sequence)
  return(dataset)
}

# Function to join positive data to corresponding negative data by chromosome and distance to gene
get_upstream_and_downstream_sets <- function(positive_df, negative_df, selected_features) {
  downstreamSet <- inner_join(negative_df, 
                              positive_df,
                              by = join_by(Chromosome, x$GeneDownstreamStart == y$Start),
                              suffix = c(".x",".y")) %>%
    select(c(paste0(selected_features,".x"), paste0(selected_features, ".y"), "DistanceGene"))
  
  upstreamSet <- inner_join(negative_df, 
                            positive_df,
                            by = join_by(Chromosome, x$GeneUpstreamEnd == y$End),
                            suffix = c(".x",".y")) %>%
    select(c(paste0(selected_features,".x"), paste0(selected_features, ".y"), "DistanceGene"))
  
  superset <- rbind(downstreamSet, upstreamSet)
  return(superset)
}

one_to_one_join <- function(negative_df, positive_df) {
  # Create a unique identifier for each row in positive_df
  positive_df <- positive_df %>%
    mutate(original_row = row_number())
  
  result <- negative_df %>%
    rowwise() %>%
    mutate(match = list(positive_df %>%
                          filter(Chromosome == .data$Chromosome & seqLength == .data$seqLength) %>%
                          slice(1))) %>%
    unnest(match, keep_empty = TRUE, names_sep = ".")
  
  # Remove duplicates from positive side, keeping only the first occurrence
  result <- result %>%
    group_by(original_row) %>%
    slice(1) %>%
    ungroup()
  
  # Rename columns to add suffixes
  neg_cols <- names(negative_df)
  pos_cols <- setdiff(names(result), c(neg_cols, "original_row"))
  
  result <- result %>%
    rename_with(~paste0(., ".neg"), all_of(neg_cols)) %>%
    rename_with(~paste0(., ".pos"), all_of(pos_cols))
  
  # Remove helper column
  result <- result %>% select(-original_row)
  
  return(result)
}

# Function to join positive data to corresponding negative data by chromosome and sequence length
get_sequence_length_sets <- function(positive_df, negative_df) {
  #lengthSet <- one_to_one_join(negative_df, positive_df)
  lengthSet <- left_join(negative_df, 
                          positive_df,
                          by = join_by(Chromosome, seqLength),
                          suffix = c(".neg",".pos"),
                          relationship = "many-to-one",
                          multiple = "any",
                          unmatched = "drop",
                          )
  
  return(lengthSet)
}

# Compute robust z-scores
get_robust_zscores <- function(subsets_df, selected_features, is_negative=FALSE) {
  list_with_zscores <- list()
  list_with_method <- list()
  for (col in selected_features) {
    col_name <- if(is_negative) paste0(col,".neg") else paste0(col,".pos")
    positive_col <- as.numeric(subsets_df[[col_name]])
    negative_col <- as.numeric(subsets_df[[paste0(col,".neg")]])

    mad_value <- mad(negative_col, na.rm = TRUE)
    median_value <- median(negative_col, na.rm = TRUE)
    Sn_value <- Sn(negative_col, na.rm = TRUE)
    Qn_value <- Qn(negative_col, na.rm = TRUE)
    tau_value <- scaleTau2(negative_col, na.rm = TRUE)
    
    mean_value <- mean(negative_col, na.rm = TRUE) 
    absolute_deviations <- abs(negative_col - mean_value)
    meanAD_value <- mean(absolute_deviations, na.rm = TRUE)
    
    if(mad_value != 0) { # MAD method
      list_with_zscores[[col]] <- (positive_col - median_value) / mad_value
      list_with_method[[col]] <- "MAD"
    } else if(Sn_value != 0) { # Sn estimator
      list_with_zscores[[col]] <- (positive_col - median_value) / Sn_value
      list_with_method[[col]] <- "Sn estimator"
    } else if(Qn_value != 0) { # # Qn estimator
      list_with_zscores[[col]] <- (positive_col - median_value) / Qn_value
      list_with_method[[col]] <- "Qn estimator"
    } else if(tau_value != 0) { # tau estimator
      list_with_zscores[[col]] <- (positive_col - median_value) / tau_value
      list_with_method[[col]] <- "tau estimator"
    } else if(meanAD_value != 0) { # meanAD estimator
      list_with_zscores[[col]] <- (positive_col - median_value) / (1.2533 * meanAD_value)
      list_with_method[[col]] <- "meanAD"
    } else { # regular z-score
      list_with_zscores[[col]] <- (positive_col - mean(negative_col)) / sd(negative_col, na.rm = TRUE)
      list_with_method[[col]] <- "regular z-score"
    }
    #list_with_zscores[[col]] <- 1.1926 * (positive_col - median_value) / Sn_value
  }
  return(list(zscores=list_with_zscores, method=list_with_method))
}

# 3. Compute z-scores:
compute_zscores_all_bins <- function(subsets_df, selected_features) {
  zscores_bin1 <- list(positive=get_robust_zscores(subsets_df$`1-5k`, selected_features)$zscores,
                       negative=get_robust_zscores(subsets_df$`1-5k`, selected_features, is_negative = TRUE)$zscores)
  zscores_bin2 <- list(positive=get_robust_zscores(subsets_df$`5k-50k`, selected_features)$zscores,
                        negative=get_robust_zscores(subsets_df$`5k-50k`, selected_features, is_negative = TRUE)$zscores)
  zscores_bin3 <- list(positive=get_robust_zscores(subsets_df$`50k-500k`, selected_features)$zscores,
                        negative=get_robust_zscores(subsets_df$`50k-500k`, selected_features, is_negative = TRUE)$zscores)
  zscores_bin4 <- list(positive=get_robust_zscores(subsets_df$`500k-5M`, selected_features)$zscores,
                        negative=get_robust_zscores(subsets_df$`500k-5M`, selected_features, is_negative = TRUE)$zscores)
  return(list(bin1=zscores_bin1,bin2=zscores_bin2,bin3=zscores_bin3,bin4=zscores_bin4))
}

# Compute deltas
calculate_delta_zscores <- function(zscores_df1, zscores_df2) {
  # Ensure both data frames have the same columns
  common_columns <- intersect(colnames(zscores_df1), colnames(zscores_df2))
  
  if (length(common_columns) == 0) {
    stop("No common features found between the two data frames")
  }
  
  # Calculate delta z-scores
  delta_zscores <- zscores_df2[common_columns] - zscores_df1[common_columns]
  
  return(delta_zscores)
}

# Plot delta z-scores against Distance to Gene bins
generate_violin_plot <- function(delta_zscores, selected_features, title, output_file_name) {
  for(col in selected_features) {
    ggplot(delta_zscores, aes(x = bin, y = !!sym(col), fill = bin)) +
      geom_violin(scale = "width", show.legend = TRUE) +
      geom_boxplot(alpha=0.0, outliers=TRUE, position = position_dodge(width = 0.9), width=0.2) +
      facet_wrap(~ "", scales = "free") + 
      labs(title = paste(title, '-', gsub("_", " ", col)), x = col, y = "Delta Z-score") +
      theme_minimal(base_size = 26) +
      theme(
        axis.text.x = element_text(size = 20, angle = 0),  # Increase x-axis text size and rotate labels
        axis.text.y = element_text(size = 20),  # Increase y-axis text size
        axis.title.x = element_text(size = 0),  # Increase x-axis title size
        axis.title.y = element_text(size = 26),  # Increase y-axis title size
        legend.position = "none",
        legend.title = element_text(size = 20),  # Increase legend title size
        legend.text = element_text(size = 18),  # Increase legend text size
        plot.title = element_text(size = 36, hjust = 0.5)  # Increase plot title size and center it
      ) +
      scale_fill_manual(values = c("dodgerblue","dodgerblue","dodgerblue","dodgerblue","dodgerblue","dodgerblue")) + # Customize colors
      ylim(-5,5)
    ggsave(paste0(output_file_name,'_',col,'.png'), path = "../results/latest1000all/violinPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)
  }
}

# Function to compute a Spearman correlation coeficient for the specified features
compute_correlation <- function(delta_zscores, selected_features) {
  corr_matrix_list <- list()
  for(feature in selected_features) {
    corr_matrix_list[[feature]] <- cor.test(delta_zscores$distance,
                                            delta_zscores[[feature]], 
                                            method = "spearman"
    )
  }
  return(corr_matrix_list)
}

# Custom K-S test function to keep sign
custom_ks_test <- function(x, y) {
  # Sort data
  x <- sort(x)
  y <- sort(y)
  
  # Calculate empirical cumulative distribution functions (ECDF)
  ecdf_x <- ecdf(x)
  ecdf_y <- ecdf(y)
  
  # Get the unique values from both samples
  unique_vals <- sort(unique(c(x, y)))
  
  # Calculate the raw differences between ECDFs
  diffs <- ecdf_x(unique_vals) - ecdf_y(unique_vals)
  
  # Find the maximum difference (positive or negative)
  max_diff <- max(diffs)
  min_diff <- min(diffs)
  
  signed_D <- 0
  if(max_diff > abs(min_diff)) {
    signed_D <- max_diff
  } else {
    signed_D <- min_diff
  }
  
  # Return the maximum and minimum differences
  return(list(signed_D = signed_D, max_diff = max_diff, min_diff = min_diff))
}

# Function to perform K-S test
run_ks_tests <- function(dataN, dataP, selected_features) {
  n <- length(selected_features)
  results <- matrix(nrow = 4, ncol = n)
  colnames(results) <- selected_features
  rownames(results) <- c("signed_D","max","min","p.val")
  
  for (i in selected_features) {
    #positive_col <- remove_outliers_IQR(dataP, i)
    #negative_col <- remove_outliers_IQR(dataN, i)
    positive_col <- dataP[[i]]
    negative_col <- dataN[[i]]
    
    ks_test <- custom_ks_test(negative_col, positive_col)
    ks_test_p <- ks.test(negative_col, positive_col)
    
    results[1,i] <- ks_test$signed_D
    results[2,i] <- ks_test$max_diff
    results[3,i] <- ks_test$min_diff
    results[4,i] <- ks_test_p$p.value
  }
  return(results)
}

################
################
# Main analisys:
################
# Load Data:
#source("load_gene_functionality_features.R")
funcProtExon2Data <- data.frame(Dataset = "protein-coding-exon2", read.csv("../data/latest1000all/functional-protein-exon2-dataset-features-w-fic.csv", header=TRUE))
funcProtExon3Data <- data.frame(Dataset = "protein-coding-exon3", read.csv("../data/latest1000all/functional-protein-exon3-dataset-features-w-fic.csv", header=TRUE))
funcLncrnaExon1Data <- data.frame(Dataset = "lncrna-exon1", read.csv("../data/latest1000all/functional-lncrna-exon1-dataset-features-w-fic.csv", header=TRUE))
funcLncrnaExon2Data <- data.frame(Dataset = "lncrna-exon2", read.csv("../data/latest1000all/functional-lncrna-exon2-dataset-features-w-fic.csv", header=TRUE))
funcSncrnaDataset <- data.frame(Dataset = "short-ncrna", read.csv("../data/latest1000all/functional-short-ncrna-dataset-features-w-fic.csv", header=TRUE))

protExon2NCData <- data.frame(Dataset = "protein-exon2-negative-control",read.csv("../data/latest1000all/protein-exon2-negative-control-dataset-features-w-fic.csv", header=TRUE))
protExon3NCData <- data.frame(Dataset = "protein-exon3-negative-control",read.csv("../data/latest1000all/protein-exon3-negative-control-dataset-features-w-fic.csv", header=TRUE))
lncrnaExon1NCData <- data.frame(Dataset = "lncrna-exon1-negative-control",read.csv("../data/latest1000all/lncrna-exon1-negative-control-dataset-features-w-fic.csv", header=TRUE))
lncrnaExon2NCData <- data.frame(Dataset = "lncrna-exon2-negative-control",read.csv("../data/latest1000all/lncrna-exon2-negative-control-dataset-features-w-fic.csv", header=TRUE))
sncrnaNCData <- data.frame(Dataset = "short-ncrna-negative-control",read.csv("../data/latest1000all/short-ncrna-negative-control-dataset-features-w-fic.csv", header=TRUE))
#str(lncrnaExon1NCData)

# Define features to analyze
loess_select_features <- c("GC_percentage",
                           "AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT",
                           "lowComplexity_density","phyloP_max","phyloP_max_100w","GERP_91_mammals_max","GERP_63_amniotes_max",
                           "RPKM_tissue","RPKM_primary.cell","copy_number","repeat_distance","Fickett_score","coding_potential","Max_covariance",
                           "MFE","accessibility","RNAalifold_score","Interaction_ave","SNP_density","MAF_avg",
                           "H3K27ac_AvgSignal","H3K36me3_AvgSignal","H3K79me2_AvgSignal","chrm_acc_AvgSignal","methylome"
                           )

loess_select_features_labels <- c("GC_content",
                                  "AA","AC","AG","AT","CA","CC","CpG","CT","GA","GC","GG","GT","TA","TC","TG","TT",
                                  "Low_complexity_density","phyloP_241_mammals","phyloP_100_vertebrates","GERP_91_eutherian_mammals","GERP_63_amniota_vertebrates",
                                  "RPKM_tissue","RPKM_primary_cell","Copy_number","Repeat_free","Fickett","RNAcode","Max_Covariance",
                                  "MFE","Accessibility","RNAalifold","Interaction_average","SNP_density","MAF",
                                  "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome")

cols1 <- colnames(funcProtExon2Data)
select_features_missing <- setdiff(cols1, loess_select_features_labels)
select_features_missing


# Rename features to new names used for paper
funcProtExon2Data <- rename_features(funcProtExon2Data, loess_select_features, loess_select_features_labels)
funcProtExon3Data <- rename_features(funcProtExon3Data, loess_select_features, loess_select_features_labels)
funcLncrnaExon1Data <- rename_features(funcLncrnaExon1Data, loess_select_features, loess_select_features_labels)
funcLncrnaExon2Data <- rename_features(funcLncrnaExon2Data, loess_select_features, loess_select_features_labels)
funcSncrnaDataset <- rename_features(funcSncrnaDataset, loess_select_features, loess_select_features_labels)

protExon2NCData <- rename_features(protExon2NCData, loess_select_features, loess_select_features_labels)
protExon3NCData <- rename_features(protExon3NCData, loess_select_features, loess_select_features_labels)
lncrnaExon1NCData <- rename_features(lncrnaExon1NCData, loess_select_features, loess_select_features_labels)
lncrnaExon2NCData <- rename_features(lncrnaExon2NCData, loess_select_features, loess_select_features_labels)
sncrnaNCData <- rename_features(sncrnaNCData, loess_select_features, loess_select_features_labels)

cols1 <- colnames(funcProtExon2Data)
select_features_missing <- setdiff(cols1, loess_select_features_labels)
select_features_missing

# add gene coords to negative data
#protExon2NCData <- get_gene_coords(protExon2NCData)
#protExon3NCData <- get_gene_coords(protExon3NCData)
#lncrnaExon1NCData <- get_gene_coords(lncrnaExon1NCData)
#lncrnaExon2NCData <- get_gene_coords(lncrnaExon2NCData)
#sncrnaNCData <- get_gene_coords(sncrnaNCData)

# add sequence length feature to data
funcProtExon2Data <- add_seq_length(funcProtExon2Data)
funcProtExon3Data <- add_seq_length(funcProtExon3Data)
funcLncrnaExon1Data <- add_seq_length(funcLncrnaExon1Data)
funcLncrnaExon2Data <- add_seq_length(funcLncrnaExon2Data)
funcSncrnaDataset <- add_seq_length(funcSncrnaDataset)

protExon2NCData <- add_seq_length(protExon2NCData, FALSE)
protExon3NCData <- add_seq_length(protExon3NCData, FALSE)
lncrnaExon1NCData <- add_seq_length(lncrnaExon1NCData, FALSE)
lncrnaExon2NCData <- add_seq_length(lncrnaExon2NCData, FALSE)
sncrnaNCData <- add_seq_length(sncrnaNCData, FALSE)


# Join positives and corresponding negatives in a single super set
proteinExon2_superset <- get_sequence_length_sets(funcProtExon2Data, protExon2NCData)
proteinExon3_superset <- get_sequence_length_sets(funcProtExon3Data, protExon3NCData)
lncrnaExon1_superset <- get_sequence_length_sets(funcLncrnaExon1Data, lncrnaExon1NCData)
lncrnaExon2_superset <- get_sequence_length_sets(funcLncrnaExon2Data, lncrnaExon2NCData)
sncrna_superset <- get_sequence_length_sets(funcSncrnaDataset, sncrnaNCData)

# Remove NAs
proteinExon2_superset_no_nas <- na.omit(proteinExon2_superset)
proteinExon3_superset_no_nas <- na.omit(proteinExon3_superset)
lncrnaExon1_superset_no_nas <- na.omit(lncrnaExon1_superset)
lncrnaExon2_superset_no_nas <- na.omit(lncrnaExon2_superset)
sncrna_superset_no_nas <- na.omit(sncrna_superset)

# Join by gene type
prot_superset <- rbind(proteinExon2_superset_no_nas,proteinExon3_superset_no_nas)
lncrna_superset <- rbind(lncrnaExon1_superset_no_nas,lncrnaExon2_superset_no_nas)

# Subset data into corresponding bins
prot_subsets <- subset_by_Distance_to_Gene(prot_superset)
lncrna_subsets <- subset_by_Distance_to_Gene(lncrna_superset)
sncrna_subsets <- subset_by_Distance_to_Gene(sncrna_superset_no_nas)

# Compute z-scores
prot_zscores_bins <- compute_zscores_all_bins(prot_subsets, loess_select_features_labels)
lncrna_zscores_bins <- compute_zscores_all_bins(lncrna_subsets, loess_select_features_labels)
sncrna_zscores_bins <- compute_zscores_all_bins(sncrna_subsets, loess_select_features_labels)

# Compute delta z-sccores for:
# Protein coding
prot_delta_zscores_bin1 <- as.data.frame(calculate_delta_zscores(as.data.frame(prot_zscores_bins$bin1$positive), 
                                                            as.data.frame(prot_zscores_bins$bin1$negative)))
prot_delta_zscores_bin2 <- as.data.frame(calculate_delta_zscores(as.data.frame(prot_zscores_bins$bin2$positive), 
                                              as.data.frame(prot_zscores_bins$bin2$negative)))
prot_delta_zscores_bin3 <- as.data.frame(calculate_delta_zscores(as.data.frame(prot_zscores_bins$bin3$positive), 
                                              as.data.frame(prot_zscores_bins$bin3$negative)))
prot_delta_zscores_bin4 <- as.data.frame(calculate_delta_zscores(as.data.frame(prot_zscores_bins$bin4$positive), 
                                              as.data.frame(prot_zscores_bins$bin4$negative)))
# lncRNA
lncrna_delta_zscores_bin1 <- as.data.frame(calculate_delta_zscores(as.data.frame(lncrna_zscores_bins$bin1$positive), 
                                              as.data.frame(lncrna_zscores_bins$bin1$negative)))
lncrna_delta_zscores_bin2 <- as.data.frame(calculate_delta_zscores(as.data.frame(lncrna_zscores_bins$bin2$positive), 
                                              as.data.frame(lncrna_zscores_bins$bin2$negative)))
lncrna_delta_zscores_bin3 <- as.data.frame(calculate_delta_zscores(as.data.frame(lncrna_zscores_bins$bin3$positive), 
                                              as.data.frame(lncrna_zscores_bins$bin3$negative)))
lncrna_delta_zscores_bin4 <- as.data.frame(calculate_delta_zscores(as.data.frame(lncrna_zscores_bins$bin4$positive), 
                                              as.data.frame(lncrna_zscores_bins$bin4$negative)))
# sncRNA
sncrna_delta_zscores_bin1 <- as.data.frame(calculate_delta_zscores(as.data.frame(sncrna_zscores_bins$bin1$positive), 
                                                                   as.data.frame(sncrna_zscores_bins$bin1$negative)))
sncrna_delta_zscores_bin2 <- as.data.frame(calculate_delta_zscores(as.data.frame(sncrna_zscores_bins$bin2$positive), 
                                                                   as.data.frame(sncrna_zscores_bins$bin2$negative)))
sncrna_delta_zscores_bin3 <- as.data.frame(calculate_delta_zscores(as.data.frame(sncrna_zscores_bins$bin3$positive), 
                                                                   as.data.frame(sncrna_zscores_bins$bin3$negative)))
sncrna_delta_zscores_bin4 <- as.data.frame(calculate_delta_zscores(as.data.frame(sncrna_zscores_bins$bin4$positive), 
                                                                   as.data.frame(sncrna_zscores_bins$bin4$negative)))

# Generate a dataframe suitable for ggplot.
prot_delta_zscores_bin1$bin <- "(1-5k)"
prot_delta_zscores_bin2$bin <- "(5k-50k)"
prot_delta_zscores_bin3$bin <- "(50k-500k)"
prot_delta_zscores_bin4$bin <- "(500k-5M)"

lncrna_delta_zscores_bin1$bin <- "(1-5k)"
lncrna_delta_zscores_bin2$bin <- "(5k-50k)"
lncrna_delta_zscores_bin3$bin <- "(50k-500k)"
lncrna_delta_zscores_bin4$bin <- "(500k-5M)"

sncrna_delta_zscores_bin1$bin <- "(1-5k)"
sncrna_delta_zscores_bin2$bin <- "(5k-50k)"
sncrna_delta_zscores_bin3$bin <- "(50k-500k)"
sncrna_delta_zscores_bin4$bin <- "(500k-5M)"



prot_delta_zscores_bin1$distance <- prot_subsets$`1-5k`$DistanceGene
prot_delta_zscores_bin2$distance <- prot_subsets$`5k-50k`$DistanceGene
prot_delta_zscores_bin3$distance <- prot_subsets$`50k-500k`$DistanceGene
prot_delta_zscores_bin4$distance <- prot_subsets$`500k-5M`$DistanceGene

lncrna_delta_zscores_bin1$distance <- lncrna_subsets$`1-5k`$DistanceGene
lncrna_delta_zscores_bin2$distance <- lncrna_subsets$`5k-50k`$DistanceGene
lncrna_delta_zscores_bin3$distance <- lncrna_subsets$`50k-500k`$DistanceGene
lncrna_delta_zscores_bin4$distance <- lncrna_subsets$`500k-5M`$DistanceGene

sncrna_delta_zscores_bin1$distance <- sncrna_subsets$`1-5k`$DistanceGene
sncrna_delta_zscores_bin2$distance <- sncrna_subsets$`5k-50k`$DistanceGene
sncrna_delta_zscores_bin3$distance <- sncrna_subsets$`50k-500k`$DistanceGene
sncrna_delta_zscores_bin4$distance <- sncrna_subsets$`500k-5M`$DistanceGene



prot_delta_zscores <- rbind(prot_delta_zscores_bin1,prot_delta_zscores_bin2,prot_delta_zscores_bin3,prot_delta_zscores_bin4)
prot_delta_zscores$bin <- factor(prot_delta_zscores$bin,levels = unique(prot_delta_zscores$bin))

lncrna_delta_zscores <- rbind(lncrna_delta_zscores_bin1,lncrna_delta_zscores_bin2,lncrna_delta_zscores_bin3,lncrna_delta_zscores_bin4)
lncrna_delta_zscores$bin <- factor(lncrna_delta_zscores$bin,levels = unique(lncrna_delta_zscores$bin))

sncrna_delta_zscores <- rbind(sncrna_delta_zscores_bin1,sncrna_delta_zscores_bin2,sncrna_delta_zscores_bin3,sncrna_delta_zscores_bin4)
sncrna_delta_zscores$bin <- factor(sncrna_delta_zscores$bin,levels = unique(sncrna_delta_zscores$bin))

# Generate violin plots for selected features
generate_violin_plot(prot_delta_zscores, loess_select_features_labels, "Protein Coding", "distance_effect_protein" )
generate_violin_plot(lncrna_delta_zscores, loess_select_features_labels, "lncRNA", "distance_effect_lncrna" )
generate_violin_plot(sncrna_delta_zscores, loess_select_features_labels, "sncRNA", "distance_effect_sncrna" )


# Compute correlation coefficient to reveal distance effect:
#Protein coding
prot_corr_matrix <- compute_correlation(prot_delta_zscores,
                                        loess_select_features_labels)

#sncRNA
sncrna_corr_matrix <- compute_correlation(sncrna_delta_zscores,
                                          loess_select_features_labels)

#lncRNA
lncrna_corr_matrix <- compute_correlation(lncrna_delta_zscores,
                                          loess_select_features_labels)


# Compute KS-test between first and last bin
# Protein coding
ks_results_prot_bin1_4 <- run_ks_tests(prot_delta_zscores_bin1, 
                                       prot_delta_zscores_bin4, 
                                       loess_select_features_labels)
# sncRNA
ks_results_sncrna_bin1_4 <- run_ks_tests(sncrna_delta_zscores_bin1, 
                                         sncrna_delta_zscores_bin4, 
                                         loess_select_features_labels)
# lncRNA
ks_results_lncrna_bin1_4 <- run_ks_tests(lncrna_delta_zscores_bin1, 
                                         lncrna_delta_zscores_bin4, 
                                         loess_select_features_labels)

################

# Compute ANOVA tests
## Set orthogonal contrasts.
op <- options(contrasts = c("contr.helmert", "contr.poly"))

one.way.prot <- aov(SNP_density ~ bin, data = prot_delta_zscores)
one.way.lincrna <- aov(SNP_density ~ bin, data = lncrna_delta_zscores)
one.way.sincrna <- aov(SNP_density ~ bin, data = sncrna_delta_zscores)

summary(one.way.prot)
summary(one.way.lincrna)
summary(one.way.sincrna)

# Perform test of "homoscedasticity"
par(mfrow=c(2,2))
plot(one.way.prot, main="SNP density - Protein coding")
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(one.way.lincrna, main="GC content - lincRNA")
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(one.way.sincrna, main="GC content - sincRNA")
par(mfrow=c(1,1))

# Post-hoc test
tukey.one.way.prot<-TukeyHSD(one.way.prot)
tukey.one.way.prot

par(mar = c(5.1, 9.2, 4.1, 2.1))
plot(tukey.one.way.prot, las = 1)
title(main = "SNP_density - Protein Coding", line = 1)

tukey.one.way.lincrna<-TukeyHSD(one.way.lincrna)
tukey.one.way.lincrna
plot(tukey.one.way.lincrna, las = 1)
title(main = "GC content - lincRNA", line = 1)

tukey.one.way.sincrna<-TukeyHSD(one.way.sincrna)
tukey.one.way.sincrna
plot(tukey.one.way.sincrna, las = 1)
title(main = "GC content - sincRNA", line = 1)

options(op)  # reset to previous contrasts
