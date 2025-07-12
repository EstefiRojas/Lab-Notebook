# Effect on Distance Analysis, delta z-scores based on distance to gene
# Try to identify interesting effects on distance for features of functional genes
# Load libraries:
library(dplyr)
library(tidyr)
library(ggplot2)
install.packages("corrr")
library('corrr')
install.packages("Hmisc")
library(Hmisc)
library(patchwork)

# Load Data:
source("load_gene_functionality_features.R")
colnames(feature_matrix)
# Define features to analyze
eod_select_features <- c("GC_percentage", "lowComplexity_density", 
                         "CpG", "GA", "GG", "GT", "TA", "AC", "CC",
                         "phyloP_max_241w", "phyloP_max_100w", 
                         "GERP_91_mammals_max", "GERP_63_amniotes_max", 
                         "RPKM_tissue", "RPKM_primary.cell", 
                         "H3K9ac_MaxScaledSignal", "H3K79me2_MaxScaledSignal", 
                         "chrm_acc_MaxScaledSignal", "methylome", 
                         "repeat_distance", "copy_number", "coding_potential", 
                         "fickett", "Max_covariance", "MFE",
                         "Interaction_ave", "SNP_density", "MAF_avg"
                         )

eod_select_features_labels <- c("GC%", "Complexity",
                                "CpG",  "GA", "GG", "GT", "TA", "AC", "CC",
                                "PhyloP-mammals", "PhyloP-vertebrates",
                                "GERP-mammals", "GERP-vertebrates",
                                "Tissue RPKM", "Primary cell RPKM", 
                                "H3K9ac", "H3K79me2",
                                "Chromatin", "Methylome",
                                "Repeat.free", "Copies", "RNAcode",
                                "Fickett", "Covariance", "MFE",
                                "Interactions", "SNPs", "MAF"
                                )


# Function to rename features to new paper names
rename_features <- function(dataframe, current_names, new_names) {
  for (i in 1:length(current_names)) {  # Iterate from 1 to length of list of names
    names(dataframe)[ names(dataframe) == current_names[i] ] <- new_names[i]
  }
  return(dataframe)
}
# Rename features to new names used for paper
names_differences <- setdiff(eod_select_features_labels, colnames(funcProtExon2Data))
names_differences
funcProtExon2Data <- rename_features(funcProtExon2Data, eod_select_features, eod_select_features_labels)
funcProtExon3Data <- rename_features(funcProtExon3Data, eod_select_features, eod_select_features_labels)
funcLncrnaExon1Data <- rename_features(funcLncrnaExon1Data, eod_select_features, eod_select_features_labels)
funcLncrnaExon2Data <- rename_features(funcLncrnaExon2Data, eod_select_features, eod_select_features_labels)
funcSncrnaDataset <- rename_features(funcSncrnaDataset, eod_select_features, eod_select_features_labels)

protExon2NCData <- rename_features(protExon2NCData, eod_select_features, eod_select_features_labels)
protExon3NCData <- rename_features(protExon3NCData, eod_select_features, eod_select_features_labels)
lncrnaExon1NCData <- rename_features(lncrnaExon1NCData, eod_select_features, eod_select_features_labels)
lncrnaExon2NCData <- rename_features(lncrnaExon2NCData, eod_select_features, eod_select_features_labels)
sncrnaNCData <- rename_features(sncrnaNCData, eod_select_features, eod_select_features_labels)

### This is not done anymore ###
# Function to add sequence length as a feature
#add_seq_length <- function(dataset, is_positive=TRUE) {
#  # Compute length of sequence by subtracting Start coord to End coord only once
#  seq_length <- dataset$End - dataset$Start
#  
#  if(!is_positive){
#    dataset_1 <- dataset
#    dataset_2 <- dataset
#    dataset_exact <- dataset
#    set_length <- dataset$End - dataset$Start
#    dataset_1$seqLength <- seq_length - 1
#    dataset_2$seqLength <- seq_length - 2
#    dataset_exact$seqLength <- seq_length
#    dataset <- rbind(dataset_exact, dataset_1, dataset_2)
#  } else
#    dataset$seqLength <- seq_length
#  return(dataset)
#}

# Select just desired features
funcProtExon2DataSelect <- funcProtExon2Data %>% select(all_of(eod_select_features))
funcProtExon3DataSelect <- funcProtExon3Data %>% select(all_of(eod_select_features))
funcLncrnaExon1DataSelect <- funcLncrnaExon1Data %>% select(all_of(eod_select_features))
funcLncrnaExon2DataSelect <- funcLncrnaExon2Data %>% select(all_of(eod_select_features))
funcSncrnaDatasetSelect <- funcSncrnaDataset %>% select(all_of(eod_select_features))

protExon2NCDataSelect <- protExon2NCData %>% select(all_of(eod_select_features))
protExon3NCDataSelect <- protExon3NCData %>% select(all_of(eod_select_features))
lncrnaExon1NCDataSelect <- lncrnaExon1NCData %>% select(all_of(eod_select_features))
lncrnaExon2NCDataSelect <- lncrnaExon2NCData %>% select(all_of(eod_select_features))
sncrnaNCDataSelect <- sncrnaNCData %>% select(all_of(eod_select_features))

# add gene id feature to data
funcProtExon2DataSelect$GeneID <- as.data.frame(read.csv("../results/functional-protein-exon2-dataset-features.csv", header = TRUE))$GeneID
funcProtExon3DataSelect$GeneID <- as.data.frame(read.csv("../results/functional-protein-exon3-dataset-features.csv", header = TRUE))$GeneID
funcLncrnaExon1DataSelect$GeneID <- as.data.frame(read.csv("../results/functional-lncrna-exon1-dataset-features.csv", header = TRUE))$GeneID
funcLncrnaExon2DataSelect$GeneID <- as.data.frame(read.csv("../results/functional-lncrna-exon2-dataset-features.csv", header = TRUE))$GeneID
funcSncrnaDatasetSelect$GeneID <- as.data.frame(read.csv("../results/functional-short-ncrna-dataset-features.csv", header = TRUE))$GeneID

#prot_exon2_negative_control_data <- as.data.frame(read.csv("../data/datasets/protein-exon2-negative-control.csv", header = TRUE))

#protExon2NCData2Select <- protExon2NCDataSelect[1:9620,]
#protExon3NCData2Select <- protExon3NCDataSelect[1:9615,]
#ncrnaExon1NCData2Select <- lncrnaExon1NCDataSelect[1:9801,]
#ncrnaExon2NCData2Select <- lncrnaExon2NCDataSelect
#sncrnaNCData2Select <- sncrnaNCDataSelect[1:9833,]

protExon2NCDataSelect$GeneID <- as.data.frame(read.csv("../data/datasets/protein-exon2-negative-control.csv", header = TRUE))$ID.1
protExon3NCDataSelect$GeneID <- as.data.frame(read.csv("../data/datasets/protein-exon3-negative-control.csv", header = TRUE))$ID.1
lncrnaExon1NCDataSelect$GeneID <- as.data.frame(read.csv("../data/datasets/lncrna-exon1-negative-control.csv", header = TRUE))$ID.1
lncrnaExon2NCDataSelect$GeneID <- as.data.frame(read.csv("../data/datasets/lncrna-exon2-negative-control.csv", header = TRUE))$ID.1[1:9798]
sncrnaNCDataSelect$GeneID <- as.data.frame(read.csv("../data/datasets/short-ncrna-negative-control-dataset.csv", header = TRUE))$ID.1

# add Chromosome feature to data
#funcProtExon2DataSelect$Chromosome <- as.data.frame(read.csv("../results/functional-protein-exon2-dataset-features.csv", header = TRUE))$Chromosome
#funcProtExon3DataSelect$Chromosome <- as.data.frame(read.csv("../results/functional-protein-exon3-dataset-features.csv", header = TRUE))$Chromosome
#uncLncrnaExon1DataSelect$Chromosome <- as.data.frame(read.csv("../results/functional-lncrna-exon1-dataset-features.csv", header = TRUE))$Chromosome
#funcLncrnaExon2DataSelect$Chromosome <- as.data.frame(read.csv("../results/functional-lncrna-exon2-dataset-features.csv", header = TRUE))$Chromosome
#uncSncrnaDatasetSelect$Chromosome <- as.data.frame(read.csv("../results/functional-short-ncrna-dataset-features.csv", header = TRUE))$Chromosome

#protExon2NCData2Select$Chromosome <- as.data.frame(read.csv("../data/datasets/protein-exon2-negative-control.csv", header = TRUE))$Chromosome
#protExon3NCData2Select$Chromosome <- as.data.frame(read.csv("../data/datasets/protein-exon3-negative-control.csv", header = TRUE))$Chromosome
#lncrnaExon1NCData2Select$Chromosome <- as.data.frame(read.csv("../data/datasets/lncrna-exon1-negative-control.csv", header = TRUE))$Chromosome
#lncrnaExon2NCData2Select$Chromosome <- as.data.frame(read.csv("../data/datasets/lncrna-exon2-negative-control.csv", header = TRUE))$Chromosome[1:9798]
#ncrnaNCData2Select$Chromosome <- as.data.frame(read.csv("../data/datasets/short-ncrna-negative-control-dataset.csv", header = TRUE))$Chromosome
# OR
funcProtExon2DataSelect$Chromosome <- funcProtExon2Data$Chromosome
funcProtExon3DataSelect$Chromosome <- funcProtExon3Data$Chromosome
funcLncrnaExon1DataSelect$Chromosome <- funcLncrnaExon1Data$Chromosome
funcLncrnaExon2DataSelect$Chromosome <- funcLncrnaExon2Data$Chromosome
funcSncrnaDatasetSelect$Chromosome <- funcSncrnaDataset$Chromosome

protExon2NCDataSelect$Chromosome <- protExon2NCData$Chromosome
protExon3NCDataSelect$Chromosome <- protExon3NCData$Chromosome
lncrnaExon1NCDataSelect$Chromosome <- lncrnaExon1NCData$Chromosome
lncrnaExon2NCDataSelect$Chromosome <- lncrnaExon2NCData$Chromosome
sncrnaNCDataSelect$Chromosome <- sncrnaNCData$Chromosome

#####
# add DistanceGene feature to zscores data
#protExon2NCData$DistanceGene <- as.data.frame(read.csv("../data/latest1000all/protein-exon2-negative-control-dataset-features.csv", header = TRUE))$DistanceGene
#protExon3NCData$DistanceGene <- as.data.frame(read.csv("../data/latest1000all/protein-exon3-negative-control-dataset-features.csv", header = TRUE))$DistanceGene
#lncrnaExon1NCData$DistanceGene <- as.data.frame(read.csv("../data/latest1000all/lncrna-exon1-negative-control-dataset-features.csv", header = TRUE))$DistanceGene
#lncrnaExon2NCData$DistanceGene <- as.data.frame(read.csv("../data/latest1000all/lncrna-exon2-negative-control-dataset-features.csv", header = TRUE))$DistanceGene[1:9798]
#sncrnaNCData$DistanceGene <- as.data.frame(read.csv("../data/latest1000all/short-ncrna-negative-control-dataset-features.csv", header = TRUE))$DistanceGene
# OR
protExon2NCDataSelect$DistanceGene <- protExon2NCData$DistanceGene
protExon3NCDataSelect$DistanceGene <- protExon3NCData$DistanceGene
lncrnaExon1NCDataSelect$DistanceGene <- lncrnaExon1NCData$DistanceGene
lncrnaExon2NCDataSelect$DistanceGene <- lncrnaExon2NCData$DistanceGene
sncrnaNCDataSelect$DistanceGene <- sncrnaNCData$DistanceGene
# OR 
list_with_zscores_protein_neg$DistanceGene <- as.data.frame(rbind(read.csv("../data/latest1000all/protein-exon2-negative-control-dataset-features.csv", header = TRUE),
                                                                  read.csv("../data/latest1000all/protein-exon3-negative-control-dataset-features.csv", header = TRUE)))$DistanceGene
list_with_zscores_lncrna_neg$DistanceGene <- as.data.frame(rbind(read.csv("../data/latest1000all/lncrna-exon1-negative-control-dataset-features.csv", header = TRUE),
                                                                 read.csv("../data/latest1000all/lncrna-exon2-negative-control-dataset-features.csv", header = TRUE)))$DistanceGene
list_with_zscores_sncrna_neg$DistanceGene <- as.data.frame(read.csv("../data/latest1000all/short-ncrna-negative-control-dataset-features.csv", header = TRUE))$DistanceGene

# OR 
testdata <- rbind(protExon2NCData, protExon3NCData)
length(testdata)
list_with_zscores_protein_neg$DistanceGene <- testdata$DistanceGene
list_with_zscores_lncrna_neg$DistanceGene <- rbind(lncrnaExon1NCData, lncrnaExon2NCData)$DistanceGene
list_with_zscores_sncrna_neg$DistanceGene <- rbind(sncrnaNCData)$DistanceGene


list_with_zscores_protein$DistanceGene <- 1
list_with_zscores_lncrna$DistanceGene <- 1
list_with_zscores_sncrna$DistanceGene <- 1

#####

write.csv(protExon2NCData2Select, "../protExon2NCData2Select.csv")
write.csv(protExon2NCData2Select, "../protExon3NCData2Select")
write.csv(lncrnaExon1NCData2Select, "../lncrnaExon1NCData2Select.csv")
write.csv(lncrnaExon2NCData2Select, "../lncrnaExon2NCData2Select")
write.csv(sncrnaNCData2Select, "../sncrnaNCData2Select")



# Function to join positive data to corresponding negative data by chromosome and gene id
join_by_chromosome_and_geneid <- function(positive_df, negative_df) {
  joinedSet <- left_join(negative_df,
                         positive_df,
                         by = join_by(Chromosome, GeneID),
                         suffix = c(".neg",".pos"),
                         relationship = "many-to-one",
                         multiple = "any",
                         unmatched = "drop",
                         )
  return(joinedSet)
}

# Join negatives and corresponding positives in a single super set
proteinExon2_superset <- join_by_chromosome_and_geneid(funcProtExon2DataSelect, protExon2NCDataSelect)
proteinExon3_superset <- join_by_chromosome_and_geneid(funcProtExon3DataSelect, protExon3NCDataSelect)
lncrnaExon1_superset <- join_by_chromosome_and_geneid(funcLncrnaExon1DataSelect, lncrnaExon1NCDataSelect)
lncrnaExon2_superset <- join_by_chromosome_and_geneid(funcLncrnaExon2DataSelect, lncrnaExon2NCDataSelect)
sncrna_superset <- join_by_chromosome_and_geneid(funcSncrnaDatasetSelect, sncrnaNCData2Select)


# Define a function to detect numeric columns based on the ability to convert values.
is_numeric_column <- function(x) {
  # Try converting to numeric and check how many NAs are produced
  num_na <- sum(is.na(as.numeric(as.character(x))))
  total <- length(x)
  # If less than 50% of the values become NA, we can consider it a numeric column
  return(num_na / total < 0.5)
}

# Apply the function to each column to get a logical vector of numeric columns
convert_non_numeric_to_na <- function(dataframe) {
  numeric_columns <- sapply(dataframe, is_numeric_column)
  
  # Step 2: For each numerical column, convert non-numeric values to NA
  for (col_name in names(dataframe)[numeric_columns]) {
    dataframe[[col_name]] <- as.numeric(as.character(dataframe[[col_name]]))
  }
  return(dataframe)
}

#feature_matrix_numeric <- as.data.frame(sapply(feature_matrix_numeric, convert_to_numeric))
proteinExon2_superset_no_nas <- na.omit(convert_non_numeric_to_na(proteinExon2_superset))
proteinExon3_superset_no_nas <- na.omit(convert_non_numeric_to_na(proteinExon3_superset))
lncrnaExon1_superset_no_nas <- na.omit(convert_non_numeric_to_na(lncrnaExon1_superset))
lncrnaExon2_superset_no_nas <- na.omit(convert_non_numeric_to_na(lncrnaExon2_superset))
sncrna_superset_no_nas <- na.omit(convert_non_numeric_to_na(sncrna_superset))


# Join by gene type
# WITH NAs:
prot_superset <- rbind(proteinExon2_superset,proteinExon3_superset)
lncrna_superset <- rbind(lncrnaExon1_superset,lncrnaExon2_superset)
#OR without NAs
prot_superset <- rbind(proteinExon2_superset_no_nas,proteinExon3_superset_no_nas)
lncrna_superset <- rbind(lncrnaExon1_superset_no_nas,lncrnaExon2_superset_no_nas)


# Function to subset negative datasets by Distance to Gene:
subset_by_Distance_to_Gene <- function(dataset) {
  list_of_subsets <- list()
  list_of_subsets[['1-5k']] <- dataset[dataset$DistanceGene > 1000 & dataset$DistanceGene <= 5000,]
  list_of_subsets[['5k-50k']] <- dataset[dataset$DistanceGene > 5000 & dataset$DistanceGene <= 50000,]
  list_of_subsets[['50k-500k']] <- dataset[dataset$DistanceGene > 50000 & dataset$DistanceGene <= 500000,]
  list_of_subsets[['500k-5M']] <- dataset[dataset$DistanceGene > 500000 & dataset$DistanceGene <= 5000000,]
  return(list_of_subsets)
}
# Subset data into corresponding bins
prot_subsets <- subset_by_Distance_to_Gene(prot_superset)
lncrna_subsets <- subset_by_Distance_to_Gene(lncrna_superset)
sncrna_subsets <- subset_by_Distance_to_Gene(sncrna_superset_no_nas)

help("mad")
# Function to compute robust z-scores
get_robust_zscores <- function(subsets_df, selected_features, is_negative=FALSE) {
  list_with_zscores <- list()
  list_with_method <- list()
  for (col in selected_features) {
    col_name <- if(is_negative) paste0(col,".neg") else paste0(col,".pos")
    positive_col <- as.numeric(subsets_df[[col_name]])
    negative_col <- as.numeric(subsets_df[[paste0(col,".neg")]])

    mad_value <- mad(negative_col, na.rm = TRUE)
    median_value <- median(negative_col, na.rm = TRUE)
    
    mean_value <- mean(negative_col, na.rm = TRUE) 
    sd_value <- sd(negative_col, na.rm = TRUE)
    absolute_deviations <- abs(negative_col - mean_value)
    meanAD_value <- mean(absolute_deviations, na.rm = TRUE)
    
    if(mad_value != 0) { # MAD method
      list_with_zscores[[col]] <- (positive_col - median_value) / mad_value
      list_with_method[[col]] <- "MAD"
    } else if(meanAD_value != 0) { # meanAD estimator
      list_with_zscores[[col]] <- (positive_col - median_value) / (1.2533 * meanAD_value)
      list_with_method[[col]] <- "meanAD"
    } else { # regular z-score
      list_with_zscores[[col]] <- (positive_col - mean_value) / sd_value
      list_with_method[[col]] <- "regular z-score"
    }
  }
  return(list(zscores=list_with_zscores, method=list_with_method))
}
# Function to compute z-scores for each bin:
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
# Compute z-scores
prot_zscores_bins <- compute_zscores_all_bins(prot_subsets, eod_select_features_labels)
lncrna_zscores_bins <- compute_zscores_all_bins(lncrna_subsets, eod_select_features_labels)
sncrna_zscores_bins <- compute_zscores_all_bins(sncrna_subsets, eod_select_features_labels)


# Function to compute deltas
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

# Function to plot delta z-scores against Distance to Gene bins
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
      scale_fill_manual(values = c("dodgerblue","dodgerblue","dodgerblue","dodgerblue")) + # Customize colors
      ylim(-5,5)
    ggsave(paste0(output_file_name,'_',col,'.png'), path = "../results/latest1000all/violinPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)
  }
}
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

# Add distance to delta z-scores for correlation coefficient
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

# Join in one dataframe for each gene type
prot_delta_zscores <- rbind(prot_delta_zscores_bin1,prot_delta_zscores_bin2,prot_delta_zscores_bin3,prot_delta_zscores_bin4)
prot_delta_zscores$bin <- factor(prot_delta_zscores$bin,levels = unique(prot_delta_zscores$bin))

lncrna_delta_zscores <- rbind(lncrna_delta_zscores_bin1,lncrna_delta_zscores_bin2,lncrna_delta_zscores_bin3,lncrna_delta_zscores_bin4)
lncrna_delta_zscores$bin <- factor(lncrna_delta_zscores$bin,levels = unique(lncrna_delta_zscores$bin))

sncrna_delta_zscores <- rbind(sncrna_delta_zscores_bin1,sncrna_delta_zscores_bin2,sncrna_delta_zscores_bin3,sncrna_delta_zscores_bin4)
sncrna_delta_zscores$bin <- factor(sncrna_delta_zscores$bin,levels = unique(sncrna_delta_zscores$bin))

# Generate violin plots for selected features
generate_violin_plot(prot_delta_zscores, eod_select_features_labels, "Protein Coding", "distance_effect_protein" )
generate_violin_plot(lncrna_delta_zscores, eod_select_features_labels, "lncRNA", "distance_effect_lncrna" )
generate_violin_plot(sncrna_delta_zscores, eod_select_features_labels, "sncRNA", "distance_effect_sncrna" )


# Function to plot delta z-scores against Distance to Gene
generate_jitter_plot <- function(delta_zscores, selected_features, title, output_file_name) {
  for(col in selected_features) {
    ggplot(delta_zscores, aes(x = distance, y = !!sym(col), color = bin)) +
      geom_jitter(width = 0.2, height = 0, show.legend = TRUE) +
      geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.01) +
      facet_wrap(~ "", scales = "free") +
      labs(title = paste(title, '-', gsub("_", " ", col)), x = "Distance to Gene", y = "Delta Z-score") +
      theme_minimal(base_size = 26) +
      theme(
        axis.text.x = element_text(size = 20, angle = 0),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 36, hjust = 0.5)
      ) +
      scale_x_log10() +
      scale_color_manual(values = c("dodgerblue", "dodgerblue", "dodgerblue", "dodgerblue"))
    ggsave(paste0(output_file_name,'_',col,'.png'), path = "../results/latest1000all/jitterPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)
  }
}

# Generate scatter plots for selected features
generate_jitter_plot(prot_delta_zscores, eod_select_features_labels, "Protein Coding", "distance_effect_protein_lowess" )
generate_jitter_plot(lncrna_delta_zscores, eod_select_features_labels, "lncRNA", "distance_effect_lncrna" )
generate_jitter_plot(sncrna_delta_zscores, eod_select_features_labels, "sncRNA", "distance_effect_sncrna" )


# Function to plot delta z-scores against Distance to Gene bins
generate_scatter_plot_log_scale <- function(delta_zscores, selected_features, title, output_file_name) {
  for(col in selected_features) {
    ggplot(delta_zscores, aes(x = distance, y = !!sym(col), color = bin)) +
      geom_point(show.legend = TRUE) +
      geom_smooth(method = "loess", se = TRUE, color = "red", size = 1) +
      facet_wrap(~ "", scales = "free") +
      labs(title = paste(title, '-', gsub("_", " ", col)), x = "log(Distance to Gene)", y = "Delta Z-score") +
      theme_minimal(base_size = 26) +
      theme(
        axis.text.x = element_text(size = 20, angle = 0),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 36, hjust = 0.5)
      ) +
      scale_x_log10() +
      scale_color_manual(values = c("dodgerblue", "dodgerblue", "dodgerblue", "dodgerblue"))
    ggsave(paste0(output_file_name,'_',col,'.png'), path = "../results/latest1000all/logScatterPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)
  }
}

# Generate scatter plots for selected features
generate_scatter_plot_log_scale(prot_delta_zscores, eod_select_features_labels, "Protein Coding", "distance_effect_protein_loess_smooth" )
generate_scatter_plot_log_scale(lncrna_delta_zscores, eod_select_features_labels, "lncRNA", "distance_effect_lncrna_loess_smooth" )
generate_scatter_plot_log_scale(sncrna_delta_zscores, eod_select_features_labels, "sncRNA", "distance_effect_sncrna_loess_smooth" )


prot_neg_zscores <- rbind(as.data.frame(prot_zscores_bins$bin1$negative),
                          as.data.frame(prot_zscores_bins$bin2$negative),
                          as.data.frame(prot_zscores_bins$bin3$negative),
                          as.data.frame(prot_zscores_bins$bin4$negative))
prot_neg_zscores$distance <- prot_delta_zscores$distance
prot_neg_zscores$type <- "neg"

prot_pos_zscores$type <- "pos"

prot_zscores <- rbind(prot_neg_zscores,prot_pos_zscores)
prot_zscores$type <- factor(prot_zscores$type, levels = c("neg","pos"))
prot_neg_zscores <- prot_zscores[prot_zscores$type == "neg",]
prot_pos_zscores <- prot_zscores[prot_zscores$type == "pos",]
prot_pos_median_zscores <- sample_n(prot_pos_zscores,20) %>% select(GC.content,SNP.density,gnomAD_MAF)
#prot_pos_median_zscores <- as.data.frame(t(apply(select(prot_pos_zscores, GC.content,SNP.density,gnomAD_MAF), 2, median, na.rm = TRUE)))
prot_pos_median_zscores$distance <- 1
prot_pos_median_zscores$type <- "pos"
prot_pos_median_zscores$X <- 1
prot_neg_zscores[prot_neg_zscores$distance==1,]$distance <- NA

prot_zscores <- rbind(prot_neg_zscores, prot_pos_median_zscores)

library(dplyr)
install.packages("LSD")
library(LSD)
 
write(prot_zscores, "prot_zscores.csv")
prot_zscores <- read.csv("prot_zscores.csv")
prot_plot <- ggplot(prot_neg_zscores, aes(x = distance, y = GC.content, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_bin2d() + 
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.1) +
  facet_wrap(~ "", scales = "free") +
  labs(title = "Protein Coding - GC content", x = "Distance to Gene", y = "Robust Z-score") +
  theme_minimal(base_size = 26) +
  coord_cartesian(xlim = c(7500, 5100000), ylim = c(-2,2)) +
  theme(
    axis.text.x = element_text(size = 20, angle = 0),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    legend.position = "none",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 36, hjust = 0.5)
  ) +
  scale_color_manual(values = c("dodgerblue", "orange")) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(1000, 10000, 100000,1000000,10000000),
                     labels = scales::label_number()
  )

prot_plot

prot_pos_zscores <- rbind(as.data.frame(prot_zscores_bins$bin1$positive),
                          as.data.frame(prot_zscores_bins$bin2$positive),
                          as.data.frame(prot_zscores_bins$bin3$positive),
                          as.data.frame(prot_zscores_bins$bin4$positive))
prot_pos_zscores$distance <- 1
prot_plot +
  geom_point(data = prot_pos_zscores, aes(x = log10(distance), y = GC.content), color = "orange", show.legend = TRUE)
ggsave(paste0("distance_effect_protein_loess_wiggle_veplus",'_',"GCcontent",'.png'), path = "../results/latest1000all/logScatterPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

prot_pos_zscores$group <- "mrna(+)"
prot_neg_zscores$group <- "mrna(-)"
prot_zscores_all <- rbind(prot_pos_zscores,prot_neg_zscores)
prot_zscores_all$group <- factor(prot_zscores_all$group, levels = c("mrna(+)","mrna(-)"))

prot_plot <- ggplot(prot_neg_zscores, aes(x = distance, y = SNP.density)) +
  geom_point(show.legend = TRUE, color = "dodgerblue") +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.01) +
  facet_wrap(~ "", scales = "free") +
  labs(title = "Protein Coding - SNP Density", x = "Distance to Gene", y = "Robust Z-score") +
  theme_minimal(base_size = 26) +
  coord_cartesian(ylim = c(-3,5)) +
  theme(
    axis.text.x = element_text(size = 20, angle = 0),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    legend.position = "none",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 36, hjust = 0.5)
  ) +
  scale_color_manual(values = c("dodgerblue", "dodgerblue", "dodgerblue", "dodgerblue"))
prot_plot
prot_plot +
  geom_point(data = prot_pos_zscores, aes(x = log10(distance), y = SNP.density), color = "orange", show.legend = TRUE)
ggsave(paste0("distance_effect_protein_loess_wiggle_veplus",'_',"SNPdensity",'.png'), path = "../results/latest1000all/logScatterPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)



# PROTEIN: #
unique(zscores_data$Dataset)
unique(protExon2NCData$Dataset)
prot_pos_zscores <- zscores_data %>%
  filter(Dataset == "protein-coding-exon2" | Dataset == "protein-coding-exon3") %>%
  mutate(distance=1)

prot_neg_zscores <- zscores_data %>%
  filter(Dataset == "protein-exon2-negative-control" | Dataset == "protein-exon3-negative-control") %>%
  mutate(distance = rbind(
    protExon2NCData %>% filter(Dataset == "protein-exon2-negative-control") %>% dplyr::select(DistanceGene),
    protExon3NCData %>% filter(Dataset == "protein-exon3-negative-control") %>% dplyr::select(DistanceGene)
  )$DistanceGene)


prot_pos_zscores$group <- "mrna(+)"
prot_neg_zscores$group <- "mrna(-)"
prot_zscores_all <- rbind(prot_pos_zscores,prot_neg_zscores)
prot_zscores_all$group <- factor(prot_zscores_all$group, levels = c("mrna(+)","mrna(-)"))

# GC%
prot_gc_plot <- ggplot(prot_neg_zscores, aes(x = distance, y = GC.content, color = group)) +
  geom_point(show.legend = TRUE) +
  #geom_bin2d() +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "mRNA") +
  xlab("Distance to Gene (Mb)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  scale_x_continuous(trans='log10', breaks = c(1000, 10000, 100000, 1000000, 5000000), 
                     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
                     limits = c(1000, 5000000)) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_color_manual(values = c("orange", "dodgerblue"))
ggsave(paste0("distance_effect_protein_loess_wiggle_veplus",'_',"GCcontent",'.png'), path = "../results/latest1000all/logScatterPlots/Paper/", prot_gc_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# SNP density
prot_snp_plot <- ggplot(prot_zscores_all, aes(x = log10(distance), y = SNP.density, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.5) +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "mRNA") +
  xlab("log(Distance to Gene)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-2,4)) +
  scale_color_manual(values = c("orange", "dodgerblue"))
ggsave(paste0("distance_effect_protein_loess_wiggle_veplus",'_',"SNPdensity",'.png'), path = "../results/latest1000all/logScatterPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)




# LNCRNA: #
unique(zscores_data$Dataset)
unique(lncrnaExon2NCData$Dataset)
lncrna_pos_zscores <- zscores_data %>%
  filter(Dataset == "lncrna-exon1" | Dataset == "lncrna-exon2") %>%
  mutate(distance=1)

lncrna_neg_zscores <- zscores_data %>%
  filter(Dataset == "lncrna-exon1-negative-control" | Dataset == "lncrna-exon2-negative-control") %>%
  mutate(distance = rbind(
    lncrnaExon1NCData %>% filter(Dataset == "lncrna-exon1-negative-control") %>% dplyr::select(DistanceGene),
    lncrnaExon2NCData %>% filter(Dataset == "lncrna-exon2-negative-control") %>% dplyr::select(DistanceGene)
  )$DistanceGene)


lncrna_pos_zscores$group <- "lncrna(+)"
lncrna_neg_zscores$group <- "lncrna(-)"
lncrna_zscores_all <- rbind(lncrna_pos_zscores,lncrna_neg_zscores)
lncrna_zscores_all$group <- factor(lncrna_zscores_all$group, levels = c("lncrna(+)","lncrna(-)"))

# GC%
lncrna_gc_plot <- ggplot(lncrna_zscores_all, aes(x = log10(distance), y = GC.content, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.5) +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "lncRNA") +
  xlab("log(Distance to Gene)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_color_manual(values = c("orange", "dodgerblue"))
ggsave(paste0("distance_effect_lncrna_loess_wiggle_veplus",'_',"GCcontent",'.png'), path = "../results/latest1000all/logScatterPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)



# SNP density
lncrna_snp_plot <- ggplot(lncrna_zscores_all, aes(x = log10(distance), y = SNP.density, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.5) +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "lncRNA") +
  xlab("log(Distance to Gene)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-2,4)) +
  scale_color_manual(values = c("orange", "dodgerblue"))
ggsave(paste0("distance_effect_lncrna_loess_wiggle_veplus",'_',"SNPdensity",'.png'), path = "../results/latest1000all/logScatterPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)



# SNCRNA: #
sncrna_pos_zscores <- zscores_data %>%
  filter(Dataset == "short-ncrna") %>%
  mutate(distance=1)

sncrna_neg_zscores <- zscores_data %>%
  filter(Dataset == "short-ncrna-negative-control") %>%
  mutate(distance = (sncrnaNCData %>% 
           filter(Dataset == "short-ncrna-negative-control") %>% 
           dplyr::select(DistanceGene))$DistanceGene
  )


sncrna_pos_zscores$group <- "sncrna(+)"
sncrna_neg_zscores$group <- "sncrna(-)"
sncrna_zscores_all <- rbind(sncrna_pos_zscores,sncrna_neg_zscores)
sncrna_zscores_all$group <- factor(sncrna_zscores_all$group, levels = c("sncrna(+)","sncrna(-)"))

# GC%
sncrna_gc_plot <- ggplot(sncrna_zscores_all, aes(x = log10(distance), y = GC.content, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.5) +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "sncRNA") +
  xlab("log(Distance to Gene)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_color_manual(values = c("orange", "dodgerblue"))
ggsave(paste0("distance_effect_sncrna_loess_wiggle_veplus",'_',"GCcontent",'.png'), path = "../results/latest1000all/logScatterPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


# SNP density
sncrna_snp_plot <- ggplot(sncrna_zscores_all, aes(x = log10(distance), y = SNP.density, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.5) +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "sncRNA") +
  xlab("log(Distance to Gene)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-2,4)) +
  scale_color_manual(values = c("orange", "dodgerblue"))
ggsave(paste0("distance_effect_sncrna_loess_wiggle_veplus",'_',"SNPdensity",'.png'), path = "../results/latest1000all/logScatterPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)



###########################
distanceGeneProtExon1 <- read_csv("../data/github/lncrna-exon1-negative-control-dataset-features.csv")
distanceGeneProtExon2 <- read_csv("../data/github/lncrna-exon2-negative-control-dataset-features.csv")
distanceGeneProt <- rbind(distanceGeneProtExon1,distanceGeneProtExon2)
distanceGeneProt$group <- "mRNA(-)"
###########################
## Generate histogram plots
ggplot(prot_neg_zscores, aes(x = DistanceGene, fill = group, colour = group)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 500) +
  labs(title = "mRNA", x = "Distane to Gene", y = "Frequency") +
  scale_fill_manual(values = c("blue")) +
  scale_color_manual(values = c("blue")) +
  #scale_x_log10(breaks = scales::log_breaks(n = 10)) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000, 100000,1000000,10000000),
                     labels = scales::label_number()
  )

ggplot(sncrna_neg_zscores, aes(x = distance, fill = group, colour = group)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 500) +
  labs(title = "sncRNA", x = "Distane to Gene", y = "Frequency") +
  scale_fill_manual(values = c("blue")) +
  scale_color_manual(values = c("blue")) +
  #scale_x_log10(breaks = scales::log_breaks(n = 10)) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000, 100000,1000000,10000000),
                     labels = scales::label_number()
  )

ggplot(lncrna_neg_zscores, aes(x = distance, fill = group, colour = group)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 500) +
  labs(title = "lncRNA", x = "Distane to Gene", y = "Frequency") +
  scale_fill_manual(values = c("blue")) +
  scale_color_manual(values = c("blue")) +
  #scale_x_log10(breaks = scales::log_breaks(n = 10)) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000, 100000,1000000,10000000),
                     labels = scales::label_number()
  )

##############################
### USE heatscatter R PACKAGE.
#reset plot
dev.off()
layout(1)

# For high-res PNG
png("effectOnDistanceJoined5k_5M.png", width=4920, height=3240, res=300)
nf <- layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
#GC%
par(mar = c(7, 5, 7, 2))  # Increase top margin (third number)
heatscatter(prot_neg_zscores$distance, 
            prot_neg_zscores$GC.content,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "mRNA", cex.main = 3)

# Filter out NA values first
valid_data <- !is.na(prot_neg_zscores$distance) & !is.na(prot_neg_zscores$GC.content)

# Sort the valid data
sorted_idx <- order(prot_neg_zscores$distance[valid_data])
sorted_distance <- prot_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- prot_neg_zscores$GC.content[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)

par(mar = c(7, 5, 7, 2))  # Increase top margin (third number)
heatscatter(sncrna_neg_zscores$distance, 
            sncrna_neg_zscores$GC.content,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "sncRNA", cex.main = 3)

# Filter out NA values first
valid_data <- !is.na(sncrna_neg_zscores$distance) & !is.na(sncrna_neg_zscores$GC.content)

# Sort the valid data
sorted_idx <- order(sncrna_neg_zscores$distance[valid_data])
sorted_distance <- sncrna_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- sncrna_neg_zscores$GC.content[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)



heatscatter(lncrna_neg_zscores$distance, 
            lncrna_neg_zscores$GC.content,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "lncRNA", cex.main = 3)

# Filter out NA values first
valid_data <- !is.na(lncrna_neg_zscores$distance) & !is.na(lncrna_neg_zscores$GC.content)

# Sort the valid data
sorted_idx <- order(lncrna_neg_zscores$distance[valid_data])
sorted_distance <- lncrna_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- lncrna_neg_zscores$GC.content[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)


# SNP density
heatscatter(prot_neg_zscores$distance, 
            prot_neg_zscores$SNP.density,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "mRNA", cex.main = 3)

# Filter out NA values first
valid_data <- !is.na(prot_neg_zscores$distance) & !is.na(prot_neg_zscores$SNP.density)

# Sort the valid data
sorted_idx <- order(prot_neg_zscores$distance[valid_data])
sorted_distance <- prot_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- prot_neg_zscores$SNP.density[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)

heatscatter(sncrna_neg_zscores$distance, 
            sncrna_neg_zscores$SNP.density,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "sncRNA", cex.main = 3)

# Filter out NA values first
valid_data <- !is.na(sncrna_neg_zscores$distance) & !is.na(sncrna_neg_zscores$SNP.density)

# Sort the valid data
sorted_idx <- order(sncrna_neg_zscores$distance[valid_data])
sorted_distance <- sncrna_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- sncrna_neg_zscores$SNP.density[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)

heatscatter(lncrna_neg_zscores$distance, 
            lncrna_neg_zscores$SNP.density,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "lncRNA", cex.main = 3)
# Filter out NA values first
valid_data <- !is.na(lncrna_neg_zscores$distance) & !is.na(lncrna_neg_zscores$SNP.density)

# Sort the valid data
sorted_idx <- order(lncrna_neg_zscores$distance[valid_data])
sorted_distance <- lncrna_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- lncrna_neg_zscores$SNP.density[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)
dev.off()


# JOIN PLOTS in ONE using patchwork package
gc_plot <- (wrap_elements((prot_gc_plot + sncrna_gc_plot + lncrna_gc_plot) + 
                                    plot_annotation(title = "GC%",
                                                    tag_levels = list(c('A','B','C'))) +
                                    plot_layout(axis_titles = "collect") &
                                    theme(plot.title = element_text(size = 56, hjust = 0.5),
                                          plot.tag.position = c(0, 1),
                                          plot.tag = element_text(size = 38, face = "bold", hjust = 0, vjust = 0, margin = margin(40,40,40,40)))))

snp_plot <- (wrap_elements((prot_snp_plot + sncrna_snp_plot + lncrna_snp_plot) + 
                             plot_annotation(title = "SNPs",
                                             tag_levels = list(c('D','E','F'))) + 
                             plot_layout(axis_titles = "collect") &
                             theme(plot.title = element_text(size = 56, hjust = 0.5),
                                   plot.tag.position = c(0, 1),
                                   plot.tag = element_text(size = 38, face = "bold", hjust = 0, vjust = 0, margin = margin(40,40,40,40)))))

joined_plot <- gc_plot / snp_plot
ggsave(paste0("distance_effect_loess_wiggle_veplus",'_',"joined",'.png'), path = "../results/latest1000all/logScatterPlots/Paper/", joined_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

#################################
# Spearman Correlation Analysis #
# Function to compute a Spearman correlation coefficient for the specified features
compute_correlation <- function(zscores, selected_features) {
  corr_matrix_list <- list()
  for(feature in selected_features) {
    print(feature)
    
    corr_obj <- rcorr(zscores$DistanceGene,
                      zscores[[feature]], 
                      type = "spearman")
    
    corr_matrix_list[[feature]] <- list(rho = corr_obj[["r"]][2,1], p.value = corr_obj[["P"]][2,1])
  }
  return(corr_matrix_list)
}
help("rcorr")
# Function to Correct p-values
adjust_pvalues <- function(corr_matrix, selected_features) {
  pvalues <- list()
  for (feature in selected_features) {
    pvalues[[feature]] <- corr_matrix[[feature]]$p.value
  }
  pvalues_adj <- p.adjust(pvalues, method="BH", n = length(pvalues) * 3)
  return(pvalues_adj)
}

# Function to Extract Spearman's rhos
extract_rhos <- function(corr_matrix, selected_features) {
  rhos <- list()
  for (feature in selected_features) {
    rhos[[feature]] <- corr_matrix[[feature]]$rho
  }
  return(rhos)
}
help("p.adjust")
prot_zscores <- list_with_zscores_protein_neg
prot_zscores <- rename_features(prot_zscores, eod_select_features, eod_select_features_labels)
sncrna_zscores <- list_with_zscores_sncrna_neg
sncrna_zscores <- rename_features(sncrna_zscores, eod_select_features, eod_select_features_labels)
lncrna_zscores <- list_with_zscores_lncrna_neg
lncrna_zscores <- rename_features(lncrna_zscores, eod_select_features, eod_select_features_labels)

# Add distanceGene feature
prot_zscores$DistanceGene <- rbind(protExon2NCData,protExon3NCData)$DistanceGene
sncrna_zscores$DistanceGene <- sncrnaNCData$DistanceGene
lncrna_zscores$DistanceGene <- rbind(lncrnaExon1NCData,lncrnaExon2NCData)$DistanceGene
# Compute correlation coefficient to reveal distance effect:
#Protein coding
prot_corr_matrix <- compute_correlation(prot_zscores, eod_select_features_labels)
prot_rhos <- as.data.frame(extract_rhos(prot_corr_matrix, eod_select_features_labels), check.names = FALSE)
rownames(prot_rhos) <- "rho"
prot_pvalues_adj <- as.data.frame(t(adjust_pvalues(prot_corr_matrix, eod_select_features_labels)))
rownames(prot_pvalues_adj) <- "p.value"
prot_corr_matrix <- rbind(prot_rhos, prot_pvalues_adj)

#sncRNA
sncrna_corr_matrix <- compute_correlation(sncrna_zscores, eod_select_features_labels)
sncrna_rhos <- as.data.frame(extract_rhos(sncrna_corr_matrix, eod_select_features_labels), check.names = FALSE)
rownames(sncrna_rhos) <- "rho"
sncrna_pvalues_adj <- as.data.frame(t(adjust_pvalues(sncrna_corr_matrix, eod_select_features_labels)))
rownames(sncrna_pvalues_adj) <- "p.value"
sncrna_corr_matrix <- rbind(sncrna_rhos, sncrna_pvalues_adj)

#lncRNA
lncrna_corr_matrix <- compute_correlation(lncrna_zscores, eod_select_features_labels)
lncrna_rhos <- as.data.frame(extract_rhos(lncrna_corr_matrix, eod_select_features_labels), check.names = FALSE)
rownames(lncrna_rhos) <- "rho"
lncrna_pvalues_adj <- as.data.frame(t(adjust_pvalues(lncrna_corr_matrix, eod_select_features_labels)))
rownames(lncrna_pvalues_adj) <- "p.value"
lncrna_corr_matrix <- rbind(lncrna_rhos, lncrna_pvalues_adj)

prot_corr_matrix$measurement <- row.names(prot_corr_matrix)
sncrna_corr_matrix$measurement <- row.names(sncrna_corr_matrix)
lncrna_corr_matrix$measurement <- row.names(lncrna_corr_matrix)

prot_corr_matrix <- prot_corr_matrix |>
  pivot_longer(cols = "GC%":"MAF",
               names_to = "feature")

prot_corr_matrix <- prot_corr_matrix |>
  pivot_wider(names_from = measurement,
              values_from = value)

sncrna_corr_matrix <- sncrna_corr_matrix |>
  pivot_longer(cols = "GC%":"MAF",
               names_to = "feature")

sncrna_corr_matrix <- sncrna_corr_matrix |>
  pivot_wider(names_from = measurement,
              values_from = value)

lncrna_corr_matrix <- lncrna_corr_matrix |>
  pivot_longer(cols = "GC%":"MAF",
               names_to = "feature")

lncrna_corr_matrix <- lncrna_corr_matrix |>
  pivot_wider(names_from = measurement,
              values_from = value)

write.csv(prot_corr_matrix, file = "../results/latest1000all/correlation_analysis/corr_matrix_protein.csv")
write.csv(sncrna_corr_matrix, file = "../results/latest1000all/correlation_analysis/corr_matrix_sncrna.csv")
write.csv(lncrna_corr_matrix, file = "../results/latest1000all/correlation_analysis/corr_matrix_lncrna.csv")


# function to calculate a custom K-S analysis to keep sign
signed_ks_test <- function(x, y) {
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

# Function to perform K-S analysis
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
    
    ks_test <- signed_ks_test(negative_col, positive_col)
    ks_test_p <- ks.test(negative_col, positive_col)
    
    results[1,i] <- ks_test$signed_D
    results[2,i] <- ks_test$max_diff
    results[3,i] <- ks_test$min_diff
    results[4,i] <- ks_test_p$p.value
  }
  return(results)
}
# Compute KS between first and last bin (don't call it ks-"test" unless returning p.value)
# Protein coding
ks_results_prot_bin1_4 <- run_ks_tests(prot_delta_zscores_bin1, 
                                       prot_delta_zscores_bin4, 
                                       eod_select_features_labels)
# sncRNA
ks_results_sncrna_bin1_4 <- run_ks_tests(sncrna_delta_zscores_bin1, 
                                         sncrna_delta_zscores_bin4, 
                                         eod_select_features_labels)
# lncRNA
ks_results_lncrna_bin1_4 <- run_ks_tests(lncrna_delta_zscores_bin1, 
                                         lncrna_delta_zscores_bin4, 
                                         eod_select_features_labels)

write.csv(ks_results_prot_bin1_4, file = "../results/latest1000all/ks_analysis/ks_protein_bin1_vs_bin4.csv")
write.csv(ks_results_sncrna_bin1_4, file = "../results/latest1000all/ks_analysis/ks_sncrna_bin1_vs_bin4.csv")
write.csv(ks_results_lncrna_bin1_4, file = "../results/latest1000all/ks_analysis/ks_lncrna_bin1_vs_bin4.csv")


################
