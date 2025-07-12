# Script to load data into R
library(tidyverse)

################################
# Load features in a dataframe #
################################
feature_matrix <- read.csv("../data/features/gene_functionality_features_latest1000all_extra.csv",
                           header = TRUE,
                           check.names = FALSE)
unique(feature_matrix$Dataset)
##############################
# Retrieve genetype matrices #
##############################
funcProtExon2Data <- feature_matrix[feature_matrix$Dataset=="protein-coding-exon2",]
funcProtExon3Data <- feature_matrix[feature_matrix$Dataset=="protein-coding-exon3",]
funcLncrnaExon1Data <- feature_matrix[feature_matrix$Dataset=="lncrna-exon1",]
funcLncrnaExon2Data <- feature_matrix[feature_matrix$Dataset=="lncrna-exon2",]
funcSncrnaDataset <- feature_matrix[feature_matrix$Dataset=="short-ncrna",]

protExon2NCData <- feature_matrix[feature_matrix$Dataset=="protein-exon2-negative-control",]
protExon3NCData <- feature_matrix[feature_matrix$Dataset=="protein-exon3-negative-control",]
lncrnaExon1NCData <- feature_matrix[feature_matrix$Dataset=="lncrna-exon1-negative-control",]
lncrnaExon2NCData <- feature_matrix[feature_matrix$Dataset=="lncrna-exon2-negative-control",]
sncrnaNCData <- feature_matrix[feature_matrix$Dataset=="short-ncrna-negative-control",]


# Function to transform features to numeric datatype
convert_to_numeric <- function(column) {
  as.numeric(as.character(column))
}

# Function to get a numeric version of the received dataframe
get_numeric_dataframe <- function(dataframe, omit_nas = FALSE) {
  # Store functional and dataset features in independent vectors
  functional_column <- dataframe$Functional
  dataset_column <- dataframe$Dataset
  # Specify features to keep in dataframe
  select_features <- c("Random number",
                       "GC content",
                       "AA","AC","AG","AT","CA","CC","CpG","CT","GA","GC","GG","GT","TA","TC","TG","TT",
                       "low_complexity_density","Conservation","phyloP max_100w","GERP_91_mammals_max","GERP_63_amniotes_max",
                       "Expression","RPKM_primary cell","Copy number","Repeat free","Fickett_score","RNAcode","RNA Covariance",
                       "MFE","accessibility","RNAalifold","RNA-RNA interactions","SNP density","gnomAD_MAF",
                       "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome")
  # Initialize dataframe to return
  # Select only specified features
  selected_features_matrix <- dataframe[, select_features]
  # If a dataframe without NAs is requested
  if(omit_nas) {
    dataframe_no_nas <- dataframe[, c(select_features,"Functional","Dataset")]
    # Remove NAs
    dataframe_no_nas <- na.omit(dataframe_no_nas)
    # Store functional and dataset features in separate vectors
    functional_column <- dataframe_no_nas$Functional
    dataset_column <- dataframe_no_nas$Dataset
    # Select only specified features replacing with no nas version
    selected_features_matrix <- dataframe_no_nas[, select_features]
  }
  
  # Convert to numeric
  feature_matrix_numeric <- as.data.frame(sapply(selected_features_matrix, convert_to_numeric))
  # Restore functional and dataset columns
  feature_matrix_numeric$Functional <- functional_column
  feature_matrix_numeric$Dataset <- dataset_column
  
  return(feature_matrix_numeric)
}

#feature_matrix_numeric_no_nas <- get_numeric_dataframe(feature_matrix, omit_nas = TRUE)
#feature_matrix_numeric <- get_numeric_dataframe(feature_matrix)

#feature_matrix_numeric[feature_matrix_numeric$Dataset == "lncrna-exon1",]$accessibility <- functional_lncrna_exon1_accessibility_feature$accessibility
#feature_matrix_numeric[feature_matrix_numeric$Dataset == "lncrna-exon2",]$accessibility <- functional_lncrna_exon2_accessibility_feature$accessibility
#feature_matrix_numeric[feature_matrix_numeric$Dataset == "lncrna-exon1-negative-control",]$accessibility <- lncrna_exon1_nc_accessibility_feature$accessibility
#feature_matrix_numeric[feature_matrix_numeric$Dataset == "lncrna-exon2-negative-control",]$accessibility <- lncrna_exon2_nc_accessibility_feature$accessibility

