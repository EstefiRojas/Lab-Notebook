library(dplyr)
install.packages("corrr")
library('corrr')
install.packages("ggcorrplot")
library(ggcorrplot)
install.packages("FactoMineR")
library("FactoMineR")
install.packages("tidyr")
library(tidyr)
library(reshape2)
library(factoextra)
install.packages("random")
library(random)
install.packages("pwr")
library(pwr)
install.packages("Hmisc")
library(Hmisc)
install.packages("robustbase")
library(robustbase) # For the Sn scale estimator


# Load Data #
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

# Rename GeneID column from functional #
names(funcProtExon2Data)[names(funcProtExon2Data) == 'GeneID'] <- 'Distance'
names(funcProtExon3Data)[names(funcProtExon3Data) == 'GeneID'] <- 'Distance'
names(funcLncrnaExon1Data)[names(funcLncrnaExon1Data) == 'GeneID'] <- 'Distance'
names(funcLncrnaExon2Data)[names(funcLncrnaExon2Data) == 'GeneID'] <- 'Distance'
names(funcSncrnaDataset)[names(funcSncrnaDataset) == 'GeneID'] <- 'Distance'

names(protExon2NCData)[names(protExon2NCData) == 'DistanceGene'] <- 'Distance'
names(protExon3NCData)[names(protExon3NCData) == 'DistanceGene'] <- 'Distance'
names(lncrnaExon1NCData)[names(lncrnaExon1NCData) == 'DistanceGene'] <- 'Distance'
names(lncrnaExon2NCData)[names(lncrnaExon2NCData) == 'DistanceGene'] <- 'Distance'
names(sncrnaNCData)[names(sncrnaNCData) == 'DistanceGene'] <- 'Distance'

as.data.frame(funcProtExon2Data[["Distance"]], nm="Distance")

funcProtExon2Data$Distance <- 0
funcProtExon3Data[["Distance"]] <- 0
funcLncrnaExon1Data[["Distance"]] <- 0
funcLncrnaExon2Data[["Distance"]] <- 0
funcSncrnaDataset[["Distance"]] <- 0


pos_distance_column <- as.data.frame(rbind(as.data.frame(funcProtExon2Data[["Distance"]], nm="Distance"),
                                           as.data.frame(funcProtExon3Data[["Distance"]], nm="Distance"),
                                           as.data.frame(funcLncrnaExon1Data[["Distance"]], nm="Distance"),
                                           as.data.frame(funcLncrnaExon2Data[["Distance"]], nm="Distance"),
                                           as.data.frame(funcSncrnaDataset[["Distance"]], nm="Distance")))


neg_distance_column <- as.data.frame(rbind(as.data.frame(protExon2NCData[["Distance"]], nm="Distance"),
                                                         as.data.frame(protExon3NCData[["Distance"]], nm="Distance"),
                                                         as.data.frame(lncrnaExon1NCData[["Distance"]], nm="Distance"),
                                                         as.data.frame(lncrnaExon2NCData[["Distance"]], nm="Distance"),
                                                         as.data.frame(sncrnaNCData[["Distance"]], nm="Distance")))

distance_column <- rbind(pos_distance_column,neg_distance_column)

funcProtExon2Data <- funcProtExon2Data[, -which(names(funcProtExon2Data) == "GeneID")]
funcProtExon3Data <- funcProtExon3Data[, -which(names(funcProtExon3Data) == "GeneID")]
funcLncrnaExon1Data <- funcLncrnaExon1Data[, -which(names(funcLncrnaExon1Data) == "GeneID")]
funcLncrnaExon2Data <- funcLncrnaExon2Data[, -which(names(funcLncrnaExon2Data) == "GeneID")]
funcSncrnaDataset <- funcSncrnaDataset[, -which(names(funcSncrnaDataset) == "GeneID")]

# Drop Distance column from negative controls #
protExon2NCData <- protExon2NCData[, -which(names(protExon2NCData) == "DistanceGene")]
protExon3NCData <- protExon3NCData[, -which(names(protExon3NCData) == "DistanceGene")]
lncrnaExon1NCData <- lncrnaExon1NCData[, -which(names(lncrnaExon1NCData) == "DistanceGene")]
lncrnaExon2NCData <- lncrnaExon2NCData[, -which(names(lncrnaExon2NCData) == "DistanceGene")]
sncrnaNCData <- sncrnaNCData[, -which(names(sncrnaNCData) == "DistanceGene")]



# Load dinucleotide frequencies feature
din_funcProtExon2Data <- data.frame(read.csv("../data/datasets/dinucleotide_feature/protein-exon2-dinucleotide-feature.csv", header=TRUE))
din_funcProtExon3Data <- data.frame(read.csv("../data/datasets/dinucleotide_feature/protein-exon3-dinucleotide-feature.csv", header=TRUE))
din_funcLncrnaExon1Data <- data.frame(read.csv("../data/datasets/dinucleotide_feature/lncrna-exon1-dinucleotide-feature.csv", header=TRUE))
din_funcLncrnaExon2Data <- data.frame(read.csv("../data/datasets/dinucleotide_feature/lncrna-exon2-dinucleotide-feature.csv", header=TRUE))
din_funcSncrnaDataset <- data.frame(read.csv("../data/datasets/dinucleotide_feature/short-ncrna-dinucleotide-feature.csv", header=TRUE))

din_protExon2NCData <- data.frame(read.csv("../data/datasets/dinucleotide_feature/protein-exon2-NC-dinucleotide-feature.csv", header=TRUE))
din_protExon3NCData <- data.frame(read.csv("../data/datasets/dinucleotide_feature/protein-exon3-NC-dinucleotide-feature.csv", header=TRUE))
din_lncrnaExon1NCData <- data.frame(read.csv("../data/datasets/dinucleotide_feature/lncrna-exon1-NC-dinucleotide-feature.csv", header=TRUE))
din_lncrnaExon2NCData <- data.frame(read.csv("../data/datasets/dinucleotide_feature/lncrna-exon2-NC-dinucleotide-feature.csv", header=TRUE))
din_sncrnaNCData <- data.frame(read.csv("../data/datasets/dinucleotide_feature/short-ncrna-NC-dinucleotide-feature.csv", header=TRUE))

# Replace old dinucleotide features
all_din <- rbind(din_funcProtExon2Data,din_funcProtExon3Data,din_funcLncrnaExon1Data,din_funcLncrnaExon2Data,din_funcSncrnaDataset,
                 din_protExon2NCData,din_protExon3NCData,din_lncrnaExon1NCData,din_lncrnaExon2NCData,din_sncrnaNCData)
feature_matrix$AA <- all_din$AA
feature_matrix$AC <- all_din$AC
feature_matrix$AG <- all_din$AG
feature_matrix$AT <- all_din$AT

feature_matrix$CA <- all_din$CA
feature_matrix$CC <- all_din$CC
feature_matrix$CpG <- all_din$CpG
feature_matrix$CT <- all_din$CT

feature_matrix$GA <- all_din$GA
feature_matrix$GC <- all_din$GC
feature_matrix$GG <- all_din$GG
feature_matrix$GT <- all_din$GT

feature_matrix$TA <- all_din$TA
feature_matrix$TC <- all_din$TC
feature_matrix$TG <- all_din$TG
feature_matrix$TT <- all_din$TT
# Drop Functional column from dinucleotide features #
#din_funcProtExon2Data <- din_funcProtExon2Data[, -which(names(din_funcProtExon2Data) == "Functional")]
#din_funcProtExon3Data <- din_funcProtExon3Data[, -which(names(din_funcProtExon3Data) == "Functional")]
#din_funcLncrnaExon1Data <- din_funcLncrnaExon1Data[, -which(names(din_funcLncrnaExon1Data) == "Functional")]
#din_funcLncrnaExon2Data <- din_funcLncrnaExon2Data[, -which(names(din_funcLncrnaExon2Data) == "Functional")]
#din_funcSncrnaDataset <- din_funcSncrnaDataset[, -which(names(din_funcSncrnaDataset) == "Functional")]
#din_protExon2NCData <- din_protExon2NCData[, -which(names(din_protExon2NCData) == "Functional")]
#din_protExon3NCData <- din_protExon3NCData[, -which(names(din_protExon3NCData) == "Functional")]
#din_lncrnaExon1NCData <- din_lncrnaExon1NCData[, -which(names(din_lncrnaExon1NCData) == "Functional")]
#din_lncrnaExon2NCData <- din_lncrnaExon2NCData[, -which(names(din_lncrnaExon2NCData) == "Functional")]
#din_sncrnaNCData <- din_sncrnaNCData[, -which(names(din_sncrnaNCData) == "Functional")]

# Join dinucleotide features to main features #
#funcProtExon2Data <- cbind(funcProtExon2Data, din_funcProtExon2Data)
#funcProtExon3Data <- cbind(funcProtExon3Data, din_funcProtExon3Data)
#funcLncrnaExon1Data <- cbind(funcLncrnaExon1Data, din_funcLncrnaExon1Data)
#funcLncrnaExon2Data <- cbind(funcLncrnaExon2Data, din_funcLncrnaExon2Data)
#funcSncrnaDataset <- cbind(funcSncrnaDataset, din_funcSncrnaDataset)
#protExon2NCData <- cbind(protExon2NCData, din_protExon2NCData)
#protExon3NCData <- cbind(protExon3NCData, din_protExon3NCData)
#lncrnaExon1NCData <- cbind(lncrnaExon1NCData, din_lncrnaExon1NCData)
#lncrnaExon2NCData <- cbind(lncrnaExon2NCData, din_lncrnaExon2NCData)
#sncrnaNCData <- cbind(sncrnaNCData, din_sncrnaNCData)



# Join data sets #
matrix_protein <- rbind(funcProtExon2Data, funcProtExon3Data)
matrix_protein[matrix_protein == "."] <- NA  # Change . to NA for missing values (GERP columns)
summary(matrix_protein)# Get the column names of both datasets
matrix_lncrna <- rbind(funcLncrnaExon1Data, funcLncrnaExon2Data)
matrix_lncrna[matrix_lncrna == "."] <- NA  # Change . to NA for missing values (GERP columns)
summary(matrix_lncrna)
matrix_short_ncrna <- funcSncrnaDataset
matrix_short_ncrna[matrix_short_ncrna == "."] <- NA  # Change . to NA for missing values (GERP columns)
summary(matrix_short_ncrna)
matrix_negative_control <- rbind(protExon2NCData,protExon3NCData,lncrnaExon1NCData,lncrnaExon2NCData,sncrnaNCData)
matrix_negative_control[matrix_negative_control == "."] <- NA  # Change . to NA for missing values (GERP columns)
summary(matrix_negative_control)


# Create negative control dataframes for aech gene type
matrix_protein_nc <- rbind(protExon2NCData, protExon3NCData)
matrix_protein_nc[matrix_protein_nc == "."] <- NA  # Change . to NA for missing values (GERP columns)
summary(matrix_protein_nc)
matrix_lncrna_nc <- rbind(lncrnaExon1NCData, lncrnaExon2NCData)
matrix_lncrna_nc[matrix_lncrna_nc == "."] <- NA  # Change . to NA for missing values (GERP columns)
summary(matrix_lncrna_nc)
matrix_short_ncrna_nc <- sncrnaNCData
matrix_short_ncrna_nc[matrix_short_ncrna_nc == "."] <- NA  # Change . to NA for missing values (GERP columns)
summary(matrix_short_ncrna_nc)

#############################
# Add random number feature #
#############################
#Protein Coding
random::randomQuota()
true_random_numbers <- random::randomNumbers(n = nrow(matrix_protein), min = 0, max = 1000, col = 1)
matrix_protein$random <- true_random_numbers[,1]
true_random_numbers <- random::randomNumbers(n = 10000, min = 0, max = 1000, col = 1)
true_random_numbers1 <- random::randomNumbers(n = 9250, min = 0, max = 1000, col = 1)
matrix_protein_nc$random <- c(true_random_numbers[,1],true_random_numbers1[,1])

#lncRNA
true_random_numbers <- random::randomNumbers(n = nrow(matrix_lncrna), min = 0, max = 1000, col = 1)
matrix_lncrna$random <- true_random_numbers[,1]
true_random_numbers <- random::randomNumbers(n = 10000, min = 0, max = 1000, col = 1)
true_random_numbers1 <- random::randomNumbers(n = 9608, min = 0, max = 1000, col = 1)
matrix_lncrna_nc$random <- c(true_random_numbers[,1],true_random_numbers1)

#Short ncRNA
true_random_numbers <- random::randomNumbers(n = nrow(matrix_short_ncrna), min = 0, max = 1000, col = 1)
matrix_short_ncrna$random <- true_random_numbers[,1]
true_random_numbers <- random::randomNumbers(n = nrow(matrix_short_ncrna_nc), min = 0, max = 1000, col = 1)
matrix_short_ncrna_nc$random <- true_random_numbers[,1]

#Negatives only
true_random_numbers <- random::randomNumbers(n = nrow(matrix_negative_control), min = 0, max = 1000, col = 1)
matrix_negative_control$random <- true_random_numbers[,1]

#########################################
# Rename columns GC. to GC_%, CG to CpG, etc... #
#########################################
names(matrix_protein)[names(matrix_protein) == 'GC.'] <- 'GC_percent'
names(matrix_protein)[names(matrix_protein) == 'CG'] <- 'CpG'
names(matrix_protein)[names(matrix_protein) == 'phyloP_meam_100w'] <- 'phyloP_mean_100w'
names(matrix_protein)[names(matrix_protein) == 'coding_potential'] <- 'RNAcode'
names(matrix_protein)[names(matrix_protein) == 'Neutral_predictor'] <- 'Random number'
names(matrix_protein)[names(matrix_protein) == 'phyloP_max'] <- 'phyloP max_241w'
names(matrix_protein)[names(matrix_protein) == 'repeat_distance'] <- 'Dfam_sum'
names(matrix_protein)[names(matrix_protein) == 'SNP_density'] <- 'gnomAD_SNP_density'
names(matrix_protein)[names(matrix_protein) == 'MAF_avg'] <- 'gnomAD_MAF'


names(matrix_protein_nc)[names(matrix_protein_nc) == 'GC.'] <- 'GC_percent'
names(matrix_protein_nc)[names(matrix_protein_nc) == 'CG'] <- 'CpG'
names(matrix_protein_nc)[names(matrix_protein_nc) == 'phyloP_meam_100w'] <- 'phyloP_mean_100w'
names(matrix_protein_nc)[names(matrix_protein_nc) == 'coding_potential'] <- 'RNAcode'
names(matrix_protein_nc)[names(matrix_protein_nc) == 'Neutral_predictor'] <- 'Random number'
names(matrix_protein_nc)[names(matrix_protein_nc) == 'phyloP_max'] <- 'phyloP max_241w'
names(matrix_protein_nc)[names(matrix_protein_nc) == 'repeat_distance'] <- 'Dfam_sum'
names(matrix_protein_nc)[names(matrix_protein_nc) == 'SNP_density'] <- 'gnomAD_SNP_density'
names(matrix_protein_nc)[names(matrix_protein_nc) == 'MAF_avg'] <- 'gnomAD_MAF'

names(matrix_lncrna)[names(matrix_lncrna) == 'GC.'] <- 'GC_percent'
names(matrix_lncrna)[names(matrix_lncrna) == 'CG'] <- 'CpG'
names(matrix_lncrna)[names(matrix_lncrna) == 'phyloP_meam_100w'] <- 'phyloP_mean_100w'
names(matrix_lncrna)[names(matrix_lncrna) == 'coding_potential'] <- 'RNAcode'
names(matrix_lncrna)[names(matrix_lncrna) == 'Neutral_predictor'] <- 'Random number'
names(matrix_lncrna)[names(matrix_lncrna) == 'phyloP_max'] <- 'phyloP max_241w'
names(matrix_lncrna)[names(matrix_lncrna) == 'repeat_distance'] <- 'Dfam_sum'
names(matrix_lncrna)[names(matrix_lncrna) == 'SNP_density'] <- 'gnomAD_SNP_density'
names(matrix_lncrna)[names(matrix_lncrna) == 'MAF_avg'] <- 'gnomAD_MAF'

names(matrix_lncrna_nc)[names(matrix_lncrna_nc) == 'GC.'] <- 'GC_percent'
names(matrix_lncrna_nc)[names(matrix_lncrna_nc) == 'CG'] <- 'CpG'
names(matrix_lncrna_nc)[names(matrix_lncrna_nc) == 'phyloP_meam_100w'] <- 'phyloP_mean_100w'
names(matrix_lncrna_nc)[names(matrix_lncrna_nc) == 'coding_potential'] <- 'RNAcode'
names(matrix_lncrna_nc)[names(matrix_lncrna_nc) == 'Neutral_predictor'] <- 'Random number'
names(matrix_lncrna_nc)[names(matrix_lncrna_nc) == 'phyloP_max'] <- 'phyloP max_241w'
names(matrix_lncrna_nc)[names(matrix_lncrna_nc) == 'repeat_distance'] <- 'Dfam_sum'
names(matrix_lncrna_nc)[names(matrix_lncrna_nc) == 'SNP_density'] <- 'gnomAD_SNP_density'
names(matrix_lncrna_nc)[names(matrix_lncrna_nc) == 'MAF_avg'] <- 'gnomAD_MAF'

names(matrix_short_ncrna)[names(matrix_short_ncrna) == 'GC.'] <- 'GC_percent'
names(matrix_short_ncrna)[names(matrix_short_ncrna) == 'CG'] <- 'CpG'
names(matrix_short_ncrna)[names(matrix_short_ncrna) == 'phyloP_meam_100w'] <- 'phyloP_mean_100w'
names(matrix_short_ncrna)[names(matrix_short_ncrna) == 'coding_potential'] <- 'RNAcode'
names(matrix_short_ncrna)[names(matrix_short_ncrna) == 'Neutral_predictor'] <- 'Random number'
names(matrix_short_ncrna)[names(matrix_short_ncrna) == 'phyloP_max'] <- 'phyloP max_241w'
names(matrix_short_ncrna)[names(matrix_short_ncrna) == 'repeat_distance'] <- 'Dfam_sum'
names(matrix_short_ncrna)[names(matrix_short_ncrna) == 'SNP_density'] <- 'gnomAD_SNP_density'
names(matrix_short_ncrna)[names(matrix_short_ncrna) == 'MAF_avg'] <- 'gnomAD_MAF'

names(matrix_short_ncrna_nc)[names(matrix_short_ncrna_nc) == 'GC.'] <- 'GC_percent'
names(matrix_short_ncrna_nc)[names(matrix_short_ncrna_nc) == 'CG'] <- 'CpG'
names(matrix_short_ncrna_nc)[names(matrix_short_ncrna_nc) == 'phyloP_meam_100w'] <- 'phyloP_mean_100w'
names(matrix_short_ncrna_nc)[names(matrix_short_ncrna_nc) == 'coding_potential'] <- 'RNAcode'
names(matrix_short_ncrna_nc)[names(matrix_short_ncrna_nc) == 'Neutral_predictor'] <- 'Random number'
names(matrix_short_ncrna_nc)[names(matrix_short_ncrna_nc) == 'phyloP_max'] <- 'phyloP max_241w'
names(matrix_short_ncrna_nc)[names(matrix_short_ncrna_nc) == 'repeat_distance'] <- 'Dfam_sum'
names(matrix_short_ncrna_nc)[names(matrix_short_ncrna_nc) == 'SNP_density'] <- 'gnomAD_SNP_density'
names(matrix_short_ncrna_nc)[names(matrix_short_ncrna_nc) == 'MAF_avg'] <- 'gnomAD_MAF'


matrix_negative_control <- rbind(protExon2NCData,protExon3NCData,lncrnaExon1NCData,lncrnaExon2NCData,sncrnaNCData)

#Join all
all_matrix <- rbind(matrix_protein, matrix_lncrna, matrix_short_ncrna, matrix_protein_nc, matrix_lncrna_nc, matrix_short_ncrna_nc)

# Check for nulls
colSums(is.na(all_matrix))
# Remove nulls
#all_matrix <- na.omit(all_matrix)

# reload matrix
feature_matrix <- read.csv("../data/features/gene_functionality_features_latest1000all.csv",header = TRUE, check.names = FALSE)
#names(feature_matrix)[names(feature_matrix) == 'Neutral_predictor'] <- 'Random_number'
# Keep only variables of interest
functional_column <- feature_matrix$Functional
dataset_column <- feature_matrix$Dataset
select_features_old_names <- c("Random number",
                     "GC content",
                     "CpG","GA","GG","TA","TC",
                     "low_complexity_density","phyloP max_241w","phyloP max_100w",
                     "RPKM_tissue","RPKM_primary cell","Copy number","Repeat free","RNAcode","Max covariance",
                     "MFE","RNAalifold","Interaction_ave","gnomAD_SNP_density","gnomAD_MAF",
                     "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome")
select_features <- c("Random number",
                     "GC content",
                     "CpG","GA","GG","TA","TC",
                     "low_complexity_density","Conservation","phyloP max_100w",
                     "Expression","RPKM_primary cell","Copy number","Repeat free","RNAcode","RNA Covariance",
                     "MFE","RNAalifold","RNA-RNA interactions","SNP density","gnomAD_MAF",
                     "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome")
select_features <- c("Random number",
                     "GC content",
                     "AA","AC","AG","AT","CA","CC","CpG","CT","GA","GC","GG","GT","TA","TC","TG","TT",
                     "low_complexity_density","Conservation","phyloP max_100w","GERP_91_mammals_max","GERP_63_amniotes_max",
                     "Expression","RPKM_primary cell","Copy number","Repeat free","Fickett_score","RNAcode","RNA Covariance",
                     "MFE","accessibility","RNAalifold","RNA-RNA interactions","SNP density","gnomAD_MAF",
                     "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome")
cols1 <- colnames(feature_matrix)
#excludeCols <- c("Dataset","ID","Functional","Chromosome","Sequence","Start","End","GERP_63_amniotes_mean","GERP_91_mammals_mean","phyloP_mean_100w")
select_features_missing <- setdiff(cols1, select_features)
select_features_missing
#all_matrix <- all_matrix[, select_features] # Order columns
#feature_matrix <- feature_matrix[, names(feature_matrix) %in% select_features]
feature_matrix_selected <- feature_matrix[, select_features]
str(feature_matrix_selected)

convert_to_numeric <- function(column) {
  as.numeric(as.character(column))
}

feature_matrix_numeric <- as.data.frame(sapply(feature_matrix_numeric, convert_to_numeric)) # Matrix for PCA by Sequence
summary(feature_matrix_numeric)

feature_matrix_numeric_no_nas <- na.omit(feature_matrix_numeric)
feature_matrix_numeric_t <- t(feature_matrix_numeric_no_nas) # Matrix for PCA by Feature
feature_matrix$Functional <- functional_column
feature_matrix$Dataset <- dataset_column
summary(feature_matrix)

# Add distance column
feature_matrix <- cbind(feature_matrix, distance_column)

# Save to file
write.csv(feature_matrix,"../data/features/gene_functionality_features_latest1000all_new_din.csv", row.names = FALSE)

# reload matrix
feature_matrix <- read.csv("../data/features/gene_functionality_features_latest1000all.csv",header = TRUE,check.names = FALSE)
names(feature_matrix)[names(feature_matrix) == 'Neutral_predictor'] <- 'Random number'
names(feature_matrix)[names(feature_matrix) == 'RPKM_tissue'] <- 'Expression'
names(feature_matrix)[names(feature_matrix) == 'phyloP max_241w'] <- 'Conservation'
names(feature_matrix)[names(feature_matrix) == 'Interaction_ave'] <- 'RNA-RNA interactions'
names(feature_matrix)[names(feature_matrix) == 'gnomAD_SNP_density'] <- 'SNP density'
names(feature_matrix)[names(feature_matrix) == 'Max covariance'] <- 'RNA Covariance'


##############################
# Retrieve genetype matrices #
##############################
funcProtExon2DataPCA <- feature_matrix_numeric[feature_matrix$Dataset=="protein-coding-exon2",]
funcProtExon3DataPCA <- feature_matrix_numeric[feature_matrix$Dataset=="protein-coding-exon3",]
funcLncrnaExon1DataPCA <- feature_matrix_numeric[feature_matrix$Dataset=="lncrna-exon1",]
funcLncrnaExon2DataPCA <- feature_matrix_numeric[feature_matrix$Dataset=="lncrna-exon2",]
funcSncrnaDatasetPCA <- feature_matrix_numeric[feature_matrix$Dataset=="lncrna-exon2",]

protExon2NCDataPCA <- feature_matrix_numeric[feature_matrix$Dataset=="protein-exon2-negative-control",]
protExon3NCDataPCA <- feature_matrix_numeric[feature_matrix$Dataset=="protein-exon3-negative-control",]
lncrnaExon1NCDataPCA <- feature_matrix_numeric[feature_matrix$Dataset=="lncrna-exon1-negative-control",]
lncrnaExon2NCDataPCA <- feature_matrix_numeric[feature_matrix$Dataset=="lncrna-exon2-negative-control",]
sncrnaNCDataPCA <- feature_matrix_numeric[feature_matrix$Dataset=="short-ncrna-negative-control",]


#################
## PCA Analysis ##
##################
# Normalize data
data_normalized_all <- as.data.frame(rbind(
  as.data.frame(get_mod_zscores_equal_cols(na.omit(rbind(feature_matrix_numeric_prot,feature_matrix_numeric_lncrna,feature_matrix_numeric_sncrna)), 
                                na.omit(rbind(feature_matrix_numeric_protein_neg,feature_matrix_numeric_lncrna_neg,feature_matrix_numeric_sncrna_neg)),
                                selected_features = select_features)$zscores),
  as.data.frame(get_mod_zscores_equal_cols(na.omit(rbind(feature_matrix_numeric_protein_neg,feature_matrix_numeric_lncrna_neg,feature_matrix_numeric_sncrna_neg)), 
                                na.omit(rbind(feature_matrix_numeric_protein_neg,feature_matrix_numeric_lncrna_neg,feature_matrix_numeric_sncrna_neg)),
                                selected_features = select_features)$zscores)
)
)
data_normalized_all <- as.data.frame(sapply(data_normalized_all, convert_to_numeric)) # Matrix for PCA by Sequence
data_normalized_all <- as.data.frame(scale(data_normalized_all))

#Protein data
data_normalized_prot <- as.data.frame(rbind(
                       as.data.frame(get_mod_zscores_equal_cols(na.omit(feature_matrix_numeric_prot), na.omit(feature_matrix_numeric_protein_neg),
                         selected_features = select_features)$zscores),
                       as.data.frame(get_mod_zscores_equal_cols(na.omit(feature_matrix_numeric_protein_neg), na.omit(feature_matrix_numeric_protein_neg),
                         selected_features = select_features)$zscores)
                     )
                   )
data_normalized_prot <- as.data.frame(sapply(data_normalized_prot, convert_to_numeric)) # Matrix for PCA by Sequence
data_normalized_prot <- as.data.frame(scale(data_normalized_prot))
#lncRNA data
data_normalized_lncrna <- as.data.frame(rbind(
  as.data.frame(get_mod_zscores_equal_cols(na.omit(feature_matrix_numeric_lncrna), na.omit(feature_matrix_numeric_lncrna_neg),
                                selected_features = select_features)$zscores),
  as.data.frame(get_mod_zscores_equal_cols(na.omit(feature_matrix_numeric_lncrna_neg), na.omit(feature_matrix_numeric_lncrna_neg),
                                selected_features = select_features)$zscores)
                     )
                   )
data_normalized_lncrna <- as.data.frame(sapply(data_normalized_lncrna, convert_to_numeric)) # Matrix for PCA by Sequence
data_normalized_lncrna <- as.data.frame(scale(data_normalized_lncrna))
#short ncRNA data
data_normalized_sncrna <- as.data.frame(rbind(
  as.data.frame(get_mod_zscores_equal_cols(na.omit(feature_matrix_numeric_sncrna), na.omit(feature_matrix_numeric_sncrna_neg),
                                selected_features = select_features)$zscores),
  as.data.frame(get_mod_zscores_equal_cols(na.omit(feature_matrix_numeric_sncrna_neg), na.omit(feature_matrix_numeric_sncrna_neg),
                                selected_features = select_features)$zscores)
                     )
                   )
data_normalized_sncrna <- as.data.frame(sapply(data_normalized_sncrna, convert_to_numeric)) # Matrix for PCA by Sequence
data_normalized_sncrna <- as.data.frame(scale(data_normalized_sncrna))

data_normalized <- as.data.frame(scale(feature_matrix_numeric_no_nas))
#data_normalized <- scale(feature_matrix_numeric_no_nas)
data_normalized_prot <- scale(na.omit(feature_matrix_numeric_prot))
data_normalized_lncrna <- scale(na.omit(feature_matrix_numeric_lncrna))
data_normalized_sncrna <- scale(na.omit(feature_matrix_numeric_sncrna))

# Compute correlation
corr_matrix <- cor(data_normalized, method = "spearman")
#corr_matrix <- cor(data_normalized_prot, method = "spearman")
#corr_matrix <- cor(data_normalized_lncrna, method = "spearman")
#corr_matrix <- cor(data_normalized_sncrna, method = "spearman")


#MDS Analisys
mds.nor <- (data_normalized) %>%
  dist() %>%
  cmdscale() %>%
  as_tibble()
colnames(mds.nor) <- c("Dim.1", "Dim.2")
ggscatter(mds.nor, x = "Dim.1", y = "Dim.2", 
          size = 1,
          label = colnames(data_normalized),
          repel = TRUE)

# Heatmap plot
heatmap(corr_matrix)

#corr_obj <- rcorr(as.matrix(data_normalized_all), type = "spearman")
corr_obj <- rcorr(as.matrix(data_normalized_prot), type = "spearman")
corr_obj <- rcorr(as.matrix(data_normalized_lncrna), type = "spearman")
corr_obj <- rcorr(as.matrix(data_normalized_sncrna), type = "spearman")

corr_matrix <- corr_obj[["r"]]
pvalues <- corr_obj[["P"]]
help("rcorr")
# Number of features
n <- ncol(data_normalized_prot)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Extract p-values
#pvalues_corrected <- pvalues[lower.tri(pvalues, diag = TRUE)]

# Adjust p-values using Bonferroni
p_values_adjusted <- p.adjust(pvalues, method="bonferroni")

p_values_adjusted <- matrix(p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(p_values_adjusted) <- diag
ggcorrplot(corr_matrix, type = "lower", method = "square", lab = TRUE, 
           lab_col = "black", lab_size = 3, ggtheme = theme_void, 
           title = "Spearman correlation Heatmap - Protein coding", show.diag = TRUE, p.mat = p_values_adjusted, sig.level = 0.05, insig = "pch",
           ) +
  theme(
    plot.title = element_text(size = 24),  # Increase title size
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    
  )
ggsave("spearman_corr_sincrna_sf.png",path = "../results/latest1000all/spearman/", scale = 3, width = 1920, height = 1080, units = "px", bg = "white")

# Apply PCA
data.pca <- princomp(data_normalized_prot)
summary(data.pca)

# Inspect variable loadings
data.pca$loadings[, 1:3]

# Scree Plot
fviz_eig(data.pca, addlabels = TRUE)

# Graph of the variables
fviz_pca_var(data.pca, col.var = "black")

# Contribution of each variable
fviz_cos2(data.pca, choice = "var", axes = 1:2)

# Biplot
fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("black", "orchid", "blue"),
             repel = TRUE)


# Extract the first two principal components
pca_data <- data.frame(PC1 = data.pca$scores[,1], 
                       PC2 = data.pca$scores[,2])

# Create the contour plot
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_density_2d() + # Calculate and plot 2D density contours
  labs(title = "2D Contour Plot of PCA - Protein coding", 
       x = "Principal Component 1", y = "Principal Component 2") 

ggplot(pca_data, aes(x = PC1, y = PC2)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
  scale_fill_viridis_c() +
  geom_point(alpha = 0.3) +
  labs(title = "2D Contour Plot of PCA - Protein coding", 
       x = "PC1", y = "PC2") 
################
## Power Test ##
################
# power test for multiple comparisons, leave r null to calculate it
pwr.r.test(n = 2641, r = NULL, sig.level = 0.05, power = 0.80,
           alternative = "two.sided")



########################
# Perform PCA by Feature
# Normalize data
data_norm_t <- t(data_normalized_all)
data_t_normalized <- data_norm_t

# Compute correlation
corr_matrix_t <- cor(data_t_normalized, method = "spearman")
#ggcorrplot(corr_matrix, type = "lower", method = "square", lab = TRUE, lab_col = "black", lab_size = 3.5)

# Apply PCA
pca_result_by_feature <- prcomp(data_t_normalized, scale. = TRUE)
summary(pca_result_by_feature)
pc_scores_by_feature <- as.data.frame(pca_result_by_feature$x[, 1:2])
pca_by_feature_df <- cbind(pc_scores_by_feature)

# Add Row names as factor to the PCA by Feature result data frame
pca_by_feature_df$feature_name <- rownames(feature_matrix_numeric_t)
pca_by_feature_df$feature_name <- factor(pca_by_feature_df$feature_name, levels = c(rownames(feature_matrix_numeric_t)))



###################
# Plot PCA by Feature
ggplot(pca_by_feature_df, aes(x = PC1, y = PC2, color = feature_name, shape = feature_name)) +
  geom_point(size = 5, position = position_dodge(width = 0.1)) +
  scale_shape_manual(values = c(10, # Neutral_predictor
                                15, # GC_percent
                                18,18,18,18,18, # Dinucleotides
                                #18,18,18,18,18,18,18,18, # Dinucleotides
                                7, # Low complexity
                                16,16, # PyloP max
                                #16,16, # GERPs
                                5,5, #RPKM
                                8,8, #Copy num, DFam sum
                                14, # RNAcode 
                                6,6,6,6, # max cov, MFE, RNAalifold, Int ave
                                11,11, # gnomADs
                                13,13,13,13,13)) +  # epigenetic
  scale_color_manual(values = c("red", # Random number
                                "green", # GC_percent
                                "lightpink1","royalblue","lightcyan4","lawngreen","indianred4", #"mediumpurple4","mediumspringgreen","mediumblue", # Dinucleotides
                                #"mediumaquamarine","maroon","lightslateblue","sienna2","snow3", "orangered3","olivedrab2","red4", # Dinucleotides
                                "darkturquoise", # Low complexity
                                "dodgerblue4","forestgreen", # Xnnnw_PP_max
                                #"turquoise3","blue", # GERPs
                                "purple","orange3", # RPKM
                                "firebrick","gold3", # copy and Dfam_sum
                                "sienna2", #RNAcode,
                                "limegreen","khaki4","orangered","maroon", # max cov, MFE, RNAalifold, int ave.
                                "purple3","darkgreen", # gnomAD
                                "lightcoral","darkorchid4","tomato4","chartreuse4","cornflowerblue" # Epigenetic
                                )) +  
  labs(title = "Principal Component Analysis by Feature", x = "PC1 (15.25%)", y = "PC2 (12.02%)") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 26),  # Increase title size
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 14),    # Increase axis label size
    legend.title = element_text(size = 24),  # Increase legend title size
    legend.text = element_text(size = 20),   # Increase legend text size
    legend.key.size = unit(2, "lines"),       # Increase legend key size
  )

########################


######################
### Z-scores Utils ###
######################
# Helper Functions:
###################
help("IQR")
# Function to remove outliers from feature col in dataset df using IQR method
remove_outliers_IQR <- function(df, col) {
  df <- as.data.frame(df[[col]], nm=col)
  df <- na.omit(df)
  Q_col <- quantile(df[[col]], probs=c(.25, .75), na.rm = FALSE)
  iqr_col <- IQR(df[[col]])
  up_col <- Q_col[2] + 1.5*iqr_col # Upper Range
  low_col <- Q_col[1] - 1.5*iqr_col # Lower Range
  
  eliminated <- subset(df, df[[col]] > low_col & df[[col]] < up_col)
  return(eliminated[[col]])
}

# Function to compute z-scores for each column in `selected_features` for functional and non functional datasets
get_zscores <-function(functional_df, negative_df, selected_features) {
  list_with_zscores <- list()
  for (col in selected_features) {
    positive_col <- remove_outliers_IQR(functional_df, col)
    negative_col <- remove_outliers_IQR(negative_df, col)
    
    if(length(positive_col) == 0) {
      positive_col <- functional_df[[col]]
    }
    
    if(length(negative_col) == 0) {
      negative_col <- negative_df[[col]]
    }
    
    mean_col <- mean(negative_col, )
    sd_col <- sd(negative_col)
    
    list_with_zscores[[col]] <- (positive_col - mean_col) / sd_col
  }
  return(list_with_zscores)
}

help(mad)
help(cMedian)
get_mod_zscores <- function(functional_df, negative_df, selected_features) {
  list_with_zscores <- list()
  list_with_method <- list()
  for (col in selected_features) {
    positive_col <- remove_outliers_IQR(functional_df, col)
    negative_col <- remove_outliers_IQR(negative_df, col)
    #positive_col <- functional_df[[col]]
    #negative_col <- negative_df[[col]]
    if(length(positive_col) == 0) {
      positive_col <- functional_df[[col]]
    }
    if(length(negative_col) == 0) {
      negative_col <- negative_df[[col]]
    }
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

get_mod_zscores_equal_cols <- function(functional_df, negative_df, selected_features) {
  list_with_zscores <- list()
  list_with_method <- list()
  for (col in selected_features) {
    #positive_col <- remove_outliers_IQR(functional_df, col)
    #negative_col <- remove_outliers_IQR(negative_df, col)
    positive_col <- functional_df[[col]]
    negative_col <- negative_df[[col]]
    if(length(positive_col) == 0) {
      positive_col <- functional_df[[col]]
    }
    if(length(negative_col) == 0) {
      negative_col <- negative_df[[col]]
    }
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

# Create a separate plot for each column in the list
plot_zscores <- function(list_with_zscores, selected_features, plot_title, fill_color) {
  plot_list <- list()
  for (col in selected_features) {
    # Generate a histogram plot for the current column
    plot <- ggplot(as.data.frame(list_with_zscores[[col]], nm = col), aes_string(x = col)) +
      geom_histogram(aes(y = ..density..)) + 
      geom_density(fill = fill_color, alpha = 0.7) +
      labs(title = paste(plot_title, "- Density Distribution of", col, "Modified Z-scores"), x = paste(col, "Modified Z-score"), y = "Density Estimate") +
      theme_minimal() 
    
    # Store the plot in the list
    plot_list[[col]] <- plot
  }
  return(plot_list)
}

# Function to save z-score plots to disk
save_plots <- function(plot_list, selected_features, prefix) {
  for (col in selected_features) {
    ggsave(paste("../results/latest1000all/z-scores/",prefix,"-",col,"_.png"), 
           plot = plot(plot_list[[col]]), 
           width = 10, 
           height = 8, 
           dpi = 300, 
           bg = "white")
  }
}

# Function to generate z-score density plots from positive and negatives in one.
generate_combined_zscore_plots <- function(list_with_zscores, 
                                           list_with_zscores_neg, 
                                           selected_features, 
                                           gene_type, 
                                           colors_to_use, 
                                           prefix){
  for (col in selected_features) {
    df_with_zscores <- as.data.frame(list_with_zscores[[col]], nm = col)
    df_with_zscores$gene_type <- gene_type
    
    df_with_zscores_neg <- as.data.frame(list_with_zscores_neg[[col]], nm = col)
    df_with_zscores_neg$gene_type <- paste(gene_type,"_negative_control")
    ## Join all in one ##
    z_scores_all_df <- rbind(df_with_zscores,df_with_zscores_neg)
    z_scores_all_df$gene_type <- factor(z_scores_all_df$gene_type, 
                                             levels = unique(z_scores_all_df$gene_type))
    
    # Initialize an empty list to store plots
    plot_list_all <- list()
    
    # Generate a histogram plot for the current column
    plot <- ggplot(z_scores_all_df, aes_string(x = col, shape = "gene_type", color = "gene_type")) +
      geom_density(alpha = 0.5) +
      #geom_point(size = 3, fill = 'white') +
      scale_shape_manual(values = c(19, 17, 15, 4)) +  #c(19, 17, 15))
      scale_color_manual(values = colors_to_use) +
      #scale_fill_manual(values = c("green","orange4","brown","blue"))
      labs(title = paste("Density Distribution of", col, "Modified Z-scores"), 
           x = paste(col, "Modified Z-score"), 
           y = "Density",
           color = "Gene Type",
           shape = "Gene Type") +
      theme_minimal() +
      xlim(-5,5)
    
    # Store the plot in the list
    plot_list_all[[col]] <- plot
    # Display individual plots
    ggsave(paste("../results/latest1000all/z-scores/",prefix,col,"_.png",sep=""), plot = plot(plot_list_all[[col]]), width = 10, height = 8, dpi = 300, bg = "white")
    
  }
}

# Function to generate z-score density plots from positive and negatives in one.
generate_combined_zscore_plots_all <- function(list_with_zscores_prot, 
                                           list_with_zscores_prot_neg,
                                           list_with_zscores_sncrna, 
                                           list_with_zscores_sncrna_neg,
                                           list_with_zscores_lncrna, 
                                           list_with_zscores_lncrna_neg,
                                           selected_features, 
                                           prefix){
  for (col in selected_features) {
    df_with_zscores_prot <- as.data.frame(list_with_zscores_prot[[col]], nm = col)
    df_with_zscores_prot$gene_type <- "protein_coding"
    
    df_with_zscores_prot_neg <- as.data.frame(list_with_zscores_prot_neg[[col]], nm = col)
    df_with_zscores_prot_neg$gene_type <- "protein_coding_negative_control"
    
    df_with_zscores_sncrna <- as.data.frame(list_with_zscores_sncrna[[col]], nm = col)
    df_with_zscores_sncrna$gene_type <- "short_ncrna"
    
    df_with_zscores_sncrna_neg <- as.data.frame(list_with_zscores_sncrna_neg[[col]], nm = col)
    df_with_zscores_sncrna_neg$gene_type <- "short_ncrna_negative_control"
    
    df_with_zscores_lncrna <- as.data.frame(list_with_zscores_lncrna[[col]], nm = col)
    df_with_zscores_lncrna$gene_type <- "lncrna"
    
    df_with_zscores_lncrna_neg <- as.data.frame(list_with_zscores_lncrna_neg[[col]], nm = col)
    df_with_zscores_lncrna_neg$gene_type <- "lncrna_negative_control"
    
    ## Join all in one ##
    z_scores_all_df <- rbind(df_with_zscores_prot,df_with_zscores_prot_neg,
                             df_with_zscores_sncrna,df_with_zscores_sncrna_neg,
                             df_with_zscores_lncrna,df_with_zscores_lncrna_neg
                             )
    z_scores_all_df$gene_type <- factor(z_scores_all_df$gene_type, 
                                        levels = unique(z_scores_all_df$gene_type))
    
    # Initialize an empty list to store plots
    plot_list_all <- list()
    
    # Generate a histogram plot for the current column
    plot <- ggplot(z_scores_all_df, aes_string(x = col, shape = "gene_type", color = "gene_type")) +
      geom_density(alpha = 0.5) +
      #geom_point(size = 3, fill = 'white') +
      scale_shape_manual(values = c(19, 17, 15, 4, 4, 4)) +  #c(19, 17, 15))
      scale_color_manual(values = c("firebrick","aquamarine","limegreen","dodgerblue","darkorange","navyblue")) +
      #scale_fill_manual(values = c("green","orange4","brown","blue"))
      labs(title = paste("Density Distribution of", col, "Modified Z-scores"), 
           x = paste(col, "Modifies Z-score"), 
           y = "Density",
           color = "Gene Type",
           shape = "Gene Type") +
      theme_minimal() +
      xlim(-5,5) 
      #scale_y_log10()   # Change the y-axis to logarithmic scale
    
    # Store the plot in the list
    plot_list_all[[col]] <- plot
  
  # Display individual plots
    ggsave(paste("../results/latest1000all/z-scores/",prefix,col,"_.png", sep = ""), plot = plot(plot_list_all[[col]]), width = 10, height = 8, dpi = 300, bg = "white")
    
  }
}

make_violin_plot <- function(list_with_zscores, selected_features, plot_title) {
  plot_list <- list()
  df_long <- data.frame(row.names = c("feature","z_score"))
  for (col in selected_features) {
    df_with_zscores <- as.data.frame(list_with_zscores[[col]], nm = col)
    
    # Melt zscores for violin plot
    df_long <- rbind(df_long, df_with_zscores %>%
      pivot_longer(cols=col, names_to = "feature", values_to = "z_score"))
    
    
  }
  #violin_plot <- ggplot(df_long, aes(x = feature, y = z_score, fill = feature)) +
  #  geom_violin() +
  #  facet_wrap(~ feature, scales = "free_y") + # One plot per feature
  #  labs(title = paste(plot_title,"- Z-Scores Violin Plot"), x = "Feature", y = "Z-score") +
  #  theme_minimal() +
  #  scale_fill_manual(values = c("blue","green","pink4","orange4","brown","purple","red","darkgreen","darkorchid4",
  #                               "gold","lightcoral","lightblue4","orange3","maroon4","purple3","khaki4","limegreen","hotpink4",
  #                               "gold4","forestgreen","firebrick","dodgerblue4","deeppink4","darkviolet","darkslategrey","darksalmon","darkolivegreen4","darkcyan",
  #                               "yellowgreen","wheat4","violetred4","turquoise4","tomato4","chartreuse4","cadetblue4")) # Customize colors
  
  #ggsave(paste("../results/latest1000all/z-scores/","violinplot_",plot_title,col,"_.png", sep = ""), plot = violin_plot, width = 10, height = 8, dpi = 300, bg = "white")
  return(df_long)
}

#########################
## Protein Coding      ##
#########################
# Select features from protein negative control dataset
matrix_protein_nc <- feature_matrix_numeric[feature_matrix$Dataset=="protein-exon2-negative-control" | feature_matrix$Dataset=="protein-exon3-negative-control",select_features]
#matrix_protein_nc <- matrix_protein_nc[, select_features] # Order columns
feature_matrix_protein_neg <- matrix_protein_nc[, names(matrix_protein_nc) %in% select_features]
feature_matrix_numeric_protein_neg <- feature_matrix_protein_neg[, sapply(feature_matrix_protein_neg, is.numeric)]
summary(feature_matrix_numeric_protein_neg)

# Select features for functional data
matrix_protein <- feature_matrix_numeric[feature_matrix$Dataset=="protein-coding-exon2" | feature_matrix$Dataset=="protein-coding-exon3",select_features]
feature_matrix_prot <- matrix_protein[, names(matrix_protein) %in% select_features]
feature_matrix_numeric_prot <- feature_matrix_prot[, sapply(feature_matrix_prot, is.numeric)]
summary(feature_matrix_numeric_prot)

# Calculate z-scores
list_prot <- get_mod_zscores(feature_matrix_numeric_prot, feature_matrix_numeric_protein_neg, select_features)
list_with_zscores_prot <- list_prot$zscores
list_with_method_prot <- list_prot$method

# Store plots in a list
plot_list_prot <- plot_zscores(list_with_zscores_prot, select_features, "Protein Coding", "firebrick")

# Save individual plots to disk
save_plots(plot_list_prot, select_features, "protein")

# Violin Plot
df_long_prot <- make_violin_plot(list_with_zscores_prot, select_features, "Protein Coding")



#########################
## lncRNA              ##
#########################
# Select features from lncrna negative control dataset
matrix_lncrna_nc <- feature_matrix_numeric[feature_matrix$Dataset=="lncrna-exon1-negative-control" | feature_matrix$Dataset=="lncrna-exon2-negative-control",select_features]
#matrix_lncrna_nc <- matrix_lncrna_nc[, select_features] # Order columns
feature_matrix_lncrna_neg <- matrix_lncrna_nc[, names(matrix_lncrna_nc) %in% select_features]
feature_matrix_numeric_lncrna_neg <- feature_matrix_lncrna_neg[, sapply(feature_matrix_lncrna_neg, is.numeric)]
summary(feature_matrix_numeric_lncrna_neg)

# Select features for functional data
matrix_lncrna <- feature_matrix_numeric[feature_matrix$Dataset=="lncrna-exon1" | feature_matrix$Dataset=="lncrna-exon2",select_features]
#matrix_lncrna <- matrix_lncrna[, select_features] # Order columns
feature_matrix_lncrna <- matrix_lncrna[, names(matrix_lncrna) %in% select_features]
feature_matrix_numeric_lncrna <- feature_matrix_lncrna[, sapply(feature_matrix_lncrna, is.numeric)]
summary(feature_matrix_numeric_lncrna)
# Calculate z-scores
list_lncrna <- get_mod_zscores(feature_matrix_numeric_lncrna, feature_matrix_numeric_lncrna_neg, select_features)
list_with_zscores_lncrna <- list_lncrna$zscores
list_with_method_lncrna <- list_lncrna$method


# Generate plots
plot_list_lncrna <- plot_zscores(list_with_zscores_lncrna, select_features, "lncRNA", "darkorange")

# Save individual plots to disk
save_plots(plot_list_lncrna, select_features, "lncrna")

# Violin Plot
df_long_lncrna <- make_violin_plot(list_with_zscores_lncrna, select_features, "lncRNA")



#########################
## Short ncRNA         ##
#########################
# Select features from short ncrna negative control dataset
matrix_short_ncrna_nc <- feature_matrix_numeric[feature_matrix$Dataset=="short-ncrna-negative-control",select_features]
#matrix_short_ncrna_nc <- matrix_short_ncrna_nc[, select_features] # Order columns
feature_matrix_sncrna_neg <- matrix_short_ncrna_nc[, names(matrix_short_ncrna_nc) %in% select_features]
feature_matrix_numeric_sncrna_neg <- feature_matrix_sncrna_neg[, sapply(feature_matrix_sncrna_neg, is.numeric)]
summary(feature_matrix_numeric_sncrna_neg)
# Select features for functional data
matrix_short_ncrna <- feature_matrix_numeric[feature_matrix$Dataset=="short-ncrna",select_features]
#matrix_short_ncrna <- matrix_short_ncrna[, select_features] # Order columns
feature_matrix_sncrna <- matrix_short_ncrna[, names(matrix_short_ncrna) %in% select_features]
feature_matrix_numeric_sncrna <- feature_matrix_sncrna[, sapply(feature_matrix_sncrna, is.numeric)]
summary(feature_matrix_numeric_sncrna)
# Calculate z-scores
list_sncrna <- get_mod_zscores(feature_matrix_numeric_sncrna, feature_matrix_numeric_sncrna_neg, select_features)
list_with_zscores_sncrna <- list_sncrna$zscores
list_with_method_sncrna <- list_sncrna$method

# Initialize an empty list to store plots
plot_list_sncrna <- plot_zscores(list_with_zscores_sncrna, select_features, "Short ncRNA", "limegreen")

# Save individual plots to disk
save_plots(plot_list_sncrna, select_features, "sncrna")

# Violin Plot
df_long_sncrna <- make_violin_plot(list_with_zscores_sncrna, select_features, "Short ncRNA")


#########################
## Negative controls   ##
#########################
# Protein Negatives
# Calculate z-scores
list_protein_neg <- get_mod_zscores(feature_matrix_numeric_protein_neg, feature_matrix_numeric_protein_neg, select_features)
list_with_zscores_protein_neg <- list_protein_neg$zscores
list_with_method_protein_neg <- list_protein_neg$method

# Plot positive and negative in one
generate_combined_zscore_plots(list_with_zscores_prot, 
                               list_with_zscores_protein_neg, 
                               select_features, 
                               "protein_coding", 
                               c("firebrick","aquamarine"), 
                               "prot-neg-")

df_long_prot_neg <- make_violin_plot(list_with_zscores_protein_neg, select_features, "Protein Coding Negative Controls")

# lncRNA Negatives
# Calculate z-scores
list_lncrna_neg <- get_mod_zscores(feature_matrix_numeric_lncrna_neg, feature_matrix_numeric_lncrna_neg, select_features)
list_with_zscores_lncrna_neg <- list_lncrna_neg$zscores
list_with_method_lncrna_neg <- list_lncrna_neg$method

# Plot positive and negative in one
generate_combined_zscore_plots(list_with_zscores_lncrna, 
                               list_with_zscores_lncrna_neg, 
                               select_features, 
                               "lncrna", 
                               c("darkorange","navyblue"), 
                               "lncrna-neg-")

df_long_lncrna_neg <- make_violin_plot(list_with_zscores_lncrna_neg, select_features, "LncRNA Negative Controls")

# short ncRNA Negatives
# Calculate z-scores
list_sncrna_neg <- get_mod_zscores(feature_matrix_numeric_sncrna_neg, feature_matrix_numeric_sncrna_neg, select_features)
list_with_zscores_sncrna_neg <- list_sncrna_neg$zscores
list_with_method_sncrna_neg <- list_sncrna_neg$method


# Plot positive and negative in one
generate_combined_zscore_plots(list_with_zscores_sncrna, 
                               list_with_zscores_sncrna_neg, 
                               select_features, 
                               "short_ncrna", 
                               c("limegreen","dodgerblue"), 
                               "sncrna-neg-")

df_long_sncrna_neg <- make_violin_plot(list_with_zscores_sncrna_neg, select_features, "Short ncRNA Negative Controls")

#####################
## All gene types  ##
#####################

generate_combined_zscore_plots_all(list_with_zscores_prot, list_with_zscores_protein_neg,
                                   list_with_zscores_sncrna, list_with_zscores_sncrna_neg,
                                   list_with_zscores_lncrna, list_with_zscores_lncrna_neg,
                                   select_features, "all-in-one")


# violin plots
df_long_prot$`Gene type` <- "protein coding (+)"
df_long_sncrna$`Gene type` <- "sncRNA (+)"
df_long_lncrna$`Gene type` <- "lncRNA (+)"

df_long_prot_neg$`Gene type` <- "protein coding (-)"
df_long_sncrna_neg$`Gene type` <- "sncRNA (-)"
df_long_lncrna_neg$`Gene type` <- "lncRNA (-)"


df_long_prot$group <- "protein coding"
df_long_sncrna$group <- "sncRNA"
df_long_lncrna$group <- "lncRNA"

df_long_prot_neg$group <- "protein coding"
df_long_sncrna_neg$group <- "sncRNA"
df_long_lncrna_neg$group <- "lncRNA"


df_long_all <- rbind(df_long_prot, df_long_prot_neg,
                     df_long_sncrna, df_long_sncrna_neg,
                     df_long_lncrna, df_long_lncrna_neg)
colnames(df_long_all)
# Convert 'feature' to a factor and specify the desired order
df_long_all$feature <- factor(df_long_all$feature, 
                              levels = c("Random.number",
                                         "AA", "AC", "AG",
                                         "AT", "CA", "CC", "CpG",         
                                         "CT", "GA", "GC", "GG",
                                         "GT", "TA", "TC", "TG",
                                         "TT", "low_complexity_density", "Conservation", "phyloP.max_100w",
                                         "GERP_91_mammals_max",    "GERP_63_amniotes_max",   "Expression",   "RPKM_primary.cell",
                                         "Copy.number",  "Repeat.free",  "Fickett_score","RNAcode",     
                                         "RNA.Covariance","MFE","accessibility","RNAalifold",  
                                         "RNA.RNA.interactions" ,  "SNP.density",  "gnomAD_MAF",   "H3K27ac",     
                                         "H3K36me3",     "H3K79me2",     "chromatin_acc","methylome"), 
                              labels = c("Random number",
                                         "GC content",
                                         "AA","AC","AG",
                                         "AT","CA","CC","CpG",
                                         "CT","GA","GC","GG",
                                         "GT","TA","TC","TG",
                                         "TT",
                                         "Low complexity density","phyloP (241 mammals)","phyloP (100 vertebrates)","GERP (91 eutherian mammals)","GERP (63 amniota vertebrates)",
                                         "Expression (tissue)","Expression (primary cell)","Copy number","Repeat free","Fickett","RNAcode","RNA Covariance",
                                         "MFE","Accessibility","RNAalifold","Interaction average","SNP density","MAF",
                                         "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome"))
df_long_all$`Gene type` <- factor(df_long_all$`Gene type`, levels = c("protein coding (+)",
                                                                  "protein coding (-)",
                                                                  "sncRNA (+)",
                                                                  "sncRNA (-)",
                                                                  "lncRNA (+)",
                                                                  "lncRNA (-)"))
df_long_all$group <- factor(df_long_all$group, levels = c("protein coding","sncRNA","lncRNA"))


df_long_conservation_paper <- df_long_all[df_long_all$feature=="PhyloP max (241 mammals)" 
                                    | df_long_all$feature=="PhyloP max (100 vertebrates)" 
                                    | df_long_all$feature=="GERP max (91 eutherian mammals)" 
                                    | df_long_all$feature=="GERP max (63 amniota vertebrates)",]
#df_long_conservation <- df_long_all[df_long_all$feature=="phyloP max_241w", ]
#summary(df_long_all[df_long_all$feature=="phyloP max_241w",])
df_long_expression_paper <- df_long_all[df_long_all$feature=="Tissue RPKM" 
                                  | df_long_all$feature=="Primary cell RPKM", ]
df_long_repeat_assoc_coding <- df_long_all[df_long_all$feature=="Repeat free region" 
                                           | df_long_all$feature=="Genomic copy number", ]
summary(df_long_all$feature)
custom_order <- c("Repeat free region", "Genomic copy number")
df_long_repeat_assoc_coding$feature <- factor(df_long_repeat_assoc_coding$feature, levels = custom_order)
df_long_repeat_assoc_coding_ordered <- df_long_repeat_assoc_coding %>%
  arrange(feature)


df_long_con_exp <- df_long_all[df_long_all$feature=="Tissue RPKM" 
                                   | df_long_all$feature=="Primary cell RPKM", ]


df_long_con_exp_epi_poster <- df_long_all[df_long_all$feature=="Expression (tissue)" 
                               | df_long_all$feature=="Conservation (241 mammals)"
                               | df_long_all$feature=="Histone mark (H3K36me3)" , ]

df_long_struct_poster <- df_long_all[df_long_all$feature=="RNA Covariance" 
                              | df_long_all$feature=="MFE" 
                              | df_long_all$feature=="RNAalifold" 
                              | df_long_all$feature=="Accessibility" 
                              | df_long_all$feature=="RNA-RNA Interactions"
                              | df_long_all$feature=="RNAcode"
                              | df_long_all$feature=="Fickett", ]

df_long_struct_paper <- df_long_all[df_long_all$feature=="RNA Covariance" 
                                     | df_long_all$feature=="MFE" 
                                     | df_long_all$feature=="RNAalifold" 
                                     | df_long_all$feature=="Accessibility" 
                                     | df_long_all$feature=="Interaction average"
                                     | df_long_all$feature=="RNAcode"
                                     | df_long_all$feature=="Fickett", ]


df_long_struct_pop_poster <- df_long_all[df_long_all$feature=="RNA Covariance" 
                              | df_long_all$feature=="RNA-RNA Interactions"
                              | df_long_all$feature=="SNP density", ]

df_long_pop_var <- df_long_all[df_long_all$feature=="MAF"
                                  | df_long_all$feature=="SNP density", ]

df_long_coding_potential <- df_long_all[df_long_all$feature=="RNAcode coding potential"
                              | df_long_all$feature=="Fickett score", ]
df_long_intrinsic1 <- df_long_all[df_long_all$feature=="GC%" 
                                  | df_long_all$feature=="Low complexity density"
                                  | df_long_all$feature=="CpG dinucleotide content"
                                  | df_long_all$feature=="GG dinucleotide content"
                                  | df_long_all$feature=="TA dinucleotide content"
                                  | df_long_all$feature=="GA dinucleotide content"
                                  | df_long_all$feature=="GT dinucleotide content"
                                  | df_long_all$feature=="AC dinucleotide content"
                                  | df_long_all$feature=="CC dinucleotide content", ]

custom_order <- c("GC%", "Low complexity density","GA dinucleotide content", "CpG dinucleotide content", 
                  "TA dinucleotide content","GG dinucleotide content","GT dinucleotide content",
                  "AC dinucleotide content","CC dinucleotide content")
df_long_intrinsic1$feature <- factor(df_long_intrinsic1$feature, levels = custom_order)
df_long_intrinsic1_ordered <- df_long_intrinsic1 %>%
  arrange(feature)

df_long_intrinsic2 <- df_long_all[df_long_all$feature=="GC"
                                  | df_long_all$feature=="TC"
                                  | df_long_all$feature=="TG"
                                  | df_long_all$feature=="TT"
                                  | df_long_all$feature=="AA" 
                                  | df_long_all$feature=="AC" 
                                  | df_long_all$feature=="AG" 
                                  | df_long_all$feature=="AT" 
                                  |  df_long_all$feature=="CA"
                                  | df_long_all$feature=="CC"
                                  | df_long_all$feature=="CT", ]

df_long_epigen_1_paper <- df_long_all[df_long_all$feature=="H3K27ac" 
                                     | df_long_all$feature=="H3K36me3" 
                                     | df_long_all$feature=="Chromatin accessibility"
                                     | df_long_all$feature=="H3K79me2"
                                     | df_long_all$feature=="Methylome", ]

df_long_epigen_2_paper <- df_long_all[df_long_all$feature=="H3K79me2"
                                      | df_long_all$feature=="Methylome", ]


df_long_random_number <- df_long_all[df_long_all$feature=="Random number",]

df_long_cons_exp_epi <- df_long_all[df_long_all$feature=="phyloP max_241w" 
                                    | df_long_all$feature=="phyloP max_100w" 
                                    | df_long_all$feature=="GERP_91_mammals_max" 
                                    | df_long_all$feature=="RPKM_tissue" 
                                    | df_long_all$feature=="RPKM_primary cell"
                                    | df_long_all$feature=="H3K36me3" 
                                    | df_long_all$feature=="H3K79me2" ,]


ggplot(df_long_all, aes(x = feature, y = z_score, fill = gene_type)) +
  geom_violin() +
  geom_boxplot(alpha=0.3, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "Z-Scores Violin Plot", x = "Feature", y = "Z-score") +
  theme_minimal() +
  scale_fill_manual(values = c("firebrick","aquamarine","limegreen","dodgerblue","darkorange","navyblue")) + # Customize colors
  ylim(-4,4)

ggplot(df_long_cons_exp_epi, aes(x = feature, y = z_score, fill = gene_type)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 26) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "right",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend text size
    plot.title = element_text(size = 36, hjust = 0.5)  # Increase plot title size and center it
  ) +
  scale_fill_manual(values = c("firebrick","aquamarine","limegreen","dodgerblue","darkorange","navyblue")) + # Customize colors
  ylim(-3,6)

ggplot(df_long_con_exp, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 29) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "right",
    legend.title = element_text(size = 24),  # Increase legend title size
    legend.text = element_text(size = 22),  # Increase legend text size
    plot.title = element_text(size = 40, hjust = 0.5)  # Increase plot title size and center it
  ) +
  scale_fill_manual(values = c("firebrick","aquamarine","limegreen","dodgerblue","darkorange","navyblue")) #+ # Customize colors
  #ylim(-2,10)

ggsave("cons_exp_epi.png",path = "../results/latest1000all/violinPlots/Poster/", scale = 3, width = 3840, height = 1080, units = "px", bg = "white", dpi = 600)


ggplot(df_long_pop_var, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 44) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "right",
    legend.title = element_text(size = 24),  # Increase legend title size
    legend.text = element_text(size = 22),  # Increase legend text size
    plot.title = element_text(size = 40, hjust = 0.5)  # Increase plot title size and center it
  ) +
  scale_fill_manual(values = c("firebrick","aquamarine","limegreen","dodgerblue","darkorange","navyblue")) + # Customize colors
  ylim(-3,6)

ggsave("pop_var.png",path = "../results/latest1000all/violinPlots/Poster/", scale = 3, width = 3840, height = 1080, units = "px", bg = "white", dpi = 600)


help("element_text")
ggplot(df_long_conservation_paper, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
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
    strip.text = element_text(margin = margin(0,0,40,0))
  ) +
  scale_fill_manual(values = c("#F4A582FF","#c9e3f6FF","#D6604DFF","#56bdfcFF","#e37b88FF","#53a4f5FF")) + # Customize colors
  coord_cartesian(ylim = c(-3, 6))
  #ylim(-3,6)

ggsave("conservation.png",path = "../results/latest1000all/violinPlots/Paper/newColors/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


ggplot(df_long_expression_paper, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "Expression", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 34) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "none",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend text size
    plot.title = element_text(size = 0, hjust = 0.5),  # Increase plot title size and center it
    axis.line.x = element_blank(),          # Remove x-axis line
    axis.ticks.x = element_blank(),         # Remove x-axis ticks
    panel.grid.major.x = element_blank(),   # Remove major grid lines along x-axis
    panel.grid.minor.x = element_blank(),    # Remove minor grid lines along x-axis
    strip.text = element_text(margin = margin(0,0,40,0))
  ) +
  scale_fill_manual(values = c("#F4A582FF","#c9e3f6FF","#D6604DFF","#56bdfcFF","#e37b88FF","#53a4f5FF")) + # Customize colors
  coord_cartesian(ylim = c(-10, 10000))
  #ylim(-2,10)
ggsave("expressionBig.png",path = "../results/latest1000all/violinPlots/Paper/newColors/", scale = 3, width = 3840, height = 1080, units = "px", bg = "white", dpi = 600)


ggplot(df_long_struct_pop, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 34) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "right",
    legend.title = element_text(size = 24),  # Increase legend title size
    legend.text = element_text(size = 22),  # Increase legend text size
    plot.title = element_text(size = 36, hjust = 0.5)  # Increase plot title size and center it
  ) +
  scale_fill_manual(values = c("firebrick","aquamarine","limegreen","dodgerblue","darkorange","navyblue"))  +# Customize colors
  ylim(-4,5)
ggsave("struct_pop.png",path = "../results/latest1000all/violinPlots/Poster/", scale = 3, width = 3840, height = 1080, units = "px", bg = "white", dpi = 600)


ggplot(df_long_repeat_assoc_coding_ordered, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin() +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "Repeat Association", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 34) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "right",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend text size
    plot.title = element_text(size = 36, hjust = 0.5)  # Increase plot title size and center it
  ) +
  scale_fill_manual(values = c("#F4A582FF","#c9e3f6FF","#D6604DFF","#56bdfcFF","#e37b88FF","#53a4f5FF")) + # Customize colors
  #coord_cartesian(ylim = c(-10, 500))
  ylim(-1,5)
ggsave("repeat_assoc_coding.png",path = "../results/latest1000all/violinPlots/Paper/newColors", scale = 3, width = 3840, height = 1080, units = "px", bg = "white", dpi = 600)


ggplot(df_long_struct_paper, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "Protein & ncRNA specific: structure, interactions & coding potential", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 34) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "right",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend text size
    plot.title = element_text(size = 36, hjust = 0.5)  # Increase plot title size and center it
  ) +
  scale_fill_manual(values = c("firebrick","aquamarine","limegreen","dodgerblue","darkorange","navyblue")) + # Customize colors
ylim(-4,6)
#help("geom_violin")
ggsave("specific.png",path = "../results/latest1000all/violinPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


ggplot(df_long_coding_potential, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 34) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "none",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend text size
    plot.title = element_text(size = 0, hjust = 0.5),  # Increase plot title size and center it
    axis.line.x = element_blank(),          # Remove x-axis line
    axis.ticks.x = element_blank(),         # Remove x-axis ticks
    panel.grid.major.x = element_blank(),   # Remove major grid lines along x-axis
    panel.grid.minor.x = element_blank(),    # Remove minor grid lines along x-axis
    strip.text = element_text(margin = margin(0,0,40,0))
  ) +
  scale_fill_manual(values = c("#F4A582FF","#c9e3f6FF","#D6604DFF","#56bdfcFF","#e37b88FF","#53a4f5FF")) #+ # Customize colors
  #coord_cartesian(ylim = c(-10, 500))
  #ylim(-4,3)
help("geom_violin")

ggplot(df_long_intrinsic1, aes(x = feature, y = z_score, fill = gene_type)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 28) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "right",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 14),  # Increase legend text size
    plot.title = element_text(size = 36, hjust = 0.5)  # Increase plot title size and center it
  ) +
  scale_fill_manual(values = c("#F4A582FF","#c9e3f6FF","#D6604DFF","#56bdfcFF","#e37b88FF","#53a4f5FF")) + # Customize colors
  ylim(-4,3)

ggplot(df_long_intrinsic1_ordered, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "Intrinsic", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 34) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "none",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend text size
    plot.title = element_text(size = 0, hjust = 0.5),  # Increase plot title size and center it
    #plot.subtitle = element_blank(),
    axis.line.x = element_blank(),          # Remove x-axis line
    axis.ticks.x = element_blank(),         # Remove x-axis ticks
    panel.grid.major.x = element_blank(),   # Remove major grid lines along x-axis
    panel.grid.minor.x = element_blank(),    # Remove minor grid lines along x-axis
    strip.text = element_text(margin = margin(0,0,40,0))
  ) +
  scale_fill_manual(values = c("#F4A582FF","#c9e3f6FF","#D6604DFF","#56bdfcFF","#e37b88FF","#53a4f5FF")) + # Customize colors
  coord_cartesian(ylim = c(-4, 3))
  #ylim(-4,3)
ggsave("intrinsic.png",path = "../results/latest1000all/violinPlots/Paper/newColors/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

help("geom_violin")
ggplot(df_long_epigen_1_paper, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "Epigenetic Signatures", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 34) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "none",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend text size
    plot.title = element_text(size = 0, hjust = 0.5),  # Increase plot title size and center it
    plot.subtitle = element_blank(),
    axis.line.x = element_blank(),          # Remove x-axis line
    axis.ticks.x = element_blank(),         # Remove x-axis ticks
    panel.grid.major.x = element_blank(),   # Remove major grid lines along x-axis
    panel.grid.minor.x = element_blank(),    # Remove minor grid lines along x-axis
    strip.text = element_text(margin = margin(0,0,40,0))
  ) +
  scale_fill_manual(values = c("#F4A582FF","#c9e3f6FF","#D6604DFF","#56bdfcFF","#e37b88FF","#53a4f5FF")) + # Customize colors
  coord_cartesian(ylim = c(-3, 600))
  #ylim(-3,6)
ggsave("epigeneticNoLims.png",path = "../results/latest1000all/violinPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


ggplot(df_long_epigen_2_paper, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "Epigenetic Signatures", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 34) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "none",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend text size
    plot.title = element_text(size = 0, hjust = 0.5),  # Increase plot title size and center it
    plot.subtitle = element_blank(),
    axis.line.x = element_blank(),          # Remove x-axis line
    axis.ticks.x = element_blank(),         # Remove x-axis ticks
    panel.grid.major.x = element_blank(),   # Remove major grid lines along x-axis
    panel.grid.minor.x = element_blank(),    # Remove minor grid lines along x-axis
    strip.text = element_text(margin = margin(0,0,40,0))
  ) +
  scale_fill_manual(values = c("#F4A582FF","#c9e3f6FF","#D6604DFF","#56bdfcFF","#e37b88FF","#53a4f5FF")) + # Customize colors
  coord_cartesian(ylim = c(-3, 6))
#ylim(-3,6)
ggsave("epigeneticMeth.png",path = "../results/latest1000all/violinPlots/Paper/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)



ggplot(df_long_random_number, aes(x = feature, y = z_score, fill = `Gene type`)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=TRUE, position = position_dodge(width = 0.9), width=0.2) +
  facet_wrap(~ feature, scales = "free") + 
  labs(title = "", x = "Feature", y = "Robust Z-score") +
  theme_minimal(base_size = 34) +
  theme(
    axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "right",
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12),  # Increase legend text size
    plot.title = element_text(size = 36, hjust = 0.5)  # Increase plot title size and center it
  ) +
  scale_fill_manual(values = c("firebrick","aquamarine","limegreen","dodgerblue","darkorange","navyblue"))  # Customize colors
  #ylim(-1.5,1.5)
ggsave("random.png",path = "../results/latest1000all/violinPlots/Paper/", scale = 3, width = 1920, height = 1080, units = "px", bg = "white", dpi = 600)


################################################################################
################################################################################
################################################################################

# Code used to correct some of the plots, change col to the desired feature and
# modify xlim and ylim to get proper scaling.

col="RNAalifold"
prefix="all-in-one"

df_with_zscores_prot <- as.data.frame(list_with_zscores_prot[[col]], nm = col)
df_with_zscores_prot$gene_type <- "protein_coding"

df_with_zscores_prot_neg <- as.data.frame(list_with_zscores_protein_neg[[col]], nm = col)
df_with_zscores_prot_neg$gene_type <- "protein_coding_negative_control"


df_with_zscores_lncrna <- as.data.frame(list_with_zscores_lncrna[[col]], nm = col)
df_with_zscores_lncrna$gene_type <- "lncrna"

df_with_zscores_lncrna_neg <- as.data.frame(list_with_zscores_lncrna_neg[[col]], nm = col)
df_with_zscores_lncrna_neg$gene_type <- "lncrna_negative_control"


df_with_zscores_sncrna <- as.data.frame(list_with_zscores_sncrna[[col]], nm = col)
df_with_zscores_sncrna$gene_type <- "short_ncrna"

df_with_zscores_sncrna_neg <- as.data.frame(list_with_zscores_sncrna_neg[[col]], nm = col)
df_with_zscores_sncrna_neg$gene_type <- "short_ncrna_negative_control"

## Join all in one ##
z_scores_all_df <- rbind(df_with_zscores_prot,df_with_zscores_prot_neg,
                         df_with_zscores_lncrna,df_with_zscores_lncrna_neg,
                         df_with_zscores_sncrna,df_with_zscores_sncrna_neg)
z_scores_all_df$gene_type <- factor(z_scores_all_df$gene_type, 
                                    levels = unique(z_scores_all_df$gene_type))
summary(z_scores_all_df)

#For log scale, use value shift to avoid negatives and 0s
shifted_z_scores_all_df <- z_scores_all_df
shifted_z_scores_all_df[col] <- z_scores_all_df[col] + 3.827
summary(shifted_z_scores_all_df)

# Initialize an empty list to store plots
plot_list_all <- list()

# Generate a histogram plot for the current column
plot <- ggplot(shifted_z_scores_all_df, aes_string(x = col, shape = "gene_type", color = "gene_type")) +
  geom_density(alpha = 0.5) +
  #geom_point(size = 3, fill = 'white') +
  scale_shape_manual(values = c(19, 17, 15, 4, 4, 4)) +  #c(19, 17, 15))
  scale_color_manual(values = c("firebrick","aquamarine","darkorange","navyblue","limegreen","dodgerblue")) +
  #scale_fill_manual(values = c("green","orange4","brown","blue"))
  labs(title = paste("Density Distribution of", col, "Z-scores"), 
       x = paste(col, "Z-score"), 
       y = "Density",
       color = "Gene Type",
       shape = "Gene Type") +
  theme_minimal() +
  #xlim(0,1) +
  ylim(0,20) +
  scale_x_log10() 
  #scale_y_log10()   # Change the y-axis to logarithmic scale

# Store the plot in the list
plot_list_all[[col]] <- plot

# Display individual plots
ggsave(paste("../results/z-scores/",prefix,col,"_.png", sep = ""), plot = plot(plot_list_all[[col]]), width = 10, height = 8, dpi = 300, bg = "white")



########################
########################

col='RNAalifold'
gene_type='protein_coding' 
colors_to_use=c("firebrick","aquamarine")
prefix='prot-neg-'
list_with_zscores=list_with_zscores_prot
list_with_zscores_neg=list_with_zscores_protein_neg

df_with_zscores <- as.data.frame(list_with_zscores[[col]], nm = col)
df_with_zscores$gene_type <- gene_type

df_with_zscores_neg <- as.data.frame(list_with_zscores_neg[[col]], nm = col)
df_with_zscores_neg$gene_type <- paste(gene_type,"_negative_control")
## Join all in one ##
z_scores_all_df <- rbind(df_with_zscores,df_with_zscores_neg)
z_scores_all_df$gene_type <- factor(z_scores_all_df$gene_type, 
                                    levels = unique(z_scores_all_df$gene_type))

# Initialize an empty list to store plots
plot_list_all <- list()

# Generate a histogram plot for the current column
plot <- ggplot(z_scores_all_df, aes_string(x = col, shape = "gene_type", color = "gene_type")) +
  geom_density(alpha = 0.5) +
  #geom_point(size = 3, fill = 'white') +
  scale_shape_manual(values = c(19, 17, 15, 4)) +  #c(19, 17, 15))
  scale_color_manual(values = colors_to_use) +
  #scale_fill_manual(values = c("green","orange4","brown","blue"))
  labs(title = paste("Density Distribution of", col, "Z-scores"), 
       x = paste(col, "Z-score"), 
       y = "Density",
       color = "Gene Type",
       shape = "Gene Type") +
  theme_minimal() +
  xlim(-1,1) +
  ylim(0,10)

# Store the plot in the list
plot_list_all[[col]] <- plot
# Display individual plots
ggsave(paste("../results/z-scores/",prefix,col,"_.png",sep=""), plot = plot(plot_list_all[[col]]), width = 10, height = 8, dpi = 300, bg = "white")
