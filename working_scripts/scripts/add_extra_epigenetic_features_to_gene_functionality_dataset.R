# Script to add extra epigenetic features of interest to the main set of gene fuctionality features.
# Load necessary libraries.
library(ggplot2)
library(tidyr)

source("scripts/config.R")

source("scripts/load_epigenetic_features.R")
source("scripts/utils.R")

# Load existing gene functionality features
funcProtExon2Data <- data.frame(read.csv(FUNC_PROT_EXON2_FEATURES_FILE, header=TRUE))
funcProtExon3Data <- data.frame(read.csv(FUNC_PROT_EXON3_FEATURES_FILE, header=TRUE))
funcLncrnaExon1Data <- data.frame(read.csv(FUNC_LNCRNA_EXON1_FEATURES_FILE, header=TRUE))
funcLncrnaExon2Data <- data.frame(read.csv(FUNC_LNCRNA_EXON2_FEATURES_FILE, header=TRUE))
funcSncrnaData <- data.frame(read.csv(FUNC_SNCRNA_FEATURES_FILE, header=TRUE))

protExon2NCData <- data.frame(read.csv(NC_PROT_EXON2_FEATURES_FILE, header=TRUE))
protExon3NCData <- data.frame(read.csv(NC_PROT_EXON3_FEATURES_FILE, header=TRUE))
lncrnaExon1NCData <- data.frame(read.csv(NC_LNCRNA_EXON1_FEATURES_FILE, header=TRUE))
lncrnaExon2NCData <- data.frame(read.csv(NC_LNCRNA_EXON2_FEATURES_FILE, header=TRUE))
sncrnaNCData <- data.frame(read.csv(NC_SNCRNA_FEATURES_FILE, header=TRUE))

# Load epigenetic features
epigenetic_matrix <- load_epigenetic_features()

# Add extra features
#H3K79me1:
funcProtExon2Data$H3K79me1_MaxScaledSignal <- data.frame(epigenetic_matrix %>%
  filter(Dataset=="protein-coding-exon2"))$H3K79me1_MaxScaledSignal
funcProtExon3Data$H3K79me1_MaxScaledSignal <- data.frame(epigenetic_matrix %>%
  filter(Dataset=="protein-coding-exon3"))$H3K79me1_MaxScaledSignal

funcLncrnaExon1Data$H3K79me1_MaxScaledSignal <- data.frame(epigenetic_matrix %>%
  filter(Dataset=="lncrna-exon1"))$H3K79me1_MaxScaledSignal
funcLncrnaExon2Data$H3K79me1_MaxScaledSignal <- data.frame(epigenetic_matrix %>%
  filter(Dataset=="lncrna-exon2"))$H3K79me1_MaxScaledSignal

funcSncrnaData$H3K79me1_MaxScaledSignal <- data.frame(epigenetic_matrix %>%
  filter(Dataset=="short-ncrna"))$H3K79me1_MaxScaledSignal

# Negative cases
protExon2NCData$H3K79me1_MaxScaledSignal <- data.frame(epigenetic_matrix %>%
  filter(Dataset=="protein-exon2-negative-control"))$H3K79me1_MaxScaledSignal
protExon3NCData$H3K79me1_MaxScaledSignal <- data.frame(epigenetic_matrix %>%
  filter(Dataset=="protein-exon3-negative-control"))$H3K79me1_MaxScaledSignal

lncrnaExon1NCData$H3K79me1_MaxScaledSignal <- data.frame(epigenetic_matrix %>%
  filter(Dataset=="lncrna-exon1-negative-control"))$H3K79me1_MaxScaledSignal
lncrnaExon2NCData$H3K79me1_MaxScaledSignal <- data.frame(epigenetic_matrix %>%
  filter(Dataset=="lncrna-exon2-negative-control"))$H3K79me1_MaxScaledSignal

sncrnaNCData$H3K79me1_MaxScaledSignal <- data.frame(epigenetic_matrix %>%
  filter(Dataset=="short-ncrna-negative-control"))$H3K79me1_MaxScaledSignal


# Write files
write.csv(funcProtExon2Data,"data/features/functional-protein-exon2-dataset-features.csv", row.names = FALSE)
write.csv(funcProtExon3Data,"data/features/functional-protein-exon3-dataset-features.csv", row.names = FALSE)
write.csv(funcLncrnaExon1Data,"data/features/functional-lncrna-exon1-dataset-features.csv", row.names = FALSE)
write.csv(funcLncrnaExon2Data,"data/features/functional-lncrna-exon2-dataset-features.csv", row.names = FALSE)
write.csv(funcSncrnaData,"data/features/functional-short-ncrna-dataset-features.csv", row.names = FALSE)

write.csv(protExon2NCData,"data/features/protein-exon2-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(protExon3NCData,"data/features/protein-exon3-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(lncrnaExon1NCData,"data/features/lncrna-exon1-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(lncrnaExon2NCData,"data/features/lncrna-exon2-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(sncrnaNCData,"data/features/short-ncrna-negative-control-dataset-features.csv", row.names = FALSE)
