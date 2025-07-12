#Load existing features
#source("load_gene_functionality_features.R")
funcProtExon2Data <- data.frame(read.csv("../results/functional-protein-exon2-dataset-features.csv", header=TRUE))
funcProtExon3Data <- data.frame(read.csv("../results/functional-protein-exon3-dataset-features.csv", header=TRUE))
funcLncrnaExon1Data <- data.frame(read.csv("../results/functional-lncrna-exon1-dataset-features.csv", header=TRUE))
funcLncrnaExon2Data <- data.frame(read.csv("../results/functional-lncrna-exon2-dataset-features.csv", header=TRUE))
funcSncrnaDataset <- data.frame(read.csv("../results/functional-short-ncrna-dataset-features.csv", header=TRUE))

protExon2NCData <- data.frame(read.csv("../results/protein-exon2-negative-control-dataset-features.csv", header=TRUE))
protExon3NCData <- data.frame(read.csv("../results/protein-exon3-negative-control-dataset-features.csv", header=TRUE))
lncrnaExon1NCData <- data.frame(read.csv("../results/lncrna-exon1-negative-control-dataset-features.csv", header=TRUE))
lncrnaExon2NCData <- data.frame(read.csv("../results/lncrna-exon2-negative-control-dataset-features.csv", header=TRUE))
sncrnaNCData <- data.frame(read.csv("../results/short-ncrna-negative-control-dataset-features.csv", header=TRUE))


# Load new methylome
functional_protein_exon2_methylome_feature <- read.csv("../data/datasets/methylome_feature/protein-exon2-methylome-feature.csv", header = TRUE)
functional_protein_exon3_methylome_feature <- read.csv("../data/datasets/methylome_feature/protein-exon3-methylome-feature.csv", header = TRUE)
functional_lncrna_exon1_methylome_feature <- read.csv("../data/datasets/methylome_feature/lncrna-exon1-methylome-feature.csv", header = TRUE)
functional_lncrna_exon2_methylome_feature <- read.csv("../data/datasets/methylome_feature/lncrna-exon2-methylome-feature.csv", header = TRUE)
functional_sncrna_methylome_feature <- read.csv("../data/datasets/methylome_feature/short-ncrna-methylome-feature.csv", header = TRUE)

methylome_column_lncrna <- rbind(as.data.frame(functional_lncrna_exon1_methylome_feature), as.data.frame(functional_lncrna_exon2_methylome_feature))

protein_exon2_nc_methylome_feature <- read.csv("../data/datasets/methylome_feature/protein-exon2-NC-methylome-feature.csv", header = TRUE)
protein_exon3_nc_methylome_feature <- read.csv("../data/datasets/methylome_feature/protein-exon3-NC-methylome-feature.csv", header = TRUE)
lncrna_exon1_nc_methylome_feature <- read.csv("../data/datasets/methylome_feature/lncrna-exon1-NC-methylome-feature.csv", header = TRUE)
lncrna_exon2_nc_methylome_feature <- read.csv("../data/datasets/methylome_feature/lncrna-exon2-NC-methylome-feature.csv", header = TRUE)
sncrna_exon2_nc_methylome_feature <- read.csv("../data/datasets/methylome_feature/short-ncrna-NC-methylome-feature.csv", header = TRUE)

methylome_column_lncrna_nc <- rbind(as.data.frame(lncrna_exon1_nc_methylome_feature),as.data.frame(lncrna_exon2_nc_methylome_feature))

# Replace methylome
funcProtExon2Data$methylome <- functional_protein_exon2_methylome_feature$methylome
funcProtExon3Data$methylome <- functional_protein_exon3_methylome_feature$methylome
funcLncrnaExon1Data$methylome <- functional_lncrna_exon1_methylome_feature$methylome
funcLncrnaExon2Data$methylome <- functional_lncrna_exon2_methylome_feature$methylome
funcSncrnaDataset$methylome <- functional_sncrna_methylome_feature$methylome

protExon2NCData$methylome <- protein_exon2_nc_methylome_feature$methylome
protExon3NCData$methylome <- protein_exon3_nc_methylome_feature$methylome
lncrnaExon1NCData$methylome <- lncrna_exon1_nc_methylome_feature$methylome
lncrnaExon2NCData$methylome <- lncrna_exon2_nc_methylome_feature$methylome
sncrnaNCData$methylome <- sncrna_exon2_nc_methylome_feature$methylome



# Load new accessibility
functional_protein_exon2_accessibility_feature <- read.csv("../data/datasets/accessibility_feature/functional-protein-exon2-accessibility-feature.csv", header = TRUE)
functional_protein_exon3_accessibility_feature <- read.csv("../data/datasets/accessibility_feature/functional-protein-exon3-accessibility-feature.csv", header = TRUE)
functional_lncrna_exon1_accessibility_feature <- read.csv("../data/datasets/accessibility_feature/functional-lncrna-exon1-accessibility-feature.csv", header = TRUE)
functional_lncrna_exon2_accessibility_feature <- read.csv("../data/datasets/accessibility_feature/functional-lncrna-exon2-accessibility-feature.csv", header = TRUE)
functional_sncrna_accessibility_feature <- read.csv("../data/datasets/accessibility_feature/functional-short-ncrna-accessibility-feature.csv", header = TRUE)

accessibility_column_lncrna <- rbind(as.data.frame(functional_lncrna_exon1_accessibility_feature), as.data.frame(functional_lncrna_exon2_accessibility_feature))

protein_exon2_nc_accessibility_feature <- read.csv("../data/datasets/accessibility_feature/protein-exon2-negative-control-accessibility-feature.csv", header = TRUE)
protein_exon3_nc_accessibility_feature <- read.csv("../data/datasets/accessibility_feature/protein-exon3-negative-control-accessibility-feature.csv", header = TRUE)
lncrna_exon1_nc_accessibility_feature <- read.csv("../data/datasets/accessibility_feature/lncrna-exon1-negative-control-accessibility-feature.csv", header = TRUE)
lncrna_exon2_nc_accessibility_feature <- read.csv("../data/datasets/accessibility_feature/lncrna-exon2-negative-control-accessibility-feature.csv", header = TRUE)
sncrna_exon2_nc_accessibility_feature <- read.csv("../data/datasets/accessibility_feature/short-ncrna-negative-control-accessibility-feature.csv", header = TRUE)

accessibility_column_lncrna_nc <- rbind(as.data.frame(lncrna_exon1_nc_accessibility_feature),as.data.frame(lncrna_exon2_nc_accessibility_feature))

# Replace Accessibility
funcProtExon2Data$accessibility <- functional_protein_exon2_accessibility_feature$accessibility
funcProtExon3Data$accessibility <- functional_protein_exon3_accessibility_feature$accessibility
funcLncrnaExon1Data$accessibility <- functional_lncrna_exon1_accessibility_feature$accessibility
funcLncrnaExon2Data$accessibility <- functional_lncrna_exon2_accessibility_feature$accessibility
funcSncrnaDataset$accessibility <- functional_sncrna_accessibility_feature$accessibility

protExon2NCData$accessibility <- protein_exon2_nc_accessibility_feature$accessibility
protExon3NCData$accessibility <- protein_exon3_nc_accessibility_feature$accessibility
lncrnaExon1NCData$accessibility <- lncrna_exon1_nc_accessibility_feature$accessibility
lncrnaExon2NCData$accessibility <- lncrna_exon2_nc_accessibility_feature$accessibility
sncrnaNCData$accessibility <- sncrna_exon2_nc_accessibility_feature$accessibility

# Load fickett
functional_protein_exon2_fickett_feature <- read.csv("../data/datasets/fickett_feature/functional-protein-exon2-fickett-feature.csv", header = TRUE)
functional_protein_exon3_fickett_feature <- read.csv("../data/datasets/fickett_feature/functional-protein-exon3-fickett-feature.csv", header = TRUE)
functional_lncrna_exon1_fickett_feature <- read.csv("../data/datasets/fickett_feature/functional-lncrna-exon1-fickett-feature.csv", header = TRUE)
functional_lncrna_exon2_fickett_feature <- read.csv("../data/datasets/fickett_feature/functional-lncrna-exon2-fickett-feature.csv", header = TRUE)
functional_sncrna_fickett_feature <- read.csv("../data/datasets/fickett_feature/functional-short-ncrna-fickett-feature.csv", header = TRUE)

protein_exon2_nc_fickett_feature <- read.csv("../data/datasets/fickett_feature/protein-exon2-negative-control-fickett-feature.csv", header = TRUE)
protein_exon3_nc_fickett_feature <- read.csv("../data/datasets/fickett_feature/protein-exon3-negative-control-fickett-feature.csv", header = TRUE)
lncrna_exon1_nc_fickett_feature <- read.csv("../data/datasets/fickett_feature/lncrna-exon1-negative-control-fickett-feature.csv", header = TRUE)
lncrna_exon2_nc_fickett_feature <- read.csv("../data/datasets/fickett_feature/lncrna-exon2-negative-control-fickett-feature.csv", header = TRUE)
sncrna_exon2_nc_fickett_feature <- read.csv("../data/datasets/fickett_feature/short-ncrna-negative-control-fickett-feature.csv", header = TRUE)

# Replace fickett
funcProtExon2Data$fickett <- functional_protein_exon2_fickett_feature$Fickett_score
funcProtExon3Data$fickett <- functional_protein_exon3_fickett_feature$Fickett_score
funcLncrnaExon1Data$fickett <- functional_lncrna_exon1_fickett_feature$Fickett_score
funcLncrnaExon2Data$fickett <- functional_lncrna_exon2_fickett_feature$Fickett_score
funcSncrnaDataset$fickett <- functional_sncrna_fickett_feature$Fickett_score

protExon2NCData$fickett <- protein_exon2_nc_fickett_feature$Fickett_score
protExon3NCData$fickett <- protein_exon3_nc_fickett_feature$Fickett_score
lncrnaExon1NCData$fickett <- lncrna_exon1_nc_fickett_feature$Fickett_score
lncrnaExon2NCData$fickett <- lncrna_exon2_nc_fickett_feature$Fickett_score
sncrnaNCData$fickett <- sncrna_exon2_nc_fickett_feature$Fickett_score

# Load gerp mammals
functional_protein_exon2_gerp_mammals_feature <- read.csv("../data/GERP_feature/functional-protein-exon2-mammals-gerp-feature.csv", header = TRUE)
functional_protein_exon3_gerp_mammals_feature <- read.csv("../data/GERP_feature/functional-protein-exon3-mammals-gerp-feature.csv", header = TRUE)
functional_lncrna_exon1_gerp_mammals_feature <- read.csv("../data/GERP_feature/functional-lncrna-exon1-mammals-gerp-feature.csv", header = TRUE)
functional_lncrna_exon2_gerp_mammals_feature <- read.csv("../data/GERP_feature/functional-lncrna-exon2-mammals-gerp-feature.csv", header = TRUE)
functional_sncrna_gerp_mammals_feature <- read.csv("../data/GERP_feature/functional-short-ncrna-mammals-gerp-feature.csv", header = TRUE)

protein_exon2_nc_gerp_mammals_feature <- read.csv("../data/GERP_feature/protein-exon2-negative-control-mammals-gerp-feature.csv", header = TRUE)
protein_exon3_nc_gerp_mammals_feature <- read.csv("../data/GERP_feature/protein-exon3-negative-control-mammals-gerp-feature.csv", header = TRUE)
lncrna_exon1_nc_gerp_mammals_feature <- read.csv("../data/GERP_feature/lncrna-exon1-negative-control-mammals-gerp-feature.csv", header = TRUE)
lncrna_exon2_nc_gerp_mammals_feature <- read.csv("../data/GERP_feature/lncrna-exon2-negative-control-mammals-gerp-feature.csv", header = TRUE)
sncrna_exon2_nc_gerp_mammals_feature <- read.csv("../data/GERP_feature/short-ncrna-negative-control-mammals-gerp-feature.csv", header = TRUE)

# Replace gerp mammals
funcProtExon2Data$GERP_91_mammals_max <- functional_protein_exon2_gerp_mammals_feature$GERP_91_mammals_max
funcProtExon3Data$GERP_91_mammals_max <- functional_protein_exon3_gerp_mammals_feature$GERP_91_mammals_max
funcLncrnaExon1Data$GERP_91_mammals_max <- functional_lncrna_exon1_gerp_mammals_feature$GERP_91_mammals_max
funcLncrnaExon2Data$GERP_91_mammals_max <- functional_lncrna_exon2_gerp_mammals_feature$GERP_91_mammals_max
funcSncrnaDataset$GERP_91_mammals_max <- functional_sncrna_gerp_mammals_feature$GERP_91_mammals_max

protExon2NCData$GERP_91_mammals_max <- protein_exon2_nc_gerp_mammals_feature$GERP_91_mammals_max
protExon3NCData$GERP_91_mammals_max <- protein_exon3_nc_gerp_mammals_feature$GERP_91_mammals_max
lncrnaExon1NCData$GERP_91_mammals_max <- lncrna_exon1_nc_gerp_mammals_feature$GERP_91_mammals_max
lncrnaExon2NCData$GERP_91_mammals_max <- lncrna_exon2_nc_gerp_mammals_feature$GERP_91_mammals_max
sncrnaNCData$GERP_91_mammals_max <- sncrna_exon2_nc_gerp_mammals_feature$GERP_91_mammals_max

# Load gerp amniotes
functional_protein_exon2_gerp_amniotes_feature <- read.csv("../data/GERP_feature/functional-protein-exon2-amniotes-gerp-feature.csv", header = TRUE)
functional_protein_exon3_gerp_amniotes_feature <- read.csv("../data/GERP_feature/functional-protein-exon3-amniotes-gerp-feature.csv", header = TRUE)
functional_lncrna_exon1_gerp_amniotes_feature <- read.csv("../data/GERP_feature/functional-lncrna-exon1-amniotes-gerp-feature.csv", header = TRUE)
functional_lncrna_exon2_gerp_amniotes_feature <- read.csv("../data/GERP_feature/functional-lncrna-exon2-amniotes-gerp-feature.csv", header = TRUE)
functional_sncrna_gerp_amniotes_feature <- read.csv("../data/GERP_feature/functional-short-ncrna-amniotes-gerp-feature.csv", header = TRUE)

protein_exon2_nc_gerp_amniotes_feature <- read.csv("../data/GERP_feature/protein-exon2-negative-control-amniotes-gerp-feature.csv", header = TRUE)
protein_exon3_nc_gerp_amniotes_feature <- read.csv("../data/GERP_feature/protein-exon3-negative-control-amniotes-gerp-feature.csv", header = TRUE)
lncrna_exon1_nc_gerp_amniotes_feature <- read.csv("../data/GERP_feature/lncrna-exon1-negative-control-amniotes-gerp-feature.csv", header = TRUE)
lncrna_exon2_nc_gerp_amniotes_feature <- read.csv("../data/GERP_feature/lncrna-exon2-negative-control-amniotes-gerp-feature.csv", header = TRUE)
sncrna_exon2_nc_gerp_amniotes_feature <- read.csv("../data/GERP_feature/short-ncrna-negative-control-amniotes-gerp-feature.csv", header = TRUE)

# Replace gerp amniotes
funcProtExon2Data$GERP_63_amniotes_max <- functional_protein_exon2_gerp_amniotes_feature$GERP_63_amniotes_max
funcProtExon3Data$GERP_63_amniotes_max <- functional_protein_exon3_gerp_amniotes_feature$GERP_63_amniotes_max
funcLncrnaExon1Data$GERP_63_amniotes_max <- functional_lncrna_exon1_gerp_amniotes_feature$GERP_63_amniotes_max
funcLncrnaExon2Data$GERP_63_amniotes_max <- functional_lncrna_exon2_gerp_amniotes_feature$GERP_63_amniotes_max
funcSncrnaDataset$GERP_63_amniotes_max <- functional_sncrna_gerp_amniotes_feature$GERP_63_amniotes_max

protExon2NCData$GERP_63_amniotes_max <- protein_exon2_nc_gerp_amniotes_feature$GERP_63_amniotes_max
protExon3NCData$GERP_63_amniotes_max <- protein_exon3_nc_gerp_amniotes_feature$GERP_63_amniotes_max
lncrnaExon1NCData$GERP_63_amniotes_max <- lncrna_exon1_nc_gerp_amniotes_feature$GERP_63_amniotes_max
lncrnaExon2NCData$GERP_63_amniotes_max <- lncrna_exon2_nc_gerp_amniotes_feature$GERP_63_amniotes_max
sncrnaNCData$GERP_63_amniotes_max <- sncrna_exon2_nc_gerp_amniotes_feature$GERP_63_amniotes_max

# Load dinucleotides
functional_protein_exon2_dinucleotide_feature <- read.csv("../data/datasets/dinucleotide_feature/protein-exon2-dinucleotide-feature.csv", header = TRUE)
functional_protein_exon3_dinucleotide_feature <- read.csv("../data/datasets/dinucleotide_feature/protein-exon3-dinucleotide-feature.csv", header = TRUE)
functional_lncrna_exon1_dinucleotide_feature <- read.csv("../data/datasets/dinucleotide_feature/lncrna-exon1-dinucleotide-feature.csv", header = TRUE)
functional_lncrna_exon2_dinucleotide_feature <- read.csv("../data/datasets/dinucleotide_feature/lncrna-exon2-dinucleotide-feature.csv", header = TRUE)
functional_sncrna_dinucleotide_feature <- read.csv("../data/datasets/dinucleotide_feature/short-ncrna-dinucleotide-feature.csv", header = TRUE)

protein_exon2_nc_dinucleotide_feature <- read.csv("../data/datasets/dinucleotide_feature/protein-exon2-NC-dinucleotide-feature.csv", header = TRUE)
protein_exon3_nc_dinucleotide_feature <- read.csv("../data/datasets/dinucleotide_feature/protein-exon3-NC-dinucleotide-feature.csv", header = TRUE)
lncrna_exon1_nc_dinucleotide_feature <- read.csv("../data/datasets/dinucleotide_feature/lncrna-exon1-NC-dinucleotide-feature.csv", header = TRUE)
lncrna_exon2_nc_dinucleotide_feature <- read.csv("../data/datasets/dinucleotide_feature/lncrna-exon2-NC-dinucleotide-feature.csv", header = TRUE)
sncrna_exon2_nc_dinucleotide_feature <- read.csv("../data/datasets/dinucleotide_feature/short-ncrna-NC-dinucleotide-feature.csv", header = TRUE)








# Replace dinucleotides
for(din in colnames(functional_protein_exon2_dinucleotide_feature)) {
  funcProtExon2Data[[din]] <- functional_protein_exon2_dinucleotide_feature[[din]]
  funcProtExon3Data[[din]] <- functional_protein_exon3_dinucleotide_feature[[din]]
  funcLncrnaExon1Data[[din]]<- functional_lncrna_exon1_dinucleotide_feature[[din]]
  funcLncrnaExon2Data[[din]] <- functional_lncrna_exon2_dinucleotide_feature[[din]]
  funcSncrnaDataset[[din]] <- functional_sncrna_dinucleotide_feature[[din]]
  
  protExon2NCData[[din]] <- protein_exon2_nc_dinucleotide_feature[[din]]
  protExon3NCData[[din]] <- protein_exon3_nc_dinucleotide_feature[[din]]
  lncrnaExon1NCData[[din]] <- lncrna_exon1_nc_dinucleotide_feature[[din]]
  lncrnaExon2NCData[[din]] <- lncrna_exon2_nc_dinucleotide_feature[[din]]
  sncrnaNCData[[din]] <- sncrna_exon2_nc_dinucleotide_feature[[din]]
}

# Write files
write.csv(funcProtExon2Data,"functional-protein-exon2-dataset-features.csv", row.names = FALSE)
write.csv(funcProtExon3Data,"functional-protein-exon3-dataset-features.csv", row.names = FALSE)
write.csv(funcLncrnaExon1Data,"functional-lncrna-exon1-dataset-features.csv", row.names = FALSE)
write.csv(funcLncrnaExon2Data,"functional-lncrna-exon2-dataset-features.csv", row.names = FALSE)
write.csv(funcSncrnaDataset,"functional-short-ncrna-dataset-features.csv", row.names = FALSE)

write.csv(protExon2NCData,"protein-exon2-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(protExon3NCData,"protein-exon3-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(lncrnaExon1NCData,"lncrna-exon1-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(lncrnaExon2NCData,"lncrna-exon2-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(sncrnaNCData,"short-ncrna-negative-control-dataset-features.csv", row.names = FALSE)

# Join all data in single csv file for convenience.
cols1 <- colnames(funcProtExon2Data)
cols2 <- colnames(protExon2NCData)
colsdiff1 <- setdiff(cols1, cols2)
colsdiff2 <- setdiff(cols2, cols1)

feature_matrix <- rbind(data.frame(Dataset="protein-coding-exon2", funcProtExon2Data[, !(names(funcProtExon2Data) %in% colsdiff1)]), 
                        data.frame(Dataset="protein-coding-exon3", funcProtExon3Data[, !(names(funcProtExon3Data) %in% colsdiff1)]),
                        data.frame(Dataset = "lncrna-exon1", funcLncrnaExon1Data[, !(names(funcLncrnaExon1Data) %in% colsdiff1)]), 
                        data.frame(Dataset = "lncrna-exon2", funcLncrnaExon2Data[, !(names(funcLncrnaExon2Data) %in% colsdiff1)]),
                        data.frame(Dataset = "short-ncrna", funcSncrnaDataset[, !(names(funcSncrnaDataset) %in% colsdiff1)]),
                        
                        data.frame(Dataset = "protein-exon2-negative-control", protExon2NCData[, !(names(protExon2NCData) %in% colsdiff2)]), 
                        data.frame(Dataset = "protein-exon3-negative-control", protExon3NCData[, !(names(protExon3NCData) %in% colsdiff2)]),
                        data.frame(Dataset = "lncrna-exon1-negative-control", lncrnaExon1NCData[, !(names(lncrnaExon1NCData) %in% colsdiff2)]), 
                        data.frame(Dataset = "lncrna-exon2-negative-control", lncrnaExon2NCData[, !(names(lncrnaExon2NCData) %in% colsdiff2)]),
                        data.frame(Dataset = "short-ncrna-negative-control", sncrnaNCData[, !(names(sncrnaNCData) %in% colsdiff2)])
                       )


write.csv(feature_matrix, "../data/features/gene_functionality_features_latest1000all_extra.csv")
