load_dinucleotide_features <- function() {
  source("scripts/config.R")
  library(dplyr)
  library(tidyr)

  # Load existing features
  funcProtExon2Data <- data.frame(read.csv(FUNC_PROT_EXON2_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-GeneID)
  funcProtExon3Data <- data.frame(read.csv(FUNC_PROT_EXON3_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-GeneID)
  funcLncrnaExon1Data <- data.frame(read.csv(FUNC_LNCRNA_EXON1_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-GeneID)
  funcLncrnaExon2Data <- data.frame(read.csv(FUNC_LNCRNA_EXON2_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-GeneID)
  funcSncrnaData <- data.frame(read.csv(FUNC_SNCRNA_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-GeneID)

  protExon2NCData <- data.frame(read.csv(NC_PROT_EXON2_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-DistanceGene)
  protExon3NCData <- data.frame(read.csv(NC_PROT_EXON3_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-DistanceGene)
  lncrnaExon1NCData <- data.frame(read.csv(NC_LNCRNA_EXON1_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-DistanceGene)
  lncrnaExon2NCData <- data.frame(read.csv(NC_LNCRNA_EXON2_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-DistanceGene)
  sncrnaNCData <- data.frame(read.csv(NC_SNCRNA_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-DistanceGene)

  # Load existing dinucleotide features
  funcProtExon2DinData <- data.frame(read.csv(PROT_EXON2_DIN_FEATURES_FILE, header=TRUE))
  funcProtExon3DinData <- data.frame(read.csv(PROT_EXON3_DIN_FEATURES_FILE, header=TRUE))
  funcLncrnaExon1DinData <- data.frame(read.csv(LNCRNA_EXON1_DIN_FEATURES_FILE, header=TRUE))
  funcLncrnaExon2DinData <- data.frame(read.csv(LNCRNA_EXON2_DIN_FEATURES_FILE, header=TRUE))
  funcSncrnaDinData <- data.frame(read.csv(SNCRNA_DIN_FEATURES_FILE, header=TRUE))

  protExon2NCDinData <- data.frame(read.csv(NC_PROT_EXON2_DIN_FEATURES_FILE, header=TRUE))
  protExon3NCDinData <- data.frame(read.csv(NC_PROT_EXON3_DIN_FEATURES_FILE, header=TRUE))
  lncrnaExon1NCDinData <- data.frame(read.csv(NC_LNCRNA_EXON1_DIN_FEATURES_FILE, header=TRUE))
  lncrnaExon2NCDinData <- data.frame(read.csv(NC_LNCRNA_EXON2_DIN_FEATURES_FILE, header=TRUE))
  sncrnaNCDinData <- data.frame(read.csv(NC_SNCRNA_DIN_FEATURES_FILE, header=TRUE))

  # Replace dinucleotides
  for(din in colnames(funcProtExon2DinData)) {
    funcProtExon2Data[[din]] <- funcProtExon2DinData[[din]]
    funcProtExon3Data[[din]] <- funcProtExon3DinData[[din]]
    funcLncrnaExon1Data[[din]]<- funcLncrnaExon1DinData[[din]]
    funcLncrnaExon2Data[[din]] <- funcLncrnaExon2DinData[[din]]
    funcSncrnaData[[din]] <- funcSncrnaDinData[[din]]
    
    protExon2NCData[[din]] <- protExon2NCDinData[[din]]
    protExon3NCData[[din]] <- protExon3NCDinData[[din]]
    lncrnaExon1NCData[[din]] <- lncrnaExon1NCDinData[[din]]
    lncrnaExon2NCData[[din]] <- lncrnaExon2NCDinData[[din]]
    sncrnaNCData[[din]] <- sncrnaNCDinData[[din]]
  }

  # Join all data in single csv file for convenience.
  cols1 <- colnames(funcProtExon2Data)
  cols2 <- colnames(protExon2NCData)
  colsdiff1 <- setdiff(cols1, cols2)
  colsdiff2 <- setdiff(cols2, cols1)

  feature_matrix <- rbind(data.frame(Dataset="protein-coding-exon2", funcProtExon2Data[, !(names(funcProtExon2Data) %in% colsdiff1)]), 
                                     data.frame(Dataset="protein-coding-exon3", funcProtExon3Data[, !(names(funcProtExon3Data) %in% colsdiff1)]),
                                     data.frame(Dataset = "lncrna-exon1", funcLncrnaExon1Data[, !(names(funcLncrnaExon1Data) %in% colsdiff1)]), 
                                     data.frame(Dataset = "lncrna-exon2", funcLncrnaExon2Data[, !(names(funcLncrnaExon2Data) %in% colsdiff1)]),
                                     data.frame(Dataset = "short-ncrna", funcSncrnaData[, !(names(funcSncrnaData) %in% colsdiff1)]),
                                     
                                     data.frame(Dataset = "protein-exon2-negative-control", protExon2NCData[, !(names(protExon2NCData) %in% colsdiff2)]), 
                                     data.frame(Dataset = "protein-exon3-negative-control", protExon3NCData[, !(names(protExon3NCData) %in% colsdiff2)]),
                                     data.frame(Dataset = "lncrna-exon1-negative-control", lncrnaExon1NCData[, !(names(lncrnaExon1NCData) %in% colsdiff2)]), 
                                     data.frame(Dataset = "lncrna-exon2-negative-control", lncrnaExon2NCData[, !(names(lncrnaExon2NCData) %in% colsdiff2)]),
                                     data.frame(Dataset = "short-ncrna-negative-control", sncrnaNCData[, !(names(sncrnaNCData) %in% colsdiff2)])
  )

  return(feature_matrix)
}