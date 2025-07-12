#Load existing histone features
#source("load_gene_functionality_features.R")
funcProtExon2HistoneData <- data.frame(read.csv("pasted_histone_outputs/combined_protein-exon2.csv", header=TRUE))
funcProtExon3HistoneData <- data.frame(read.csv("pasted_histone_outputs/combined_protein-exon3.csv", header=TRUE))
funcLncrnaExon1HistoneData <- data.frame(read.csv("pasted_histone_outputs/combined_lncrna-exon1.csv", header=TRUE))
funcLncrnaExon2HistoneData <- data.frame(read.csv("pasted_histone_outputs/combined_lncrna-exon2.csv", header=TRUE))
funcSncrnaHistoneData <- data.frame(read.csv("pasted_histone_outputs/combined_short-ncrna.csv", header=TRUE))

protExon2NCHistoneData <- data.frame(read.csv("pasted_histone_outputs/combined_protein-exon2-NC.csv", header=TRUE))
protExon3NCHistoneData <- data.frame(read.csv("pasted_histone_outputs/combined_protein-exon3-NC.csv", header=TRUE))
lncrnaExon1NCHistoneData <- data.frame(read.csv("pasted_histone_outputs/combined_lncrna-exon1-NC.csv", header=TRUE))
lncrnaExon2NCHistoneData <- data.frame(read.csv("pasted_histone_outputs/combined_lncrna-exon2-NC.csv", header=TRUE))
sncrnaNCHistoneData <- data.frame(read.csv("pasted_histone_outputs/combined_short-ncrna-NC.csv", header=TRUE))

# Replace histones
#H3K27ac
funcProtExon2Data$H3K27ac_MaxScaledSignal <- funcProtExon2HistoneData$H3K27ac_MaxScaledSignal
funcProtExon3Data$H3K27ac_MaxScaledSignal <- funcProtExon3HistoneData$H3K27ac_MaxScaledSignal
funcLncrnaExon1Data$H3K27ac_MaxScaledSignal <- funcLncrnaExon1HistoneData$H3K27ac_MaxScaledSignal
funcLncrnaExon2Data$H3K27ac_MaxScaledSignal <- funcLncrnaExon2HistoneData$H3K27ac_MaxScaledSignal
funcSncrnaDataset$H3K27ac_MaxScaledSignal <- funcSncrnaHistoneData$H3K27ac_MaxScaledSignal

protExon2NCData$H3K27ac_MaxScaledSignal <- protExon2NCHistoneData$H3K27ac_MaxScaledSignal
protExon3NCData$H3K27ac_MaxScaledSignal <- protExon3NCHistoneData$H3K27ac_MaxScaledSignal
lncrnaExon1NCData$H3K27ac_MaxScaledSignal <- lncrnaExon1NCHistoneData$H3K27ac_MaxScaledSignal
lncrnaExon2NCData$H3K27ac_MaxScaledSignal <- lncrnaExon2NCHistoneData$H3K27ac_MaxScaledSignal
sncrnaNCData$H3K27ac_MaxScaledSignal <- sncrnaNCHistoneData$H3K27ac_MaxScaledSignal

# H3K36me3
funcProtExon2Data$H3K36me3_MaxScaledSignal <- funcProtExon2HistoneData$H3K36me3_MaxScaledSignal
funcProtExon3Data$H3K36me3_MaxScaledSignal <- funcProtExon3HistoneData$H3K36me3_MaxScaledSignal
funcLncrnaExon1Data$H3K36me3_MaxScaledSignal <- funcLncrnaExon1HistoneData$H3K36me3_MaxScaledSignal
funcLncrnaExon2Data$H3K36me3_MaxScaledSignal <- funcLncrnaExon2HistoneData$H3K36me3_MaxScaledSignal
funcSncrnaDataset$H3K36me3_MaxScaledSignal <- funcSncrnaHistoneData$H3K36me3_MaxScaledSignal

protExon2NCData$H3K36me3_MaxScaledSignal <- protExon2NCHistoneData$H3K36me3_MaxScaledSignal
protExon3NCData$H3K36me3_MaxScaledSignal <- protExon3NCHistoneData$H3K36me3_MaxScaledSignal
lncrnaExon1NCData$H3K36me3_MaxScaledSignal <- lncrnaExon1NCHistoneData$H3K36me3_MaxScaledSignal
lncrnaExon2NCData$H3K36me3_MaxScaledSignal <- lncrnaExon2NCHistoneData$H3K36me3_MaxScaledSignal
sncrnaNCData$H3K36me3_MaxScaledSignal <- sncrnaNCHistoneData$H3K36me3_MaxScaledSignal

#H3K79me2
funcProtExon2Data$H3K79me2_MaxScaledSignal <- funcProtExon2HistoneData$H3K79me2_MaxScaledSignal
funcProtExon3Data$H3K79me2_MaxScaledSignal <- funcProtExon3HistoneData$H3K79me2_MaxScaledSignal
funcLncrnaExon1Data$H3K79me2_MaxScaledSignal <- funcLncrnaExon1HistoneData$H3K79me2_MaxScaledSignal
funcLncrnaExon2Data$H3K79me2_MaxScaledSignal <- funcLncrnaExon2HistoneData$H3K79me2_MaxScaledSignal
funcSncrnaDataset$H3K79me2_MaxScaledSignal <- funcSncrnaHistoneData$H3K79me2_MaxScaledSignal

protExon2NCData$H3K79me2_MaxScaledSignal <- protExon2NCHistoneData$H3K79me2_MaxScaledSignal
protExon3NCData$H3K79me2_MaxScaledSignal <- protExon3NCHistoneData$H3K79me2_MaxScaledSignal
lncrnaExon1NCData$H3K79me2_MaxScaledSignal <- lncrnaExon1NCHistoneData$H3K79me2_MaxScaledSignal
lncrnaExon2NCData$H3K79me2_MaxScaledSignal <- lncrnaExon2NCHistoneData$H3K79me2_MaxScaledSignal
sncrnaNCData$H3K79me2_MaxScaledSignal <- sncrnaNCHistoneData$H3K79me2_MaxScaledSignal


#H3K9ac
funcProtExon2Data$H3K9ac_MaxScaledSignal <- funcProtExon2HistoneData$H3K9ac_MaxScaledSignal
funcProtExon3Data$H3K9ac_MaxScaledSignal <- funcProtExon3HistoneData$H3K9ac_MaxScaledSignal
funcLncrnaExon1Data$H3K9ac_MaxScaledSignal <- funcLncrnaExon1HistoneData$H3K9ac_MaxScaledSignal
funcLncrnaExon2Data$H3K9ac_MaxScaledSignal <- funcLncrnaExon2HistoneData$H3K9ac_MaxScaledSignal
funcSncrnaDataset$H3K9ac_MaxScaledSignal <- funcSncrnaHistoneData$H3K9ac_MaxScaledSignal

protExon2NCData$H3K9ac_MaxScaledSignal <- protExon2NCHistoneData$H3K9ac_MaxScaledSignal
protExon3NCData$H3K9ac_MaxScaledSignal <- protExon3NCHistoneData$H3K9ac_MaxScaledSignal
lncrnaExon1NCData$H3K9ac_MaxScaledSignal <- lncrnaExon1NCHistoneData$H3K9ac_MaxScaledSignal
lncrnaExon2NCData$H3K9ac_MaxScaledSignal <- lncrnaExon2NCHistoneData$H3K9ac_MaxScaledSignal
sncrnaNCData$H3K9ac_MaxScaledSignal <- sncrnaNCHistoneData$H3K9ac_MaxScaledSignal


# Load chromatin accessibility feature
funcProtExon2chrmacc <- read.csv("../data/datasets/chrm_acc_feature/narrowPeak/chrm_acc_protein-exon2-chrm_acc-feature.csv", header = TRUE)
funcProtExon3chrmacc <- read.csv("../data/datasets/chrm_acc_feature/narrowPeak/chrm_acc_protein-exon3-chrm_acc-feature.csv", header = TRUE)
funcLncrnaExon1chrmacc <- read.csv("../data/datasets/chrm_acc_feature/narrowPeak/chrm_acc_lncrna-exon1-chrm_acc-feature.csv", header = TRUE)
funcLncrnaExon2chrmacc <- read.csv("../data/datasets/chrm_acc_feature/narrowPeak/chrm_acc_lncrna-exon2-chrm_acc-feature.csv", header = TRUE)
funcSncrnachrmacc <- read.csv("../data/datasets/chrm_acc_feature/narrowPeak/chrm_acc_short-ncrna-chrm_acc-feature.csv", header = TRUE)
negProtExon2chrmacc <- read.csv("../data/datasets/chrm_acc_feature/narrowPeak/chrm_acc_protein-exon2-NC-chrm_acc-feature.csv", header = TRUE)
negProtExon3chrmacc <- read.csv("../data/datasets/chrm_acc_feature/narrowPeak/chrm_acc_protein-exon3-NC-chrm_acc-feature.csv", header = TRUE)
negLncrnaExon1chrmacc <- read.csv("../data/datasets/chrm_acc_feature/narrowPeak/chrm_acc_lncrna-exon1-NC-chrm_acc-feature.csv", header = TRUE)
negLncrnaExon2chrmacc <- read.csv("../data/datasets/chrm_acc_feature/narrowPeak/chrm_acc_lncrna-exon2-NC-chrm_acc-feature.csv", header = TRUE)
negSncrnachrmacc <- read.csv("../data/datasets/chrm_acc_feature/narrowPeak/chrm_acc_short-ncrna-NC-chrm_acc-feature.csv", header = TRUE)

funcProtExon2HistoneData$chrm_acc_MaxScaledSignal <- funcProtExon2chrmacc$chrm_acc_MaxScaledSignal
funcProtExon3HistoneData$chrm_acc_MaxScaledSignal <- funcProtExon3chrmacc$chrm_acc_MaxScaledSignal
funcLncrnaExon1HistoneData$chrm_acc_MaxScaledSignal <- funcLncrnaExon1chrmacc$chrm_acc_MaxScaledSignal
funcLncrnaExon2HistoneData$chrm_acc_MaxScaledSignal <- funcLncrnaExon2chrmacc$chrm_acc_MaxScaledSignal
funcSncrnaHistoneData$chrm_acc_MaxScaledSignal <- funcSncrnachrmacc$chrm_acc_MaxScaledSignal

protExon2NCHistoneData$chrm_acc_MaxScaledSignal <- negProtExon2chrmacc$chrm_acc_MaxScaledSignal
protExon3NCHistoneData$chrm_acc_MaxScaledSignal <- negProtExon3chrmacc$chrm_acc_MaxScaledSignal
lncrnaExon1NCHistoneData$chrm_acc_MaxScaledSignal <- negLncrnaExon1chrmacc$chrm_acc_MaxScaledSignal
lncrnaExon2NCHistoneData$chrm_acc_MaxScaledSignal <- negLncrnaExon2chrmacc$chrm_acc_MaxScaledSignal
sncrnaNCHistoneData$chrm_acc_MaxScaledSignal <- negSncrnachrmacc$chrm_acc_MaxScaledSignal

# Load methylome feature
funcProtExon2methylome <- read.csv("../data/datasets/methylome_feature/protein-exon2-methylome-feature.csv", header = TRUE)
funcProtExon3methylome <- read.csv("../data/datasets/methylome_feature/protein-exon3-methylome-feature.csv", header = TRUE)
funcLncrnaExon1methylome <- read.csv("../data/datasets/methylome_feature/lncrna-exon1-methylome-feature.csv", header = TRUE)
funcLncrnaExon2methylome <- read.csv("../data/datasets/methylome_feature/lncrna-exon2-methylome-feature.csv", header = TRUE)
funcSncrnamethylome <- read.csv("../data/datasets/methylome_feature/short-ncrna-methylome-feature.csv", header = TRUE)
negProtExon2methylome <- read.csv("../data/datasets/methylome_feature/protein-exon2-NC-methylome-feature.csv", header = TRUE)
negProtExon3methylome <- read.csv("../data/datasets/methylome_feature/protein-exon3-NC-methylome-feature.csv", header = TRUE)
negLncrnaExon1methylome <- read.csv("../data/datasets/methylome_feature/lncrna-exon1-NC-methylome-feature.csv", header = TRUE)
negLncrnaExon2methylome <- read.csv("../data/datasets/methylome_feature/lncrna-exon2-NC-methylome-feature.csv", header = TRUE)
negSncrnamethylome <- read.csv("../data/datasets/methylome_feature/short-ncrna-NC-methylome-feature.csv", header = TRUE)

funcProtExon2HistoneData$methylome <- funcProtExon2methylome$methylome
funcProtExon3HistoneData$methylome <- funcProtExon3methylome$methylome
funcLncrnaExon1HistoneData$methylome <- funcLncrnaExon1methylome$methylome
funcLncrnaExon2HistoneData$methylome <- funcLncrnaExon2methylome$methylome
funcSncrnaHistoneData$methylome <- funcSncrnamethylome$methylome

protExon2NCHistoneData$methylome <- negProtExon2methylome$methylome
protExon3NCHistoneData$methylome <- negProtExon3methylome$methylome
lncrnaExon1NCHistoneData$methylome <- negLncrnaExon1methylome$methylome
lncrnaExon2NCHistoneData$methylome <- negLncrnaExon2methylome$methylome
sncrnaNCHistoneData$methylome <- negSncrnamethylome$methylome

# Join all data in single csv file for convenience.
cols1 <- colnames(funcProtExon2HistoneData)
cols2 <- colnames(protExon2NCHistoneData)
colsdiff1 <- setdiff(cols1, cols2)
colsdiff2 <- setdiff(cols2, cols1)

feature_matrix_histones <- rbind(data.frame(Dataset="protein-coding-exon2", funcProtExon2HistoneData[, !(names(funcProtExon2HistoneData) %in% colsdiff1)]), 
                        data.frame(Dataset="protein-coding-exon3", funcProtExon3HistoneData[, !(names(funcProtExon3HistoneData) %in% colsdiff1)]),
                        data.frame(Dataset = "lncrna-exon1", funcLncrnaExon1HistoneData[, !(names(funcLncrnaExon1HistoneData) %in% colsdiff1)]), 
                        data.frame(Dataset = "lncrna-exon2", funcLncrnaExon2HistoneData[, !(names(funcLncrnaExon2HistoneData) %in% colsdiff1)]),
                        data.frame(Dataset = "short-ncrna", funcSncrnaHistoneData[, !(names(funcSncrnaHistoneData) %in% colsdiff1)]),
                        
                        data.frame(Dataset = "protein-exon2-negative-control", protExon2NCHistoneData[, !(names(protExon2NCHistoneData) %in% colsdiff2)]), 
                        data.frame(Dataset = "protein-exon3-negative-control", protExon3NCHistoneData[, !(names(protExon3NCHistoneData) %in% colsdiff2)]),
                        data.frame(Dataset = "lncrna-exon1-negative-control", lncrnaExon1NCHistoneData[, !(names(lncrnaExon1NCHistoneData) %in% colsdiff2)]), 
                        data.frame(Dataset = "lncrna-exon2-negative-control", lncrnaExon2NCHistoneData[, !(names(lncrnaExon2NCHistoneData) %in% colsdiff2)]),
                        data.frame(Dataset = "short-ncrna-negative-control", sncrnaNCHistoneData[, !(names(sncrnaNCHistoneData) %in% colsdiff2)])
)


write.csv(feature_matrix, "pasted_histone_outputs/histone_features_matrix.csv")
write.csv(feature_matrix, "../data/features/gene_functionality_features_latest1000all_extra.csv")


# Save individual dataset files:
#features to select:
FEATURES<-c("ID","Functional","Chromosome","Start","End","Sequence","GeneID","GC_percentage","GA","CpG","GG","TA",
            "phyloP_max_241w","phyloP_max_100w","RPKM_tissue","RPKM_primary.cell","copy_number","repeat_distance",
            "Interaction_ave","coding_potential","RNAalifold_score","Max_covariance","MFE","SNP_density","MAF_avg",
            "Random","accessibility","fickett","GERP_91_mammals_max","GERP_63_amniotes_max","AA","AC","AG","AT","CA",
            "CC","CT","GC","GT","TC","TG","TT","methylome","lowComplexity_density","H3K9ac_MaxScaledSignal",
            "H3K79me2_MaxScaledSignal","chrm_acc_MaxScaledSignal")
FEATURES_NC<-c("ID","Functional","Chromosome","Start","End","Sequence","DistanceGene","GC_percentage","GA","CpG","GG","TA",
               "phyloP_max_241w","phyloP_max_100w","RPKM_tissue","RPKM_primary.cell","copy_number","repeat_distance",
               "Interaction_ave","coding_potential","RNAalifold_score","Max_covariance","MFE","SNP_density","MAF_avg",
               "Random","accessibility","fickett","GERP_91_mammals_max","GERP_63_amniotes_max","AA","AC","AG","AT","CA",
               "CC","CT","GC","GT","TC","TG","TT","methylome","lowComplexity_density","H3K9ac_MaxScaledSignal",
               "H3K79me2_MaxScaledSignal","chrm_acc_MaxScaledSignal")

FEATURES_Z<-c("GC_percentage","GA","CpG","GG","TA",
            "phyloP_max_241w","phyloP_max_100w","RPKM_tissue","RPKM_primary.cell","copy_number","repeat_distance",
            "Interaction_ave","coding_potential","RNAalifold_score","Max_covariance","MFE","SNP_density","MAF_avg",
            "Random","accessibility","fickett","GERP_91_mammals_max","GERP_63_amniotes_max","AA","AC","AG","AT","CA",
            "CC","CT","GC","GT","TC","TG","TT","methylome","lowComplexity_density","H3K9ac_MaxScaledSignal",
            "H3K79me2_MaxScaledSignal","chrm_acc_MaxScaledSignal")
FEATURES_NC_Z<-c("GC_percentage","GA","CpG","GG","TA",
               "phyloP_max_241w","phyloP_max_100w","RPKM_tissue","RPKM_primary.cell","copy_number","repeat_distance",
               "Interaction_ave","coding_potential","RNAalifold_score","Max_covariance","MFE","SNP_density","MAF_avg",
               "Random","accessibility","fickett","GERP_91_mammals_max","GERP_63_amniotes_max","AA","AC","AG","AT","CA",
               "CC","CT","GC","GT","TC","TG","TT","methylome","lowComplexity_density","H3K9ac_MaxScaledSignal",
               "H3K79me2_MaxScaledSignal","chrm_acc_MaxScaledSignal")

setwd("/Volumes/ADATA HD710 PRO/Downloads/Estefi/Otago University/Experiment/newGeneFunctionality/gene-functionality/scripts")
write.csv(funcProtExon2Data %>% dplyr::select(all_of(FEATURES)), 
          file = "../results/functional-protein-exon2-dataset-features.csv", row.names = FALSE)
write.csv(funcProtExon3Data %>% dplyr::select(all_of(FEATURES)), 
          file = "../results/functional-protein-exon3-dataset-features.csv", row.names = FALSE)
write.csv(funcLncrnaExon1Data %>% dplyr::select(all_of(FEATURES)), 
          file = "../results/functional-lncrna-exon1-dataset-features.csv", row.names = FALSE)
write.csv(funcLncrnaExon2Data %>% dplyr::select(all_of(FEATURES)), 
          file = "../results/functional-lncrna-exon2-dataset-features.csv", row.names = FALSE)
write.csv(funcSncrnaDataset %>% dplyr::select(all_of(FEATURES)), 
          file = "../results/functional-short-ncrna-dataset-features.csv", row.names = FALSE)


write.csv(protExon2NCData %>% dplyr::select(all_of(FEATURES_NC)), 
          file = "../results/protein-exon2-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(protExon3NCData %>% dplyr::select(all_of(FEATURES_NC)), 
          file = "../results/protein-exon3-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(lncrnaExon1NCData %>% dplyr::select(all_of(FEATURES_NC)), 
          file = "../results/lncrna-exon1-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(lncrnaExon2NCData %>% dplyr::select(all_of(FEATURES_NC)), 
          file = "../results/lncrna-exon2-negative-control-dataset-features.csv", row.names = FALSE)
write.csv(sncrnaNCData %>% dplyr::select(all_of(FEATURES_NC)), 
          file = "../results/short-ncrna-negative-control-dataset-features.csv", row.names = FALSE)

unique(feature_matrix$Dataset)
# Save individual z-score data-set files:
# Unpack zscores
funcProtExon2Zscores <- protein_functional_z_scores %>% 
  filter(Dataset=="protein-coding-exon2") %>%
  dplyr::select(-Dataset,-H3K27ac_MaxScaledSignal,-H3K36me3_MaxScaledSignal)
funcProtExon3Zscores <- protein_functional_z_scores %>% 
    filter(Dataset=="protein-coding-exon3") %>%
  dplyr::select(-Dataset,-H3K27ac_MaxScaledSignal,-H3K36me3_MaxScaledSignal)
funcLncrnaExon1Zscores <- lncrna_functional_z_scores %>% 
    filter(Dataset=="lncrna-exon1") %>%
  dplyr::select(-Dataset,-H3K27ac_MaxScaledSignal,-H3K36me3_MaxScaledSignal)
funcLncrnaExon2Zscores <- lncrna_functional_z_scores %>% 
    filter(Dataset=="lncrna-exon2") %>%
  dplyr::select(-Dataset,-H3K27ac_MaxScaledSignal,-H3K36me3_MaxScaledSignal)
funcSncrnaZscores <- sncrna_functional_z_scores %>% 
  filter(Dataset=="short-ncrna") %>%
  dplyr::select(-Dataset,-H3K27ac_MaxScaledSignal,-H3K36me3_MaxScaledSignal)

protExon2NCZscores <- protein_negative_z_scores %>% 
  filter(Dataset=="protein-exon2-negative-control") %>%
  dplyr::select(-Dataset,-H3K27ac_MaxScaledSignal,-H3K36me3_MaxScaledSignal)
protExon3NCZscores <- protein_negative_z_scores %>% 
  filter(Dataset=="protein-exon3-negative-control") %>%
  dplyr::select(-Dataset,-H3K27ac_MaxScaledSignal,-H3K36me3_MaxScaledSignal)
lncrnaExon1NCZscores <- lncrna_negative_z_scores %>% 
  filter(Dataset=="lncrna-exon1-negative-control") %>%
  dplyr::select(-Dataset,-H3K27ac_MaxScaledSignal,-H3K36me3_MaxScaledSignal)
lncrnaExon2NCZscores <- lncrna_negative_z_scores %>% 
  filter(Dataset=="lncrna-exon2-negative-control") %>%
  dplyr::select(-Dataset,-H3K27ac_MaxScaledSignal,-H3K36me3_MaxScaledSignal)
sncrnaNCZscores <- sncrna_negative_z_scores %>% 
  filter(Dataset=="short-ncrna-negative-control") %>%
  dplyr::select(-Dataset,-H3K27ac_MaxScaledSignal,-H3K36me3_MaxScaledSignal)

# Save Z-scores:
write.csv(funcProtExon2Zscores %>% dplyr::select(all_of(FEATURES_Z)), 
          file = "../results/z-scores/functional-protein-exon2-dataset-zscores.csv", row.names = FALSE)
write.csv(funcProtExon3Zscores %>% dplyr::select(all_of(FEATURES_Z)), 
          file = "../results/z-scores/functional-protein-exon3-dataset-zscores.csv", row.names = FALSE)
write.csv(funcLncrnaExon1Zscores %>% dplyr::select(all_of(FEATURES_Z)), 
          file = "../results/z-scores/functional-lncrna-exon1-dataset-zscores.csv", row.names = FALSE)
write.csv(funcLncrnaExon2Zscores %>% dplyr::select(all_of(FEATURES_Z)), 
          file = "../results/z-scores/functional-lncrna-exon2-dataset-zscores.csv", row.names = FALSE)
write.csv(funcSncrnaZscores %>% dplyr::select(all_of(FEATURES_Z)), 
          file = "../results/z-scores/functional-short-ncrna-dataset-zscores.csv", row.names = FALSE)


write.csv(protExon2NCZscores %>% dplyr::select(all_of(FEATURES_NC_Z)), 
          file = "../results/z-scores/protein-exon2-negative-control-dataset-zscores.csv", row.names = FALSE)
write.csv(protExon3NCZscores %>% dplyr::select(all_of(FEATURES_NC_Z)), 
          file = "../results/z-scores/protein-exon3-negative-control-dataset-zscores.csv", row.names = FALSE)
write.csv(lncrnaExon1NCZscores %>% dplyr::select(all_of(FEATURES_NC_Z)), 
          file = "../results/z-scores/lncrna-exon1-negative-control-dataset-zscores.csv", row.names = FALSE)
write.csv(lncrnaExon2NCZscores %>% dplyr::select(all_of(FEATURES_NC_Z)), 
          file = "../results/z-scores/lncrna-exon2-negative-control-dataset-zscores.csv", row.names = FALSE)
write.csv(sncrnaNCZscores %>% dplyr::select(all_of(FEATURES_NC_Z)), 
          file = "../results/z-scores/short-ncrna-negative-control-dataset-zscores.csv", row.names = FALSE)
