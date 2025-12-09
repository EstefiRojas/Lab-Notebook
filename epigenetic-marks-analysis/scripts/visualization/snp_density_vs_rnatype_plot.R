# Load libraries
library(dplyr)
library(ggplot2)

# Load gene functionality data
source("load_gene_functionality_features.R")

# Load RNA type feature
functional_lncrna_exon1_rna_type_feature <- read.csv("../data/datasets/rnatype_feature/functional-lncrna-exon1-rnatype-hgncid-feature.csv", header = TRUE)
functional_lncrna_exon2_rna_type_feature <- read.csv("../data/datasets/rnatype_feature/functional-lncrna-exon2-rnatype-hgncid-feature.csv", header = TRUE)
functional_sncrna_rna_type_feature <- read.csv("../data/datasets/rnatype_feature/functional-short-ncrna-rnatype-hgncid-feature.csv", header = TRUE)

funcSncrnaDataset$hgnc_id <- functional_sncrna_rna_type_feature$hgnc_id

funcLncrnaExon1Data$hgnc_id <- functional_lncrna_exon1_rna_type_feature$hgnc_id

funcLncrnaExon2Data$hgnc_id <- functional_lncrna_exon2_rna_type_feature$hgnc_id

summary(funcSncrnaDataset$gnomAD_SNP_density)
# Fill RNA type
funcProtExon2Data$rna_type <- "mRNA"
funcProtExon3Data$rna_type <- "mRNA"
funcLncrnaExon1Data$rna_type <- functional_lncrna_exon1_rna_type_feature$rna_type
funcLncrnaExon1Data$rna_type <- "lncRNA"
funcLncrnaExon2Data$rna_type <- functional_lncrna_exon2_rna_type_feature$rna_type
funcLncrnaExon2Data$rna_type <- "lncRNA"
funcSncrnaDataset$rna_type <- functional_sncrna_rna_type_feature$rna_type
unique(funcLncrnaExon1Data$rna_type)
sncrnaNCData$rna_type <- "sncRNA(-)"
unique(funcSncrnaDataset$rna_type)

funcLncrnaExon1Data$GeneName <- functional_lncrna_exon1_rna_type_feature$gene_name
funcLncrnaExon2Data$GeneName <- functional_lncrna_exon2_rna_type_feature$gene_name
funcSncrnaDataset$GeneName <- functional_sncrna_rna_type_feature$gene_name

#Filter intereseting sncRNAs
# ncRNAs
ncrna_dataset <- funcSncrnaDataset %>%
  filter(rna_type=="ncRNA")
summary(ncrna_dataset$gnomAD_SNP_density)

# pre_miRNA
pre_mirna_dataset <- funcSncrnaDataset %>%
  filter(rna_type=="pre_miRNA")
pre_mirna_dataset

# tRNA
trna_dataset <- funcSncrnaDataset %>%
  filter(rna_type=="tRNA")
summary(trna_dataset$gnomAD_SNP_density)

# scaRNA
trna_dataset <- funcSncrnaDataset %>%
  filter(rna_type=="scaRNA")
trna_dataset

# snRNA
snrna_dataset <- funcSncrnaDataset %>%
  filter(rna_type=="snRNA")
summary(snrna_dataset$gnomAD_SNP_density)

# vault_RNA
vaultrna_dataset <- funcSncrnaDataset %>%
  filter(rna_type=="vault_RNA")
summary(vaultrna_dataset$gnomAD_SNP_density)

# SNP density greater that 1 (more SNPs than nucleotides in the sequence)
hyper_snp_counts <- funcSncrnaDataset %>% 
  select(Chromosome,Start,End,hgnc_id,GeneName,SNP_density,rna_type) %>%
  filter(SNP_density > 1)
write.csv(hyper_snp_counts, "../data/hyper_SNP_count_sncRNA_dataset.csv", row.names = FALSE)

hyper_snp_counts_lncrnaex1 <- funcLncrnaExon1Data %>%
  select(Chromosome,Start,End,hgnc_id,GeneName,SNP_density,rna_type) %>%
  filter(SNP_density > 1)
write.csv(hyper_snp_counts_lncrnaex1, "../data/hyper_SNP_count_lncRNAex1_dataset.csv", row.names = FALSE)

hyper_snp_counts_lncrnaex2 <- funcLncrnaExon2Data %>%
  select(Chromosome,Start,End,hgnc_id,GeneName,SNP_density,rna_type) %>%
  filter(SNP_density > 1)
write.csv(hyper_snp_counts_lncrnaex2, "../data/hyper_SNP_count_lncRNAex2_dataset.csv")



allDatasets <- as.data.frame(rbind(funcProtExon2Data[,c("SNP density","rna_type")],
                                   funcProtExon3Data[,c("SNP density","rna_type")],
                                   funcLncrnaExon1Data[,c("SNP density","rna_type")],
                                   funcLncrnaExon2Data[,c("SNP density","rna_type")],
                                   funcSncrnaDataset[,c("SNP density","rna_type")],
                                   sncrnaNCData[,c("SNP density","rna_type")]))

unique(allDatasets$rna_type)
allDatasets$rna_type <- factor(allDatasets$rna_type,
                               levels = unique(allDatasets$rna_type),
                               labels = c("mRNA", 
                                          "lncRNA",
                                          "pre-miRNA",
                                          "snoRNA",
                                          "snRNA",
                                          "scaRNA",
                                          "tRNA",
                                          "snaR RNA",
                                          "vault RNA",
                                          "sncRNA(-)"))

funcSncrnaDataset$rna_type <- factor(funcSncrnaDataset$rna_type,
                                     levels = unique(funcSncrnaDataset$rna_type))
summary(allDatasets$rna_type)


# Calculate the median SNP density for each rna_type
medians <- allDatasets %>%
  group_by(rna_type) %>%
  summarise(median_density = median(`SNP density`, na.rm = TRUE)) %>%
  arrange(median_density)

# Reorder the rna_type factor based on the calculated medians
allDatasets$rna_type <- factor(allDatasets$rna_type, levels = medians$rna_type)


ggplot(allDatasets, aes(x = rna_type, y = `SNP density`, fill = `rna_type`)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.0, outliers=FALSE, position = position_dodge(width = 0.9), width=0.2, size = 1) +
  #facet_wrap(~ rna_type, scales = "free") + 
  labs(title = "SNP density Distribution by RNA type", y = "SNP Density") +
  theme_minimal(base_size = 56) +
  theme(
    axis.text.x = element_text(size = 20, hjust = 0.5),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 36),  # Increase y-axis title size
    legend.position = "none",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend text size
    plot.title = element_text(size = 52, hjust = 0.5, margin = margin(b = 60)),  # Increase plot title size and center it, increase margin
    plot.subtitle = element_blank(),
    axis.line.x = element_blank(),          # Remove x-axis line
    axis.ticks.x = element_blank(),         # Remove x-axis ticks
    panel.grid.major.x = element_blank(),   # Remove major grid lines along x-axis
    panel.grid.minor.x = element_blank(),    # Remove minor grid lines along x-axis
    strip.text = element_blank()            # Remove facet titles
  ) +
  scale_fill_manual(values = c("#F4A582FF","#FCF2F1","#FAECEA","#e37b88FF","#56bdfcFF", "#F8E6E3", "#F7DFDC", "#F5D9D4","#EBB0A6","#E4988B")) +
  coord_cartesian(ylim = c(0, 3)) # Zoom into an specific area without removing data points

ggsave("RNATypeVsSnpDensityDistribution.png",path = "../results/latest1000all/rna_type/", scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

