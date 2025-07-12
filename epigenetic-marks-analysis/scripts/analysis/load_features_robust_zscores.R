# Load robust-z-score datasets
prot_zscores <- data.frame(Dataset="protein-coding",read.csv("../data/z-scores/functional_protein_z_scores.csv", header = TRUE))
sncrna_zscores <- data.frame(Dataset="short-ncrna",read.csv("../data/z-scores/functional_sncrna_z_scores.csv", header = TRUE))
lncrna_zscores <- data.frame(Dataset="lncrna",read.csv("../data/z-scores/functional_lncrna_z_scores.csv", header = TRUE))

prot_neg_zscores <- data.frame(Dataset="protein-negative-control",read.csv("../data/z-scores/protein_negative_control_z_scores.csv", header = TRUE))
sncrna_neg_zscores <- data.frame(Dataset="short-ncrna-negative-control",read.csv("../data/z-scores/sncrna_negative_control_z_scores.csv", header = TRUE))
lncrna_neg_zscores <- data.frame(Dataset="lncrna-negative-control",read.csv("../data/z-scores/lncrna_negative_control_z_scores.csv", header = TRUE))

prot_zscores_all <- rbind(prot_zscores, prot_neg_zscores)
sncrna_zscores_all <- rbind(sncrna_zscores, sncrna_neg_zscores)
lncrna_zscores_all <- rbind(lncrna_zscores, lncrna_neg_zscores)

library(dplyr)
prot_zscores %>% 
  summarise(across(everything(), ~sum(!is.na(.))))

prot_neg_zscores %>% 
  summarise(across(everything(), ~sum(!is.na(.))))

sncrna_zscores %>% 
  summarise(across(everything(), ~sum(!is.na(.))))

sncrna_neg_zscores %>% 
  summarise(across(everything(), ~sum(!is.na(.))))

lncrna_zscores %>% 
  summarise(across(everything(), ~sum(!is.na(.))))

lncrna_neg_zscores %>% 
  summarise(across(everything(), ~sum(!is.na(.))))
