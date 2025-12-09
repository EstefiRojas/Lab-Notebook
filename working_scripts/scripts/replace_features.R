library(tidyverse)

#####################
# Load raw features #
lncrna_data <- rbind(data.frame(rna_type = "lncrna-exon1",read.csv("data/features/functional-lncrna-exon1-dataset-features.csv")),
                     data.frame(rna_type = "lncrna-exon2",read.csv("data/features/functional-lncrna-exon2-dataset-features.csv")))
sncrna_data <- data.frame(rna_type = "short-ncrna",read.csv("data/features/functional-short-ncrna-dataset-features.csv"))
mrna_data <- rbind(data.frame(rna_type = "protein-exon2",read.csv("data/features/functional-protein-exon2-dataset-features.csv")),
                   data.frame(rna_type = "protein-exon3",read.csv("data/features/functional-protein-exon3-dataset-features.csv")))

lncrna_data_nc <- rbind(data.frame(rna_type = "lncrna-exon1",read.csv("data/features/lncrna-exon1-negative-control-dataset-features.csv")),
                        data.frame(rna_type = "lncrna-exon2",read.csv("data/features/lncrna-exon2-negative-control-dataset-features.csv")))
sncrna_data_nc <- data.frame(rna_type = "short-ncrna",read.csv("data/features/short-ncrna-negative-control-dataset-features.csv"))
mrna_data_nc <- rbind(data.frame(rna_type = "protein-exon2",read.csv("data/features/protein-exon2-negative-control-dataset-features.csv")),
                      data.frame(rna_type = "protein-exon3",read.csv("data/features/protein-exon3-negative-control-dataset-features.csv")))


# Get new features computed on 30 Nov 2025 after code fix (coding potential, max covariance)
# Coding Potential
lncrna_coding_potential <- read.csv("data/features/20251130/robust_z_scores_positive_lncrna_coding_potential_20251130.csv")
sncrna_coding_potential <- read.csv("data/features/20251130/robust_z_scores_positive_sncrna_coding_potential_20251130.csv")
mrna_coding_potential <- read.csv("data/features/20251130/robust_z_scores_positive_mrna_coding_potential_20251130.csv")
lncrna_coding_potential_nc <- read.csv("data/features/20251130/robust_z_scores_negative_lncrna_coding_potential_20251130.csv")
sncrna_coding_potential_nc <- read.csv("data/features/20251130/robust_z_scores_negative_sncrna_coding_potential_20251130.csv")
mrna_coding_potential_nc <- read.csv("data/features/20251130/robust_z_scores_negative_mrna_coding_potential_20251130.csv")

# Max covariance
lncrna_max_cov <- read.csv("data/features/20251130/robust_z_scores_positive_lncrna_Max_covariance_20251130.csv")
sncrna_max_cov <- read.csv("data/features/20251130/robust_z_scores_positive_sncrna_Max_covariance_20251130.csv")
mrna_max_cov <- read.csv("data/features/20251130/robust_z_scores_positive_mrna_Max_covariance_20251130.csv")
lncrna_max_cov_nc <- read.csv("data/features/20251130/robust_z_scores_negative_lncrna_Max_covariance_20251130.csv")
sncrna_max_cov_nc <- read.csv("data/features/20251130/robust_z_scores_negative_sncrna_Max_covariance_20251130.csv")
mrna_max_cov_nc <- read.csv("data/features/20251130/robust_z_scores_negative_mrna_Max_covariance_20251130.csv")

###########################
# replace with new values #
# Positive cases
lncrna_data_new <- lncrna_data %>%
  # Join coding potential
  left_join(lncrna_coding_potential, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z")) %>%
  # Join Max covariance
  left_join(lncrna_max_cov, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z"))

sncrna_data_new <- sncrna_data %>%
  # Join coding potential
  left_join(sncrna_coding_potential, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z")) %>%
  # Join Max covariance
  left_join(sncrna_max_cov, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z"))

mrna_data_new <- mrna_data %>%
  # Join coding potential
  left_join(mrna_coding_potential, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z")) %>%
  # Join Max covariance
  left_join(mrna_max_cov, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z"))

# Negative cases
lncrna_data_new_nc <- lncrna_data_nc %>%
  # Join coding potential
  left_join(lncrna_coding_potential_nc, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z")) %>%
  # Join Max covariance
  left_join(lncrna_max_cov_nc, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z"))

sncrna_data_new_nc <- sncrna_data_nc %>%
  # Join coding potential
  left_join(sncrna_coding_potential_nc, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z")) %>%
  # Join Max covariance
  left_join(sncrna_max_cov_nc, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z"))

mrna_data_new_nc <- mrna_data_nc %>%
  # Join coding potential
  left_join(mrna_coding_potential_nc, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z")) %>%
  # Join Max covariance
  left_join(mrna_max_cov_nc, by = c("ID", "rna_type"), suffix = c("_OLD", "")) %>%
  # Remove the old column 
  select(-ends_with("_OLD"),-starts_with("robust_z"))

#################################
# Write new data into csv files #
write.csv(lncrna_data_new %>% filter(rna_type=="lncrna-exon1") %>% select(-rna_type), "data/features/functional-lncrna-exon1-dataset-features-30-11-2025.csv", row.names = FALSE)
write.csv(lncrna_data_new %>% filter(rna_type=="lncrna-exon2") %>% select(-rna_type), "data/features/functional-lncrna-exon2-dataset-features-30-11-2025.csv", row.names = FALSE)
write.csv(sncrna_data_new %>% filter(rna_type=="short-ncrna") %>% select(-rna_type), "data/features/functional-short-ncrna-dataset-features-30-11-2025.csv", row.names = FALSE)
write.csv(mrna_data_new %>% filter(rna_type=="protein-exon2") %>% select(-rna_type), "data/features/functional-protein-exon2-dataset-features-30-11-2025.csv", row.names = FALSE)
write.csv(mrna_data_new %>% filter(rna_type=="protein-exon3") %>% select(-rna_type), "data/features/functional-protein-exon3-dataset-features-30-11-2025.csv", row.names = FALSE)

write.csv(lncrna_data_new_nc %>% filter(rna_type=="lncrna-exon1") %>% select(-rna_type), "data/features/lncrna-exon1-negative-control-dataset-features-30-11-2025.csv", row.names = FALSE)
write.csv(lncrna_data_new_nc %>% filter(rna_type=="lncrna-exon2") %>% select(-rna_type), "data/features/lncrna-exon2-negative-control-dataset-features-30-11-2025.csv", row.names = FALSE)
write.csv(sncrna_data_new_nc %>% filter(rna_type=="short-ncrna") %>% select(-rna_type), "data/features/short-ncrna-negative-control-dataset-features-30-11-2025.csv", row.names = FALSE)
write.csv(mrna_data_new_nc %>% filter(rna_type=="protein-exon2") %>% select(-rna_type), "data/features/protein-exon2-negative-control-dataset-features-30-11-2025.csv", row.names = FALSE)
write.csv(mrna_data_new_nc %>% filter(rna_type=="protein-exon3") %>% select(-rna_type), "data/features/protein-exon3-negative-control-dataset-features-30-11-2025.csv", row.names = FALSE)

