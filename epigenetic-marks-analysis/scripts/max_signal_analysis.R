#! /usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Install necessary packages
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

setwd("/Volumes/ADATA HD710 PRO/Downloads/Estefi/Otago University/RsyncCopy/epigenetic-marks-analysis/scripts")

source("./custom_ks_test_signed_d_improved.R")

#################################
# K-S analysis for lncRNA exon1 #
lncrna_data <- read.csv("../data/datasets/histone_feature/H3K79me2/H3K79me2_lncrna-exon1-histone-feature.csv", 
                        header = TRUE, 
                        sep = ",")

lncrna_nc_data <- read.csv("../data/datasets/histone_feature/H3K79me2/H3K79me2_lncrna-exon1-NC-histone-feature.csv", 
                           header = TRUE, 
                           sep = ",")
subset_l <- lncrna_data %>%
  select(AvgSignal,MaxSignal,MaxScaledSignal)
subset_nc_l <- lncrna_nc_data %>%
  select(AvgSignal,MaxSignal,MaxScaledSignal)

# Run the K-S analysis
ks_results_l <- run_ks_tests(subset_nc_l, subset_l)

###############################
# K-S analysis for mRNA exon2 #
mrna_data <- read.csv("../data/datasets/histone_feature/H3K79me2/H3K79me2_protein-exon2-histone-feature.csv", 
                         header = TRUE, 
                         sep = ",")
mrna_nc_data <- read.csv("../data/datasets/histone_feature/H3K79me2/H3K79me2_protein-exon2-NC-histone-feature.csv", 
                            header = TRUE, 
                            sep = ",")

subset_m <- mrna_data %>%
  select(AvgSignal,MaxSignal,MaxScaledSignal)
subset_nc_m <- mrna_nc_data %>%
  select(AvgSignal,MaxSignal,MaxScaledSignal)

# Run the K-S tests
ks_results_m <- run_ks_tests(subset_nc_m, subset_m)

###########################
# K-S analysis for sncRNA #
sncrna_data <- read.csv("../data/datasets/histone_feature/H3K79me2/H3K79me2_short-ncrna-histone-feature.csv", 
                         header = TRUE, 
                         sep = ",")
sncrna_nc_data <- read.csv("../data/datasets/histone_feature/H3K79me2/H3K79me2_short-ncrna-NC-histone-feature.csv", 
                            header = TRUE, 
                            sep = ",")

subset_s <- sncrna_data %>%
  select(AvgSignal,MaxSignal,MaxScaledSignal)
subset_nc_s <- sncrna_nc_data %>%
  select(AvgSignal,MaxSignal,MaxScaledSignal)

# Run the K-S tests
ks_results_s <- run_ks_tests(subset_nc_s, subset_s)
