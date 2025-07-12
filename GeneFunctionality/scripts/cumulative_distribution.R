# Load libraries
library(dplyr)
library(ggplot2)
library(Hmisc)
library(ggsignif)
library(scales)  # For percentage formatting

# Load gwas association counts and minimum p-values file
lncrna_gwas_data <- read.csv("../data/merged_coordinates_hg38.csv", header = TRUE)
