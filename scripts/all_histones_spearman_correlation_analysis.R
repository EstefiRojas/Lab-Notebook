# Load necessary library
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
library(boot) # Used to obtain confidence intervals
install.packages("psych")
library(psych)

help("corr.test")
###################################################################
# Sample data vectors for positive and negative controls LONG NCRNA
lncrna_file <- "../data/histone_feature/H3k36me3/H3k36me3_lncrna_matrix.csv"
data_lncrna <- read.csv(lncrna_file, header=TRUE, sep = "\t")

# Calculate basic statistics
summary(data_lncrna$AvgSignal)
# Statistics split by 'Functional' category
by(lncrna_file1$AvgSignal, data_lncrna$Functional, summary) 

# Split data by 'Functional' for t-test
negative_control <- data_lncrna[data_lncrna$Functional == "No", "AvgSignal"]
positive_data <- data_lncrna[data_lncrna$Functional == "Yes", "AvgSignal"]

# Kolmogorov-Smirnov Test (Normality)
ks.test(negative_control, "pnorm", mean(negative_control), sd(negative_control))
ks.test(positive_data, "pnorm", mean(positive_data), sd(positive_data))

# Student's t-test (assuming equal variances)
#t.test(negative_control, positive_data, var.equal = TRUE)

# Mann-Whitney U Test (aka Wilcoxon Rank-Sum Test)
wilcox.test(negative_control, positive_data)
############################################



####################################################################
# Sample data vectors for positive and negative controls SHORT NCRNA
sncrna_file <- "../data/histone_feature/H3k36me3/H3k36me3_short_ncrna_matrix.csv"
data_sncrna <- read.csv(sncrna_file, header=TRUE, sep="\t")

# Calculate basic statistics
summary(data_sncrna$AvgSignal)
# Statistics split by 'Functional' category
by(data_sncrna$AvgSignal, data_sncrna$Functional, summary) 

# Split data by 'Functional' for t-test
negative_control <- data_sncrna[data_sncrna$Functional == "No", "AvgSignal"]
positive_data <- data_sncrna[data_sncrna$Functional == "Yes", "AvgSignal"]

# Kolmogorov-Smirnov Test (Normality)
ks.test(negative_control, "pnorm", mean(negative_control), sd(negative_control))
ks.test(positive_data, "pnorm", mean(positive_data), sd(positive_data))

# Student's t-test (assuming equal variances)
#t.test(negative_control, positive_data, var.equal = TRUE)

# Mann-Whitney U Test (aka Wilcoxon Rank-Sum Test)
wilcox.test(negative_control, positive_data)
############################################


#######################################################################
# Sample data vectors for positive and negative controls PROTEIN CODING
protein_file <- "../data/histone_feature/H3k36me3/H3k36me3_protein_matrix.csv"
data_prot <- read.csv(protein_file, header=TRUE, sep="\t")

# Calculate basic statistics
summary(data_prot$AvgSignal)
# Statistics split by 'Functional' category
by(data_prot$AvgSignal, data_prot$Functional, summary) 

# Split data by 'Functional' for t-test
negative_control <- data_prot[data_prot$Functional == "No", "AvgSignal"]
positive_data <- data_prot[data_prot$Functional == "Yes", "AvgSignal"]

# Kolmogorov-Smirnov Test (Normality)
ks.test(negative_control, "pnorm", mean(negative_control), sd(negative_control))
ks.test(positive_data, "pnorm", mean(positive_data), sd(positive_data))

# Student's t-test (assuming equal variances)
#t.test(negative_control, positive_data, var.equal = TRUE)

# Mann-Whitney U Test (aka Wilcoxon Rank-Sum Test)
wilcox.test(negative_control, positive_data)
############################################


##########################
# Fisher's Test for LNCRNA
contingency_table <- table(data_lncrna$Functional, data_lncrna$HistonePresence)
contingency_table

fisher.test(contingency_table)
##############################

###############################
# Fisher's Test for SHORT NCRNA
contingency_table <- table(data_sncrna$Functional, data_sncrna$HistonePresence)
contingency_table

fisher.test(contingency_table)
##############################

##################################
# Fisher's Test for PROTEIN CODING
contingency_table <- table(data_prot$Functional, data_prot$HistonePresence)
contingency_table

fisher.test(contingency_table)
##############################

####################################################################
# Define a function to calculate Spearman correlation from a dataset
calc_spearman <- function(data, indices) {
  d <- data[indices,] # Sample data based on indices
  cor.test(d[,1], d[,2], method = "spearman")$estimate
}

# Set random seed for reproducibility
set.seed(123) 

# Bootstrap with, 1000 replications
boot_prot_ci_results <- boot(data = data.frame(data_prot[,4], data_prot[,5]), 
                     statistic = calc_spearman, 
                     R = 1000)
boot_lncrna_ci_results <- boot(data = data.frame(data_lncrna[,4], data_lncrna[,5]), 
                             statistic = calc_spearman, 
                             R = 1000)
boot_shortncrna_ci_results <- boot(data = data.frame(data_sncrna[,4], data_sncrna[,5]), 
                               statistic = calc_spearman, 
                               R = 1000)
# Get confidence intervals
boot.ci(boot_prot_ci_results, type = "perc") # Use percentile-based confidence intervals 
boot.ci(boot_lncrna_ci_results, type = "perc") # Use percentile-based confidence intervals 
boot.ci(boot_shortncrna_ci_results, type = "perc") # Use percentile-based confidence intervals 
##################################################



#args <- commandArgs(trailingOnly = TRUE)
#file <- args[1]


options(repos = "https://cloud.r-project.org")
install.packages("RVAideMemoire")
library(RVAideMemoire)
install.packages("randomForest")
library(randomForest)

# 1. Read the data from the file into a data frame
spearman_correlation <- function(file){
  
  data <- read.csv(file, header=TRUE, sep = "\t")
  
  # Convert Functional "Yes" "No" into numeric values 
  data$Functional <- factor(data$Functional, levels=c("No","Yes")) 
  functional_numeric <- as.numeric(factor(data$Functional)) - 1 
  functional_numeric
  # Approx NA values (CHECK method) 
  #data <- na.roughfix(data)
  
  # 2. Extract the other numeric columns
  numeric_data <- data[, !names(data) %in% c("Chromosome","ID","Functional")]
  
  # Calculate the Spearman correlation between "functional" column and each numeric column
  correlations <- sapply(numeric_data, function(col) {
    spearman <- spearman.ci(functional_numeric, col, conf.level = 0.95, nrep = 1000)
    cor_matrix <- cor.test(functional_numeric, col, method = "spearman", alpha = 0.05)
    rho <- spearman$estimate[[1]]
    ci_inf <- spearman$conf.int[[1]] 
    ci_sup <- spearman$conf.int[[2]]
    pval <- cor_matrix$p.value
    return(c( cor = rho, ci_inf = ci_inf, ci_sup = ci_sup, pval = pval))
  })
  
  # Return correlation coefficients and confidence intervals as a named list with variable name as the name
  return(correlations)
}

load_and_run <- function(histone_name, correlation_list){
  #H2AFZ - cell line
  lncrna_file <- paste0("../data/histone_feature/",histone_name,"/",histone_name,"_lncrna_matrix.csv")
  #output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"
  
  sncrna_file <- paste0("../data/histone_feature/",histone_name,"/",histone_name,"_short_ncrna_matrix.csv")
  #output_file <- "../data/histone_feature/protein_spearman_matrix.csv"
  
  protein_file <- paste0("../data/histone_feature/",histone_name,"/",histone_name,"_protein_matrix.csv")
  #output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"
  
  correlation_list[[lncrna_file]] <- spearman_correlation(lncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[sncrna_file]] <- spearman_correlation(sncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[protein_file]] <- spearman_correlation(protein_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  return(correlation_list)
}

load_and_run_by_biosample <- function(histone_name, biosample_class, correlation_list){
  #H2AFZ - cell line
  lncrna_file <- paste0("../data/histone_feature/",histone_name,"_",biosample_class,"/",histone_name,"_",biosample_class,"_lncrna_matrix.csv")
  #output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"
  
  sncrna_file <- paste0("../data/histone_feature/",histone_name,"_",biosample_class,"/",histone_name,"_",biosample_class,"_short_ncrna_matrix.csv")
  #output_file <- "../data/histone_feature/protein_spearman_matrix.csv"
  
  protein_file <- paste0("../data/histone_feature/",histone_name,"_",biosample_class,"/",histone_name,"_",biosample_class,"_protein_matrix.csv")
  #output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"
  
  correlation_list[[lncrna_file]] <- spearman_correlation(lncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[sncrna_file]] <- spearman_correlation(sncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[protein_file]] <- spearman_correlation(protein_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  return(correlation_list)
}

load_and_run_chrm_by_biosample <- function(chrm_name, correlation_list){
  #H2AFZ - cell line
  lncrna_file <- paste0("../data/chrm_acc_feature/by_biosample_class/",chrm_name,"/chrm_acc_lncrna_matrix.csv")
  #output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"
  
  sncrna_file <- paste0("../data/chrm_acc_feature/by_biosample_class/",chrm_name,"/chrm_acc_short_ncrna_matrix.csv")
  #output_file <- "../data/histone_feature/protein_spearman_matrix.csv"
  
  protein_file <- paste0("../data/chrm_acc_feature/by_biosample_class/",chrm_name,"/chrm_acc_protein_matrix.csv")
  #output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"
  
  correlation_list[[lncrna_file]] <- spearman_correlation(lncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[sncrna_file]] <- spearman_correlation(sncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[protein_file]] <- spearman_correlation(protein_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  return(correlation_list)
}

load_and_run_chrm <- function(correlation_list){
  #H2AFZ - cell line
  lncrna_file <- paste0("../data/chrm_acc_feature/by_biosample_class/chrm_acc_lncrna_matrix.csv")
  #output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"
  
  sncrna_file <- paste0("../data/chrm_acc_feature/by_biosample_class/chrm_acc_short_ncrna_matrix.csv")
  #output_file <- "../data/histone_feature/protein_spearman_matrix.csv"
  
  protein_file <- paste0("../data/chrm_acc_feature/by_biosample_class/chrm_acc_protein_matrix.csv")
  #output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"
  
  correlation_list[[lncrna_file]] <- spearman_correlation(lncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[sncrna_file]] <- spearman_correlation(sncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[protein_file]] <- spearman_correlation(protein_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  return(correlation_list)
}


load_and_run_methyl_by_biosample <- function(methyl_name, correlation_list){
  #H2AFZ - cell line
  lncrna_file <- paste0("../data/methylome_feature/by_biosample_class/",methyl_name,"/lncrna_matrix.csv")
  #output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"
  
  sncrna_file <- paste0("../data/methylome_feature/by_biosample_class/",methyl_name,"/short_ncrna_matrix.csv")
  #output_file <- "../data/histone_feature/protein_spearman_matrix.csv"
  
  protein_file <- paste0("../data/methylome_feature/by_biosample_class/",methyl_name,"/protein_matrix.csv")
  #output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"
  
  correlation_list[[lncrna_file]] <- spearman_correlation(lncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[sncrna_file]] <- spearman_correlation(sncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[protein_file]] <- spearman_correlation(protein_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  return(correlation_list)
}

load_and_run_methyl <- function(correlation_list){
  #H2AFZ - cell line
  lncrna_file <- paste0("../data/methylome_feature/lncrna_matrix.csv")
  #output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"
  
  sncrna_file <- paste0("../data/methylome_feature/short_ncrna_matrix.csv")
  #output_file <- "../data/histone_feature/protein_spearman_matrix.csv"
  
  protein_file <- paste0("../data/methylome_feature/protein_matrix.csv")
  #output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"
  
  correlation_list[[lncrna_file]] <- spearman_correlation(lncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[sncrna_file]] <- spearman_correlation(sncrna_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  correlation_list[[protein_file]] <- spearman_correlation(protein_file)
  #correlation_matrix_combined <- do.call(rbind, correlation_list)
  
  return(correlation_list)
}


list_of_histones <- list("H2AFZ",      "H2BK15ac",  "H3K18ac",   "H3K27me3",  "H3K56ac",   "H3K9me2",  "H4K5ac",
                         "H2AK5ac",    "H2BK20ac",  "H3K20me1",  "H3K36me3",  "H3K79me1",  "H3K9me3",  "H4K8ac",
                         "H2AK9ac",    "H2BK5ac",   "H3K23ac",   "H3K4ac",    "H3K79me2",  "H3T11ph",  "H4K91ac",
                         "H2BK120ac",  "H3F3A",     "H3K23me2",  "H3K4me1",   "H3K9ac",    "H4K12ac",
                         "H2BK12ac",   "H3K14ac",   "H3K27ac",   "H3K4me3",   "H3K9me1",   "H4K4me2")

list_of_chrm_acc_experiments <-list("ATAC-seq_cell-line","ATAC-seq_IVDC","ATAC-seq_primary-cell","ATAC-seq_tissue",
                                    "DNAse-seq_cell-line","DNAse-seq_IVDC","DNAse-seq_primary-cell","DNAse-seq_tissue")

list_of_methylome_experiments <-list("WGBS_cell-line","WGBS_IVDC","WGBS_primary-cell","WGBS_tissue")

# Initialize empty list to store correlation vectors
correlation_list <- list()



correlation_list2 <- list()
correlation_list2 <- load_and_run("H2AFZ", correlation_list2)
for (hn in list_of_histones) {
  correlation_list <-load_and_run(hn, correlation_list)
}
corr_list_hist <- correlation_list


for (hn in list_of_chrm_acc_experiments) {
  correlation_list <-load_and_run_chrm_by_biosample(hn, correlation_list)
}
corr_list_chrmacc_by_biosample <- correlation_list

correlation_list <- load_and_run_chrm(correlation_list)

for (hn in list_of_methylome_experiments) {
  correlation_list <-load_and_run_methyl_by_biosample(hn, correlation_list)
}
coor_list_methylome_by_biosample <- correlation_list

correlation_list <- load_and_run_methyl(correlation_list)


correlation_list_all_hist_chrm_methyl <- correlation_list

#ChrmAcc - cell line
lncrna_file <- paste0("../data/chrm_acc_feature/by_biosample_class/chrm_acc_lncrna_matrix.csv")
#output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"

sncrna_file <- paste0("../data/chrm_acc_feature/by_biosample_class/chrm_acc_short_ncrna_matrix.csv")
#output_file <- "../data/histone_feature/protein_spearman_matrix.csv"

protein_file <- paste0("../data/chrm_acc_feature/by_biosample_class/chrm_acc_protein_matrix.csv")
#output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"

correlation_list[[lncrna_file]] <- spearman_correlation(lncrna_file)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[sncrna_file]] <- spearman_correlation(sncrna_file)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[protein_file]] <- spearman_correlation(protein_file)


correlation_list <- load_and_run("H3K79me2", correlation_list)



correlation_list <-load_and_run("H2AFZ","cell-line", correlation_list)
correlation_list <-load_and_run("H2AFZ","IVDC", correlation_list)
correlation_list <-load_and_run("H2AFZ","primary-cell", correlation_list)
correlation_list <-load_and_run("H2AFZ","tissue", correlation_list)

correlation_list <-load_and_run("H2AK5ac","cell-line", correlation_list)
correlation_list <-load_and_run("H2AK5ac","IVDC", correlation_list)

correlation_list <-load_and_run("H2AK9ac","cell-line", correlation_list)

correlation_list <-load_and_run("H2BK5ac","cell-line", correlation_list)
correlation_list <-load_and_run("H2BK5ac","IVDC", correlation_list)

correlation_list <-load_and_run("H2BK12ac","cell-line", correlation_list)
correlation_list <-load_and_run("H2BK12ac","IVDC", correlation_list)

correlation_list <-load_and_run("H2BK15ac","cell-line", correlation_list)
correlation_list <-load_and_run("H2BK15ac","IVDC", correlation_list)

correlation_list <-load_and_run("H2BK20ac","cell-line", correlation_list)

correlation_list <-load_and_run("H2BK120ac","cell-line", correlation_list)
correlation_list <-load_and_run("H2BK120ac","IVDC", correlation_list)

correlation_list <-load_and_run("H3F3A","cell-line", correlation_list)
correlation_list <-load_and_run("H3F3A","IVDC", correlation_list)

correlation_list <-load_and_run("H3K4ac","cell-line", correlation_list)
correlation_list <-load_and_run("H3K4ac","IVDC", correlation_list)

correlation_list <-load_and_run("H3K4me1","cell-line", correlation_list)
correlation_list <-load_and_run("H3K4me1","IVDC", correlation_list)
correlation_list <-load_and_run("H3K4me1","primary-cell", correlation_list)
correlation_list <-load_and_run("H3K4me1","tissue", correlation_list)

correlation_list <-load_and_run("H3K4me3","cell-line", correlation_list)
correlation_list <-load_and_run("H3K4me3","IVDC", correlation_list)
correlation_list <-load_and_run("H3K4me3","primary-cell", correlation_list)
correlation_list <-load_and_run("H3K4me3","tissue", correlation_list)

correlation_list <-load_and_run("H3K9ac","cell-line", correlation_list)
correlation_list <-load_and_run("H3K9ac","IVDC", correlation_list)
correlation_list <-load_and_run("H3K9ac","primary-cell", correlation_list)
correlation_list <-load_and_run("H3K9ac","tissue", correlation_list)

correlation_list <-load_and_run("H3K9me1","cell-line", correlation_list)
correlation_list <-load_and_run("H3K9me1","primary-cell", correlation_list)

correlation_list <-load_and_run("H3K9me2","cell-line", correlation_list)
correlation_list <-load_and_run("H3K9me2","IVDC", correlation_list)

correlation_list <-load_and_run("H3K9me3","cell-line", correlation_list)
correlation_list <-load_and_run("H3K9me3","IVDC", correlation_list)
correlation_list <-load_and_run("H3K9me3","primary-cell", correlation_list)
correlation_list <-load_and_run("H3K9me3","tissue", correlation_list)

correlation_list <-load_and_run("H3K14ac","cell-line", correlation_list)
correlation_list <-load_and_run("H3K14ac","IVDC", correlation_list)

correlation_list <-load_and_run("H3K18ac","cell-line", correlation_list)
correlation_list <-load_and_run("H3K18ac","IVDC", correlation_list)

correlation_list <-load_and_run("H3K20me1","cell-line", correlation_list)
correlation_list <-load_and_run("H3K20me1","IVDC", correlation_list)
correlation_list <-load_and_run("H3K20me1","primary-cell", correlation_list)

correlation_list <-load_and_run("H3K23ac","cell-line", correlation_list)
correlation_list <-load_and_run("H3K23ac","IVDC", correlation_list)

correlation_list <-load_and_run("H3K23me2","cell-line", correlation_list)

correlation_list <-load_and_run("H3K27ac","cell-line", correlation_list)
correlation_list <-load_and_run("H3K27ac","IVDC", correlation_list)
correlation_list <-load_and_run("H3K27ac","primary-cell", correlation_list)
correlation_list <-load_and_run("H3K27ac","tissue", correlation_list)

correlation_list <-load_and_run("H3K27me3","cell-line", correlation_list)
correlation_list <-load_and_run("H3K27me3","IVDC", correlation_list)
correlation_list <-load_and_run("H3K27me3","primary-cell", correlation_list)
correlation_list <-load_and_run("H3K27me3","tissue", correlation_list)

correlation_list <-load_and_run("H3K36me3","cell-line", correlation_list)
correlation_list <-load_and_run("H3K36me3","IVDC", correlation_list)
correlation_list <-load_and_run("H3K36me3","primary-cell", correlation_list)
correlation_list <-load_and_run("H3K36me3","tissue", correlation_list)

correlation_list <-load_and_run("H3K56ac","cell-line", correlation_list)
correlation_list <-load_and_run("H3K56ac","IVDC", correlation_list)

correlation_list <-load_and_run("H3K79me1","cell-line", correlation_list)
correlation_list <-load_and_run("H3K79me1","IVDC", correlation_list)

correlation_list <-load_and_run("H3K79me2","cell-line", correlation_list)
correlation_list <-load_and_run("H3K79me2","IVDC", correlation_list)
correlation_list <-load_and_run("H3K79me2","primary-cell", correlation_list)
correlation_list <-load_and_run("H3K79me2","tissue", correlation_list)

correlation_list <-load_and_run("H3T11ph","cell-line", correlation_list)

correlation_list <-load_and_run("H4K4me2","cell-line", correlation_list)
correlation_list <-load_and_run("H4K4me2","IVDC", correlation_list)
correlation_list <-load_and_run("H4K4me2","primary-cell", correlation_list)
correlation_list <-load_and_run("H4K4me2","tissue", correlation_list)

correlation_list <-load_and_run("H4K5ac","cell-line", correlation_list)

correlation_list <-load_and_run("H4K8ac","cell-line", correlation_list)
correlation_list <-load_and_run("H4K8ac","IVDC", correlation_list)

correlation_list <-load_and_run("H4K12ac","IVDC", correlation_list)

correlation_list <-load_and_run("H4K91ac","cell-line", correlation_list)
correlation_list <-load_and_run("H4K91ac","IVDC", correlation_list)

corr_list_hist_by_biosample <- correlation_list



correlation_list[[3]]
######################################
# Create a graph to visualice features
df <- data.frame(sequence_type = c("lncRNA Exons","Short ncRNAs","Protein-coding Exons"),
                 correlation = c(correlation_list[[1]][1,3], correlation_list[[2]][1,3], correlation_list[[3]][1,3],
                                 correlation_list[[4]][1,3], correlation_list[[5]][1,3], correlation_list[[6]][1,3],
                                 correlation_list[[7]][1,3], correlation_list[[8]][1,3], correlation_list[[9]][1,3],
                                 correlation_list[[10]][1,3], correlation_list[[11]][1,3], correlation_list[[12]][1,3],
                                 correlation_list[[13]][1,3], correlation_list[[14]][1,3], correlation_list[[15]][1,3],
                                 correlation_list[[16]][1,3], correlation_list[[17]][1,3], correlation_list[[18]][1,3],
                                 correlation_list[[19]][1,3], correlation_list[[20]][1,3], correlation_list[[21]][1,3],
                                 correlation_list[[22]][1,3], correlation_list[[23]][1,3], correlation_list[[24]][1,3],
                                 correlation_list[[25]][1,3], correlation_list[[26]][1,3], correlation_list[[27]][1,3],
                                 correlation_list[[28]][1,3], correlation_list[[29]][1,3], correlation_list[[30]][1,3],
                                 correlation_list[[31]][1,3], correlation_list[[32]][1,3], correlation_list[[33]][1,3],
                                 correlation_list[[34]][1,3], correlation_list[[35]][1,3], correlation_list[[36]][1,3],
                                 correlation_list[[37]][1,3], correlation_list[[38]][1,3], correlation_list[[39]][1,3],
                                 correlation_list[[40]][1,3], correlation_list[[41]][1,3], correlation_list[[42]][1,3],
                                 correlation_list[[43]][1,3], correlation_list[[44]][1,3], correlation_list[[45]][1,3],
                                 correlation_list[[46]][1,3], correlation_list[[47]][1,3], correlation_list[[48]][1,3],
                                 correlation_list[[49]][1,3], correlation_list[[50]][1,3], correlation_list[[51]][1,3],
                                 correlation_list[[52]][1,3], correlation_list[[53]][1,3], correlation_list[[54]][1,3],
                                 correlation_list[[55]][1,3], correlation_list[[56]][1,3], correlation_list[[57]][1,3],
                                 correlation_list[[58]][1,3], correlation_list[[59]][1,3], correlation_list[[60]][1,3],
                                 correlation_list[[61]][1,3], correlation_list[[62]][1,3], correlation_list[[63]][1,3],
                                 correlation_list[[64]][1,3], correlation_list[[65]][1,3], correlation_list[[66]][1,3],
                                 correlation_list[[67]][1,3], correlation_list[[68]][1,3], correlation_list[[69]][1,3],
                                 correlation_list[[70]][1,3], correlation_list[[71]][1,3], correlation_list[[72]][1,3],
                                 correlation_list[[73]][1,3], correlation_list[[74]][1,3], correlation_list[[75]][1,3],
                                 correlation_list[[76]][1,3], correlation_list[[77]][1,3], correlation_list[[78]][1,3],
                                 correlation_list[[79]][1,3], correlation_list[[80]][1,3], correlation_list[[81]][1,3],
                                 correlation_list[[82]][1,3], correlation_list[[83]][1,3], correlation_list[[84]][1,3],
                                 correlation_list[[85]][1,3], correlation_list[[86]][1,3], correlation_list[[87]][1,3],
                                 correlation_list[[88]][1,3], correlation_list[[89]][1,3], correlation_list[[90]][1,3],
                                 correlation_list[[91]][1,3], correlation_list[[92]][1,3], correlation_list[[93]][1,3],
                                 correlation_list[[94]][1,3], correlation_list[[95]][1,3], correlation_list[[96]][1,3],
                                 correlation_list[[97]][1,3], correlation_list[[98]][1,3], correlation_list[[99]][1,3],
                                 correlation_list[[100]][1,3], correlation_list[[101]][1,3], correlation_list[[102]][1,3],
                                 correlation_list[[103]][1,3], correlation_list[[104]][1,3], correlation_list[[105]][1,3]
                 ),
                 lower_ci = c(correlation_list[[1]][2,3], correlation_list[[2]][2,3], correlation_list[[3]][2,3],
                              correlation_list[[4]][2,3], correlation_list[[5]][2,3], correlation_list[[6]][2,3],
                              correlation_list[[7]][2,3], correlation_list[[8]][2,3], correlation_list[[9]][2,3],
                              correlation_list[[10]][2,3], correlation_list[[11]][2,3], correlation_list[[12]][2,3],
                              correlation_list[[13]][2,3], correlation_list[[14]][2,3], correlation_list[[15]][2,3],
                              correlation_list[[16]][2,3], correlation_list[[17]][2,3], correlation_list[[18]][2,3],
                              correlation_list[[19]][2,3], correlation_list[[20]][2,3], correlation_list[[21]][2,3],
                              correlation_list[[22]][2,3], correlation_list[[23]][2,3], correlation_list[[24]][2,3],
                              correlation_list[[25]][2,3], correlation_list[[26]][2,3], correlation_list[[27]][2,3],
                              correlation_list[[28]][2,3], correlation_list[[29]][2,3], correlation_list[[30]][2,3],
                              correlation_list[[31]][2,3], correlation_list[[32]][2,3], correlation_list[[33]][2,3],
                              correlation_list[[34]][2,3], correlation_list[[35]][2,3], correlation_list[[36]][2,3],
                              correlation_list[[37]][2,3], correlation_list[[38]][2,3], correlation_list[[39]][2,3],
                              correlation_list[[40]][2,3], correlation_list[[41]][2,3], correlation_list[[42]][2,3],
                              correlation_list[[43]][2,3], correlation_list[[44]][2,3], correlation_list[[45]][2,3],
                              correlation_list[[46]][2,3], correlation_list[[47]][2,3], correlation_list[[48]][2,3],
                              correlation_list[[49]][2,3], correlation_list[[50]][2,3], correlation_list[[51]][2,3],
                              correlation_list[[52]][2,3], correlation_list[[53]][2,3], correlation_list[[54]][2,3],
                              correlation_list[[55]][2,3], correlation_list[[56]][2,3], correlation_list[[57]][2,3],
                              correlation_list[[58]][2,3], correlation_list[[59]][2,3], correlation_list[[60]][2,3],
                              correlation_list[[61]][2,3], correlation_list[[62]][2,3], correlation_list[[63]][2,3],
                              correlation_list[[64]][2,3], correlation_list[[65]][2,3], correlation_list[[66]][2,3],
                              correlation_list[[67]][2,3], correlation_list[[68]][2,3], correlation_list[[69]][2,3],
                              correlation_list[[70]][2,3], correlation_list[[71]][2,3], correlation_list[[72]][2,3],
                              correlation_list[[73]][2,3], correlation_list[[74]][2,3], correlation_list[[75]][2,3],
                              correlation_list[[76]][2,3], correlation_list[[77]][2,3], correlation_list[[78]][2,3],
                              correlation_list[[79]][2,3], correlation_list[[80]][2,3], correlation_list[[81]][2,3],
                              correlation_list[[82]][2,3], correlation_list[[83]][2,3], correlation_list[[84]][2,3],
                              correlation_list[[85]][2,3], correlation_list[[86]][2,3], correlation_list[[87]][2,3],
                              correlation_list[[88]][2,3], correlation_list[[89]][2,3], correlation_list[[90]][2,3],
                              correlation_list[[91]][2,3], correlation_list[[92]][2,3], correlation_list[[93]][2,3],
                              correlation_list[[94]][2,3], correlation_list[[95]][2,3], correlation_list[[96]][2,3],
                              correlation_list[[97]][2,3], correlation_list[[98]][2,3], correlation_list[[99]][2,3],
                              correlation_list[[100]][2,3], correlation_list[[101]][2,3], correlation_list[[102]][2,3],
                              correlation_list[[103]][2,3], correlation_list[[104]][2,3], correlation_list[[105]][2,3]
                 ),
                 upper_ci = c(correlation_list[[1]][3,3], correlation_list[[2]][3,3], correlation_list[[3]][3,3],
                              correlation_list[[4]][3,3], correlation_list[[5]][3,3], correlation_list[[6]][3,3],
                              correlation_list[[7]][3,3], correlation_list[[8]][3,3], correlation_list[[9]][3,3],
                              correlation_list[[10]][3,3], correlation_list[[11]][3,3], correlation_list[[12]][3,3],
                              correlation_list[[13]][3,3], correlation_list[[14]][3,3], correlation_list[[15]][3,3],
                              correlation_list[[16]][3,3], correlation_list[[17]][3,3], correlation_list[[18]][3,3],
                              correlation_list[[19]][3,3], correlation_list[[20]][3,3], correlation_list[[21]][3,3],
                              correlation_list[[22]][3,3], correlation_list[[23]][3,3], correlation_list[[24]][3,3],
                              correlation_list[[25]][3,3], correlation_list[[26]][3,3], correlation_list[[27]][3,3],
                              correlation_list[[28]][3,3], correlation_list[[29]][3,3], correlation_list[[30]][3,3],
                              correlation_list[[31]][3,3], correlation_list[[32]][3,3], correlation_list[[33]][3,3],
                              correlation_list[[34]][3,3], correlation_list[[35]][3,3], correlation_list[[36]][3,3],
                              correlation_list[[37]][3,3], correlation_list[[38]][3,3], correlation_list[[39]][3,3],
                              correlation_list[[40]][3,3], correlation_list[[41]][3,3], correlation_list[[42]][3,3],
                              correlation_list[[43]][3,3], correlation_list[[44]][3,3], correlation_list[[45]][3,3],
                              correlation_list[[46]][3,3], correlation_list[[47]][3,3], correlation_list[[48]][3,3],
                              correlation_list[[49]][3,3], correlation_list[[50]][3,3], correlation_list[[51]][3,3],
                              correlation_list[[52]][3,3], correlation_list[[53]][3,3], correlation_list[[54]][3,3],
                              correlation_list[[55]][3,3], correlation_list[[56]][3,3], correlation_list[[57]][3,3],
                              correlation_list[[58]][3,3], correlation_list[[59]][3,3], correlation_list[[60]][3,3],
                              correlation_list[[61]][3,3], correlation_list[[62]][3,3], correlation_list[[63]][3,3],
                              correlation_list[[64]][3,3], correlation_list[[65]][3,3], correlation_list[[66]][3,3],
                              correlation_list[[67]][3,3], correlation_list[[68]][3,3], correlation_list[[69]][3,3],
                              correlation_list[[70]][3,3], correlation_list[[71]][3,3], correlation_list[[72]][3,3],
                              correlation_list[[73]][3,3], correlation_list[[74]][3,3], correlation_list[[75]][3,3],
                              correlation_list[[76]][3,3], correlation_list[[77]][3,3], correlation_list[[78]][3,3],
                              correlation_list[[79]][3,3], correlation_list[[80]][3,3], correlation_list[[81]][3,3],
                              correlation_list[[82]][3,3], correlation_list[[83]][3,3], correlation_list[[84]][3,3],
                              correlation_list[[85]][3,3], correlation_list[[86]][3,3], correlation_list[[87]][3,3],
                              correlation_list[[88]][3,3], correlation_list[[89]][3,3], correlation_list[[90]][3,3],
                              correlation_list[[91]][3,3], correlation_list[[92]][3,3], correlation_list[[93]][3,3],
                              correlation_list[[94]][3,3], correlation_list[[95]][3,3], correlation_list[[96]][3,3],
                              correlation_list[[97]][3,3], correlation_list[[98]][3,3], correlation_list[[99]][3,3],
                              correlation_list[[100]][3,3], correlation_list[[101]][3,3], correlation_list[[102]][3,3],
                              correlation_list[[103]][3,3], correlation_list[[104]][3,3], correlation_list[[105]][3,3]
                 ),
                 pval = c(correlation_list[[1]][4,3], correlation_list[[2]][4,3], correlation_list[[3]][4,3],
                          correlation_list[[4]][4,3], correlation_list[[5]][4,3], correlation_list[[6]][4,3],
                          correlation_list[[7]][4,3], correlation_list[[8]][4,3], correlation_list[[9]][4,3],
                          correlation_list[[10]][4,3], correlation_list[[11]][4,3], correlation_list[[12]][4,3],
                          correlation_list[[13]][4,3], correlation_list[[14]][4,3], correlation_list[[15]][4,3],
                          correlation_list[[16]][4,3], correlation_list[[17]][4,3], correlation_list[[18]][4,3],
                          correlation_list[[19]][4,3], correlation_list[[20]][4,3], correlation_list[[21]][4,3],
                          correlation_list[[22]][4,3], correlation_list[[23]][4,3], correlation_list[[24]][4,3],
                          correlation_list[[25]][4,3], correlation_list[[26]][4,3], correlation_list[[27]][4,3],
                          correlation_list[[28]][4,3], correlation_list[[29]][4,3], correlation_list[[30]][4,3],
                          correlation_list[[31]][4,3], correlation_list[[32]][4,3], correlation_list[[33]][4,3],
                          correlation_list[[34]][4,3], correlation_list[[35]][4,3], correlation_list[[36]][4,3],
                          correlation_list[[37]][4,3], correlation_list[[38]][4,3], correlation_list[[39]][4,3],
                          correlation_list[[40]][4,3], correlation_list[[41]][4,3], correlation_list[[42]][4,3],
                          correlation_list[[43]][4,3], correlation_list[[44]][4,3], correlation_list[[45]][4,3],
                          correlation_list[[46]][4,3], correlation_list[[47]][4,3], correlation_list[[48]][4,3],
                          correlation_list[[49]][4,3], correlation_list[[50]][4,3], correlation_list[[51]][4,3],
                          correlation_list[[52]][4,3], correlation_list[[53]][4,3], correlation_list[[54]][4,3],
                          correlation_list[[55]][4,3], correlation_list[[56]][4,3], correlation_list[[57]][4,3],
                          correlation_list[[58]][4,3], correlation_list[[59]][4,3], correlation_list[[60]][4,3],
                          correlation_list[[61]][4,3], correlation_list[[62]][4,3], correlation_list[[63]][4,3],
                          correlation_list[[64]][4,3], correlation_list[[65]][4,3], correlation_list[[66]][4,3],
                          correlation_list[[67]][4,3], correlation_list[[68]][4,3], correlation_list[[69]][4,3],
                          correlation_list[[70]][4,3], correlation_list[[71]][4,3], correlation_list[[72]][4,3],
                          correlation_list[[73]][4,3], correlation_list[[74]][4,3], correlation_list[[75]][4,3],
                          correlation_list[[76]][4,3], correlation_list[[77]][4,3], correlation_list[[78]][4,3],
                          correlation_list[[79]][4,3], correlation_list[[80]][4,3], correlation_list[[81]][4,3],
                          correlation_list[[82]][4,3], correlation_list[[83]][4,3], correlation_list[[84]][4,3],
                          correlation_list[[85]][4,3], correlation_list[[86]][4,3], correlation_list[[87]][4,3],
                          correlation_list[[88]][4,3], correlation_list[[89]][4,3], correlation_list[[90]][4,3],
                          correlation_list[[91]][4,3], correlation_list[[92]][4,3], correlation_list[[93]][4,3],
                          correlation_list[[94]][4,3], correlation_list[[95]][4,3], correlation_list[[96]][4,3],
                          correlation_list[[97]][4,3], correlation_list[[98]][4,3], correlation_list[[99]][4,3],
                          correlation_list[[100]][4,3], correlation_list[[101]][4,3], correlation_list[[102]][4,3],
                          correlation_list[[103]][4,3], correlation_list[[104]][4,3], correlation_list[[105]][4,3]
                          ),
                 feature_type = c("H2AFZ","H2AFZ","H2AFZ",      
                                  "H2BK15ac","H2BK15ac","H2BK15ac",  
                                  "H3K18ac","H3K18ac","H3K18ac",   
                                  "H3K27me3","H3K27me3","H3K27me3",  
                                  "H3K56ac","H3K56ac","H3K56ac",   
                                  "H3K9me2","H3K9me2","H3K9me2",  
                                  "H4K5ac","H4K5ac","H4K5ac",
                                  "H2AK5ac","H2AK5ac","H2AK5ac",    
                                  "H2BK20ac","H2BK20ac","H2BK20ac",  
                                  "H3K20me1","H3K20me1","H3K20me1", 
                                  "H3K36me3","H3K36me3","H3K36me3", 
                                  "H3K79me1","H3K79me1","H3K79me1", 
                                  "H3K9me3","H3K9me3","H3K9me3",
                                  "H4K8ac","H4K8ac","H4K8ac",
                                  "H2AK9ac","H2AK9ac","H2AK9ac",  
                                  "H2BK5ac","H2BK5ac","H2BK5ac",  
                                  "H3K23ac","H3K23ac","H3K23ac",
                                  "H3K4ac","H3K4ac","H3K4ac",   
                                  "H3K79me2","H3K79me2", "H3K79me2",
                                  "H3T11ph","H3T11ph","H3T11ph",  
                                  "H4K91ac","H4K91ac","H4K91ac",
                                  "H2BK120ac","H2BK120ac","H2BK120ac",  
                                  "H3F3A","H3F3A","H3F3A",    
                                  "H3K23me2","H3K23me2","H3K23me2", 
                                  "H3K4me1","H3K4me1","H3K4me1",  
                                  "H3K9ac","H3K9ac","H3K9ac", 
                                  "H4K12ac","H4K12ac","H4K12ac",
                                  "H2BK12ac","H2BK12ac","H2BK12ac", 
                                  "H3K14ac","H3K14ac","H3K14ac",  
                                  "H3K27ac","H3K27ac","H3K27ac",  
                                  "H3K4me3","H3K4me3","H3K4me3",  
                                  "H3K9me1","H3K9me1","H3K9me1",  
                                  "H4K4me2","H4K4me2","H4K4me2",
                                  "Chromatin Accessibility","Chromatin Accessibility","Chromatin Accessibility",
                                  "Methylome","Methylome","Methylome"))

#df <- data.frame(sequence_type = c("lncRNA Exons","Short ncRNAs","Protein-coding Exons"),
#                 correlation = c(correlation_list[[1]][1,3], correlation_list[[2]][1,3], correlation_list[[3]][1,3]
#                                
#                 ),
#                 lower_ci = c(correlation_list[[1]][2,3], correlation_list[[2]][2,3], correlation_list[[3]][2,3]
#                            
#                 ),
#                 upper_ci = c(correlation_list[[1]][3,3], correlation_list[[2]][3,3], correlation_list[[3]][3,3]
#                              
#                 ),
#                 feature_type = c("Chromatin Accessibility","Chromatin Accessibility","Chromatin Accessibility"
#                                  ))

# Create a factor with the correct order for the x-axis
df$feature_name <- factor(df$sequence_type, levels = unique(df$sequence_type))

# Create a factor for sequence_type with the correct order for shapes
df$sequence_type <- factor(df$sequence_type, 
                           levels = c("lncRNA Exons","Short ncRNAs","Protein-coding Exons"))

# Create a factor for feature_type with the correct order for colors
df$feature_type <- factor(df$feature_type, 
                          levels = unique(df$feature_type))

ggplot(df, aes(x = feature_type, y = correlation, color = feature_type, shape = sequence_type, alpha = ifelse(pval > 0.05, 0.4, 1))) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.05, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(19, 17, 15)) +  #c(19, 17, 15))
  scale_color_manual(values = c("blue","green","pink4","orange4","brown","purple","red","darkgreen","darkorchid4",
                                "gold","lightcoral","lightblue4","orange3","maroon4","purple3","khaki4","limegreen","hotpink4",
                                "gold4","forestgreen","firebrick","dodgerblue4","deeppink4","darkviolet","darkslategrey","darksalmon","darkolivegreen4","darkcyan",
                                "yellowgreen","wheat4","violetred4","turquoise4","tomato4","chartreuse4","cadetblue4")) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),  # Rotate x-axis labels
        legend.position = "right",
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(title = "Spearman correlation Histone and DNA modifications Analysis",
       x = "",
       y = "Spearman correlation",
       color = "Feature Type",
       shape = "Sequence Type") + 
  ylim(-1, 1) + 
  scale_alpha_continuous(range = c(0.4, 1),
                         guide = "none")
  #geom_hline(yintercept = 0.1, linetype = "dotted") + 
  #geom_hline(yintercept = -0.1, linetype = "dotted") 
###############################



# Find out range for Histone marks
histone_data <- read.csv("temp_histone_marks.bed", header=FALSE, sep = "\t")
summary(histone_data$V7)

histone_data <- read.csv("../data/histone_feature/lncrna-histone-feature-matrix.csv", header=TRUE)
summary(histone_data$AvgSignal)


#######################################
# Statistics for methylation files:
methyl_file <- "../data/ENCFF913ZNZ.bed"
methyl_data <- read.csv(methyl_file, header = FALSE, sep = "\t", nrows = 100000)
summary(methyl_data)

