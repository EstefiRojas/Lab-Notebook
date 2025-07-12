# Load necessary library
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
library(boot) # Used to obtain confidence intervals

args <- commandArgs(trailingOnly = TRUE)
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


#POL2R2A
lncrna_file9 <- "../data/pol2r2a_feature/lncrna_matrix.csv"
#output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"

sncrna_file9 <- "../data/pol2r2a_feature/short_ncrna_matrix.csv"
#output_file <- "../data/histone_feature/protein_spearman_matrix.csv"

protein_file9 <- "../data/pol2r2a_feature/protein_matrix.csv"
#output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"

# Initialize empty list to store correlation vectors
correlation_list <- list()

correlation_list[[lncrna_file9]] <- spearman_correlation(lncrna_file9)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[sncrna_file9]] <- spearman_correlation(sncrna_file9)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[protein_file9]] <- spearman_correlation(protein_file9)
#correlation_matrix_combined <- do.call(rbind, correlation_list)


correlation_list[[lncrna_file9]]
######################################
# Create a graph to visualice features
df <- data.frame(sequence_type = c("lncRNA Exons","Short ncRNAs","Protein-coding Exons"),
                 correlation = c(correlation_list[[lncrna_file9]][1,1], correlation_list[[sncrna_file9]][1,1], correlation_list[[protein_file9]][1,1]
                                 ),
                 lower_ci = c(correlation_list[[lncrna_file9]][2,1], correlation_list[[sncrna_file9]][2,1], correlation_list[[protein_file9]][2,1]
                              ),
                 upper_ci = c(correlation_list[[lncrna_file9]][3,1], correlation_list[[sncrna_file9]][3,1], correlation_list[[protein_file9]][3,1]
                              ),
                 pval = c(correlation_list[[lncrna_file9]][4,1], correlation_list[[sncrna_file9]][4,1], correlation_list[[protein_file9]][4,1]
                          ),
                 feature_type = c("POL2R2A", "POL2R2A", "POL2R2A"
                                  ))

# Create a factor with the correct order for the x-axis
df$feature_name <- factor(df$sequence_type, levels = unique(df$sequence_type))

# Create a factor for sequence_type with the correct order for shapes
df$sequence_type <- factor(df$sequence_type, 
                           levels = c("lncRNA Exons","Short ncRNAs","Protein-coding Exons"))

# Create a factor for feature_type with the correct order for colors
df$feature_type <- factor(df$feature_type, 
                          levels = c("POL2R2A"))

ggplot(df, aes(x = feature_type, y = correlation, color = feature_type, shape = sequence_type, alpha = ifelse(pval > 0.05, 0.3, 1))) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.05, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(19, 17, 15)) +  #c(19, 17, 15))
  scale_color_manual(values = c("green4","pink4","orange4","brown","orchid4","red4","darkviolet","cyan4","darkgreen","salmon4","chartreuse4","aquamarine4","darkblue","darkred","darkslategrey")) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels
        legend.position = "right") +
  labs(title = "Spearman correlation Intrinsic Sequence Analysis",
       x = "",
       y = "Spearman correlation",
       color = "Feature Type",
       shape = "Sequence Type") +
  scale_alpha_continuous(range = c(0.3, 1),
                         guide = "none")
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

