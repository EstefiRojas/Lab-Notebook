###############
# Script to Compute z-scores for all gene functionality features separated by Gene Type.
# The aim is to explore the resulting z-scores and determine the reason for some 
# extremely big values obtained for some features.
# Can be executed from the command line using:
#
#  Rscript compute_zscores_review02.R ../data/features/gene_functionality_features_latest1000all_extra.csv --verbose --subset=100
###############

library(dplyr)
library(ggplot2)

# Function to compute robust z-scores
compute_robust_zscores_optimized <- function(positive_matrix, negative_matrix, verbose = FALSE) {
  if (verbose) cat("Computing robust z-scores...\n")  # Debug log
  
  # Use NEGATIVE CONTROLS median absolute deviations (MAD) and medians to compute robust z-scores
  mad_values <- apply(negative_matrix, 2, mad, na.rm = TRUE)  # Calculate MAD for each column in the negative control set
  median_values <- apply(negative_matrix, 2, median, na.rm = TRUE)  # Calculate median for each column in the negative control set
  
  # Log MAD and median values for debugging
  if (verbose) {
    cat("MAD values:\n")
    print(mad_values)
    cat("Median values:\n")
    print(median_values)
  }
  
  # Handle case where MAD is 0 or very small by using meanAD instead
  mad_values <- apply(negative_matrix, 2, function(col) {
    mad_value <- mad(col, na.rm = TRUE)
    if (mad_value <= 1e-6) {  # If MAD is 0 or very small, calculate mean absolute deviation (meanAD)
      meanAD <- mean(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE)
      if (meanAD <= 1e-6) {  # If meanAD is also very small, set a floor value to avoid division by extremely small numbers
        warning(sprintf("Column %s with very small MAD and meanAD encountered: Setting value to 1e-6.\n", col))
        return(1e-6)  # Set MAD to a small constant value to prevent instability in z-score calculation
      } else {
        warning(sprintf("Column %s with very small MAD of %f encountered: Using meanAD value of %f.\n", col, mad_value, meanAD))
        return(1.2533 * meanAD)  # Use meanAD if it is greater than the threshold
      }
    } else {
      return(mad_value)  # Use MAD if it is sufficient
    }
  })
  
  # Log final MAD values after adjustments
  if (verbose) {
    cat("Adjusted MAD values:\n")
    print(mad_values)
  }
  
  # Calculate z-scores: subtract the median and divide by MAD for each feature
  zscores <- sweep(positive_matrix, 2, median_values, "-")  # Center each column by subtracting the median
  zscores <- sweep(zscores, 2, mad_values, "/")  # Scale each column by dividing by MAD
  
  # Log computed z-scores for debugging
  if (verbose) {
    cat("Computed z-scores:\n")
    print(head(zscores))
  }
  
  zscores
}

# Function to visualize data
visualize_data <- function(df, col_name, title, plot_type = "hist", verbose = FALSE, 
                           graph_color = NA, zoom = FALSE, y_axis_title = "Robust z-score") {
  if (verbose) cat(sprintf("Visualizing data for column: %s\n", col_name))  # Debug log
  
  # Create a dataframe with the selected column and remove any NA values
  original_data <- na.omit(data.frame(Value = df[[col_name]], Status = "Original"))
  
  # Log the number of values being visualized
  if (verbose) cat(sprintf("Number of non-NA values in column %s: %d\n", col_name, nrow(original_data)))
  
  # Generate plot based on the specified plot type (histogram or density plot)
  p <- ggplot(original_data, aes(x = Value, fill = Status)) +
    {if (plot_type == "hist") geom_histogram(alpha = 0.6, bins = 100) 
      else if(plot_type == "density") geom_density(alpha = 0.6)
      else geom_boxplot(aes(x = Status, y = Value))} +
    {if(plot_type != "boxplot") labs(x = y_axis_title) else labs(x = "", y = y_axis_title)} +
    {if(!is.na(graph_color))  scale_fill_manual(values = c(graph_color))} +
    ggtitle(paste(title, col_name)) +
    {if(zoom & plot_type == "boxplot") coord_cartesian(ylim =  c(0, 200)) } +
    theme_minimal(base_size = 64) +
    theme(axis.text.x = element_text(size = 58),
          axis.text.y = element_text(size = 58),  # Increase y-axis text size
          axis.title.x = element_text(size = 64),  # Increase x-axis title size
          axis.title.y = element_text(size = 64),  # Increase y-axis title size
          legend.position = "none",
          legend.title = element_text(size = 62),  # Increase legend title size
          legend.text = element_text(size = 68),  # Increase legend text size
          plot.title = element_text(size = 74, hjust = 0.5), # Increase plot title size and center it
          plot.subtitle = element_blank(),
          strip.text = element_text(size = 0, margin = margin(0,0,40,0))) # Remove the legend for a cleaner plot
  
  # Save the plot to a file to ensure it is rendered in non-interactive environments
  plot_file <- paste0(gsub(" ", "_", title), "_", col_name, "_", plot_type, "_", nrow(original_data), ".png")
  plot_path <-  "../results/latest1000all/codeReview2/"
  #ggsave(plot_file, plot = p)
  ggsave(plot_file, plot = p, path = plot_path, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 300)
  if (verbose) cat(sprintf("Plot saved to %s\n", plot_file))
}

# Main function
main <- function() {
  # Capture command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Please provide the path to the data file as an argument.")  # Ensure a file path is provided
  }
  
  # Check if verbose option is enabled
  verbose <- "--verbose" %in% args
  
  # Load and preprocess data
  data_file <- args[1]  # Get the file path from the command line arguments
  if (verbose) cat(sprintf("Loading data from file: %s\n", data_file))  # Debug log
  data <- read.csv(data_file, header = TRUE, check.names = TRUE)# Load the data from the specified CSV file
  
  # Check if subset option is provided
  subset_value <- as.numeric(gsub("--subset=", "", args[grep("--subset=", args)]))
  if (is.na(subset_value)) {
    subset_value <- nrow(data)  # Use nrow to indicate no subsetting by default
  }
  
  # Log the first few rows of the loaded data
  if (verbose) {
    cat("Loaded data preview:\n")
    print(head(data$RPKM_tissue))
  }
  
  # Convert to numeric
  data_numeric <- data %>% 
    dplyr::select(-X, -Dataset, -ID, -Functional, -Chromosome, -Start, -End, -Sequence) %>%
    sapply(function(feature) as.numeric(as.character(feature))) %>%
    as.data.frame()
  
  # Restore Dataset column
  data_numeric$Dataset <- data$Dataset

  # Log the first few rows of the loaded data
  if (verbose) {
    cat("\n#######################\n")
    cat("Numeric data preview:\n")
    print(head(data_numeric$RPKM_tissue))
  }
  
  # Separate the data into positive and negative feature matrices
  if (verbose) cat("\n#######################\n")
  if (verbose) cat("Separating data into positive and negative feature matrices...\n")  # Debug log
  protein_positive_feature_matrix <- data_numeric %>% 
    filter(Dataset == "protein-coding-exon2" | Dataset == "protein-coding-exon3") %>% 
    select(RPKM_tissue, -Dataset)
  protein_negative_feature_matrix <- data_numeric %>% 
    filter(Dataset == "protein-exon2-negative-control" | Dataset == "protein-exon3-negative-control") %>% 
    select(RPKM_tissue, -Dataset)
  
  # Select random samples from positive and negative datasets
  set.seed(197)  # Set seed for reproducibility
  protein_positive_feature_matrix <- protein_positive_feature_matrix %>% sample_n(min(subset_value, nrow(protein_positive_feature_matrix)))  # Subset rows from positive data
  protein_negative_feature_matrix <- protein_negative_feature_matrix %>% sample_n(min(10 * subset_value, nrow(protein_negative_feature_matrix)))  # Subset rows from negative data

  # Log the dimensions of the matrices
  if (verbose) {
    cat("\n#######################\n")
    cat(sprintf("Dimensions of positive feature matrix: %d rows, %d columns\n", nrow(protein_positive_feature_matrix), ncol(protein_positive_feature_matrix)))
    cat(sprintf("Dimensions of negative feature matrix: %d rows, %d columns\n", nrow(protein_negative_feature_matrix), ncol(protein_negative_feature_matrix)))
  }
  
  # Compute summary statistics without outlier removal
  if (verbose) cat("\n#######################\n")
  if (verbose) cat("Summary Statistics:\n")
  list(positive = protein_positive_feature_matrix, negative = protein_negative_feature_matrix) %>%
    lapply(summary) %>% print  # Print summary statistics for both positive and negative matrices
  
  # Visualize RPKM_tissue
  if ("RPKM_tissue" %in% colnames(protein_positive_feature_matrix)) {
    if (verbose) cat("\n#######################\n")
    if (verbose) cat("Visualizing RPKM_tissue...\n")
    print(visualize_data(protein_positive_feature_matrix, "RPKM_tissue", "Protein Functional Raw", plot_type = "density", verbose, y_axis_title = "RPKM tissue"))  # Generate and print the density plot 
    print(visualize_data(protein_positive_feature_matrix, "RPKM_tissue", "Protein Functional Raw", plot_type = "hist", verbose, y_axis_title = "RPKM tissue"))  # Generate and print the histogram plot
    print(visualize_data(protein_positive_feature_matrix, "RPKM_tissue", "Protein Functional Raw", plot_type = "boxplot", verbose, y_axis_title = "RPKM tissue"))  # Generate and print the boxplot plot
    
    print(visualize_data(protein_negative_feature_matrix, "RPKM_tissue", "Protein Negative Control Raw", plot_type = "density", verbose, graph_color="blue", y_axis_title = "RPKM tissue"))  # Generate and print the density plot 
    print(visualize_data(protein_negative_feature_matrix, "RPKM_tissue", "Protein Negative Control Raw", plot_type = "hist", verbose, graph_color="blue", y_axis_title = "RPKM tissue"))  # Generate and print the histogram plot
    print(visualize_data(protein_negative_feature_matrix, "RPKM_tissue", "Protein Negative Control Raw", plot_type = "boxplot", verbose, graph_color="blue", y_axis_title = "RPKM tissue"))  # Generate and print the boxplot plot
  } else {
    warning("RPKM_tissue column not found in the z-scores dataframe.")  # Warn if the column is missing
  }
  
  # Compute z-scores
  if (verbose) cat("\n#######################\n")
  if (verbose) cat("Computing z-scores for positive feature matrix...\n")
  protein_functional_z_scores <- compute_robust_zscores_optimized(protein_positive_feature_matrix, protein_negative_feature_matrix, verbose)
  
  if (verbose) cat("\n#######################\n")
  if (verbose) cat("Computing z-scores for negative feature matrix...\n")
  protein_negative_z_scores <- compute_robust_zscores_optimized(protein_negative_feature_matrix, protein_negative_feature_matrix, verbose)
  
  # Print z-score summary statistics
  if (verbose) {
    cat("\n#######################\n")
    cat("Functional Z-Scores Summary:\n")
    print(summary(protein_functional_z_scores))  # Print summary statistics for the computed z-scores
    cat("\n#######################\n")
    cat("Negative Control Z-Scores Summary:\n")
    print(summary(protein_negative_z_scores))
  }
  
  # Save z-scores to CSV
  z_scores_file <- "z_scores/protein_z_scores_output.csv"
  write.csv(rbind(protein_functional_z_scores, protein_negative_z_scores), z_scores_file, row.names = FALSE)
  if (verbose) cat(sprintf("Z-scores saved to %s\n", z_scores_file))
  
  
  # Visualize RPKM_tissue
  if ("RPKM_tissue" %in% colnames(protein_functional_z_scores)) {
    if (verbose) cat("\n#######################\n")
    if (verbose) cat("Visualizing RPKM_tissue...\n")
    print(visualize_data(protein_functional_z_scores, "RPKM_tissue", "Protein Functional", plot_type = "density", verbose))  # Generate and print the density plot 
    print(visualize_data(protein_functional_z_scores, "RPKM_tissue", "Protein Functional", plot_type = "hist", verbose))  # Generate and print the density plot
    print(visualize_data(protein_functional_z_scores, "RPKM_tissue", "Protein Functional", plot_type = "boxplot", verbose))  # Generate and print the density plot
    
    print(visualize_data(protein_negative_z_scores, "RPKM_tissue", "Protein Negative Control", plot_type = "density", verbose, graph_color="blue"))  # Generate and print the density plot 
    print(visualize_data(protein_negative_z_scores, "RPKM_tissue", "Protein Negative Control", plot_type = "hist", verbose, graph_color="blue"))  # Generate and print the density plot
    print(visualize_data(protein_negative_z_scores, "RPKM_tissue", "Protein Negative Control", plot_type = "boxplot", verbose, graph_color="blue"))  # Generate and print the density plot
  } else {
    warning("RPKM_tissue column not found in the z-scores dataframe.")  # Warn if the GC.content column is missing
  }
}

# Run the main function
main()

