# Script to load functional genetic features data into a dataframe.

################################
# Load features in a dataframe #
################################
file <- "../data/features/gene_functionality_features_latest1000all_new_din.csv"
feature_matrix <- read.csv(file, header = TRUE, check.names = FALSE)

##############################
# Retrieve genetype matrices #
##############################
func_prot_exon2_data <- feature_matrix[feature_matrix$Dataset ==
                                         "protein-coding-exon2", ]
func_prot_exon3_data <- feature_matrix[feature_matrix$Dataset ==
                                         "protein-coding-exon3", ]
func_lncrna_exon1_data <- feature_matrix[feature_matrix$Dataset ==
                                           "lncrna-exon1", ]
func_lncrna_exon2_data <- feature_matrix[feature_matrix$Dataset ==
                                           "lncrna-exon2", ]
func_sncrna_dataset <- feature_matrix[feature_matrix$Dataset ==
                                        "short-ncrna", ]

prot_exon2_nc_data <- feature_matrix[feature_matrix$Dataset ==
                                       "protein-exon2-negative-control", ]
prot_exon3_nc_data <- feature_matrix[feature_matrix$Dataset ==
                                       "protein-exon3-negative-control", ]
lncrna_exon1_nc_data <- feature_matrix[feature_matrix$Dataset ==
                                         "lncrna-exon1-negative-control", ]
lncrna_exon2_nc_data <- feature_matrix[feature_matrix$Dataset ==
                                         "lncrna-exon2-negative-control", ]
sncrna_nc_data <- feature_matrix[feature_matrix$Dataset ==
                                   "short-ncrna-negative-control", ]


# Function to transform features to numeric datatype
convert_to_numeric <- function(column) {
  as.numeric(as.character(column))
}

# Function to get a numeric version of the received dataframe
get_numeric_dataframe <- function(dataframe, select_features,
                                  omit_nas = FALSE) {
  # Store functional and dataset features in independent vectors
  functional_column <- dataframe$Functional
  dataset_column <- dataframe$Dataset

  # Select only specified features
  selected_features_matrix <- dataframe[, select_features]
  # If a dataframe without NAs is requested
  if (omit_nas) {
    dataframe_no_nas <- dataframe[, c(select_features, "Functional", "Dataset")]
    # Remove NAs
    dataframe_no_nas <- na.omit(dataframe_no_nas)
    # Store functional and dataset features in separate vectors
    functional_column <- dataframe_no_nas$Functional
    dataset_column <- dataframe_no_nas$Dataset
    # Select only specified features replacing with no nas version
    selected_features_matrix <- dataframe_no_nas[, select_features]
  }

  # Convert to numeric
  feature_matrix_numeric <- as.data.frame(sapply(selected_features_matrix,
                                                 convert_to_numeric))
  # Restore functional and dataset columns
  feature_matrix_numeric$Functional <- functional_column
  feature_matrix_numeric$Dataset <- dataset_column

  return(feature_matrix_numeric)
}

get_numeric_features <- function(selected_features) {
  feature_matrix_numeric_no_nas <- get_numeric_dataframe(feature_matrix,
                                                         selected_features,
                                                         omit_nas = TRUE)
  feature_matrix_numeric <- get_numeric_dataframe(feature_matrix,
                                                  selected_features)

  return(list("feature_matrix_numeric_no_nas" = feature_matrix_numeric_no_nas,
              "feature_matrix_numeric" = feature_matrix_numeric))
}
