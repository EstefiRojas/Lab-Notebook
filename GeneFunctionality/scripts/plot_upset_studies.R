# Load necessary library
if (!require("UpSetR")) {
    install.packages("UpSetR", repos = "http://cran.us.r-project.org")
}
library(UpSetR)

# Define paths
input_file <- "/Volumes/ADATA HD710 PRO/Downloads/Estefi/Otago University/Lab-Notebook/GeneFunctionality/results/plots2/upset_binary_matrix.csv"
output_file <- "/Volumes/ADATA HD710 PRO/Downloads/Estefi/Otago University/Lab-Notebook/GeneFunctionality/results/plots2/upset_plot_studies.png"

# Check if file exists
if (!file.exists(input_file)) {
    stop(paste("Input file not found:", input_file))
}

# Read data
print(paste("Reading matrix from:", input_file))
upset_data <- read.csv(input_file, row.names = 1)

# Identify studies and aggregate
print("Aggregating data by Study...")

# Extract unique study names (assumes format "Study_Category")
# We take the first part of the string before the first underscore
all_cols <- colnames(upset_data)
studies <- unique(sapply(strsplit(all_cols, "_"), `[`, 1))

print(paste("Found studies:", paste(studies, collapse = ", ")))

# Create a new dataframe for aggregated data
study_data <- data.frame(row.names = rownames(upset_data))

for (study in studies) {
    # Find columns belonging to this study
    study_cols <- grep(paste0("^", study, "_"), all_cols, value = TRUE)

    if (length(study_cols) > 0) {
        # Check if row has any non-zero value in these columns
        if (length(study_cols) == 1) {
            # If only one column, just take it (though binary check ensures 0/1)
            study_data[[study]] <- ifelse(upset_data[[study_cols]] > 0, 1, 0)
        } else {
            # Row sum > 0 means present in at least one category
            study_data[[study]] <- ifelse(rowSums(upset_data[, study_cols]) > 0, 1, 0)
        }
    }
}

# Generate UpSet Plot
print("Generating Study-Level UpSet plot...")
png(filename = output_file, width = 1600, height = 900, res = 150)

upset(study_data,
    nsets = ncol(study_data),
    nintersects = NA,
    order.by = "freq",
    show.numbers = "yes",
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    # Using the same text scale as the main plot
    text.scale = c(1.5, 1.3, 1.5, 1.3, 1.3, 1),
    set_size.show = TRUE,
    set_size.numbers_size = 6,
    set_size.scale_max = 6400 # Adding a bit more buffer as sums will be larger
)

dev.off()
print(paste("Plot saved to:", output_file))
