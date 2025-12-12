# Load necessary library
if (!require("UpSetR")) {
    install.packages("UpSetR", repos = "http://cran.us.r-project.org")
}
library(UpSetR)

# Define paths
input_file <- "/Volumes/ADATA HD710 PRO/Downloads/Estefi/Otago University/Lab-Notebook/GeneFunctionality/results/plots2/upset_binary_matrix.csv"
output_file <- "/Volumes/ADATA HD710 PRO/Downloads/Estefi/Otago University/Lab-Notebook/GeneFunctionality/results/plots2/upset_plot_R.png"

# Check if file exists
if (!file.exists(input_file)) {
    stop(paste("Input file not found:", input_file))
}

# Read data
# Read data
print(paste("Reading matrix from:", input_file))
upset_data <- read.csv(input_file, row.names = 1, check.names = FALSE)

# Clean column names (replace "_" and "." with " - ")
colnames(upset_data) <- gsub("_", " - ", colnames(upset_data))
colnames(upset_data) <- gsub("\\.", " - ", colnames(upset_data))

# Define Desired Order
studies <- c("Huang", "Liang", "Liu", "Montero") # Order of studies
categories <- c("Core", "Common", "Rare", "Non-essential") # Order of categories

ordered_sets <- c()
# Iterate Categories FIRST, then Studies
for (cat in categories) {
    for (study in studies) {
        set_name <- paste(study, "-", cat)
        if (set_name %in% colnames(upset_data)) {
            ordered_sets <- c(ordered_sets, set_name)
        }
    }
}

# Inverse Order as requested
ordered_sets <- rev(ordered_sets)
# Generate UpSet Plot
print("Generating UpSet plot...")
png(filename = output_file, width = 2000, height = 900, res = 150)

upset(upset_data,
    sets = ordered_sets,
    keep.order = TRUE,
    nintersects = NA,
    order.by = "freq",
    show.numbers = "yes",
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    # Adjusted text scale: c(intersection size title, intersection size tick
    # labels,set size title, set size tick labels, set names, numbers above
    # bars)
    text.scale = c(1.5, 1.3, 1.5, 1.3, 1.3, 1),
    set_size.show = TRUE,
    set_size.numbers_size = 6,
    set_size.scale_max = 6400
)

dev.off()
print(paste("Plot saved to:", output_file))
