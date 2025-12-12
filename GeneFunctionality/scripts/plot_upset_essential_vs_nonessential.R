# Load necessary library
if (!require("UpSetR")) {
    install.packages("UpSetR", repos = "http://cran.us.r-project.org")
}
library(UpSetR)

# Define paths
input_file <- "/Volumes/ADATA HD710 PRO/Downloads/Estefi/Otago University/Lab-Notebook/GeneFunctionality/results/plots2/upset_binary_matrix.csv"
output_file <- "/Volumes/ADATA HD710 PRO/Downloads/Estefi/Otago University/Lab-Notebook/GeneFunctionality/results/plots2/upset_essential_vs_nonessential.png"

# Check if file exists
if (!file.exists(input_file)) {
    stop(paste("Input file not found:", input_file))
}

# Read data
print(paste("Reading matrix from:", input_file))
upset_data <- read.csv(input_file, row.names = 1, check.names = FALSE)

# Targeted studies
target_studies <- c("Huang", "Liang", "Montero")

# Initialize new dataframe for plot data
plot_data <- data.frame(row.names = rownames(upset_data))

print("Processing Essential and Non-essential groups for target studies...")

all_cols <- colnames(upset_data)

for (study in target_studies) {
    # 1. Essential: Combine Rare, Common, Core
    # Find columns for this study that match the essential categories
    essential_cols <- grep(paste0("^", study, "_(Rare|Common|Core)$"), all_cols, value = TRUE)

    if (length(essential_cols) > 0) {
        # If any of these are 1, then "Study - Essential" is 1
        col_name <- paste(study, "- Essential")
        if (length(essential_cols) == 1) {
            plot_data[[col_name]] <- upset_data[[essential_cols]]
        } else {
            plot_data[[col_name]] <- ifelse(rowSums(upset_data[, essential_cols]) > 0, 1, 0)
        }
        print(paste("Created:", col_name, "from", paste(essential_cols, collapse = ", ")))
    } else {
        warning(paste("No essential columns found for study:", study))
    }

    # 2. Non-essential: Just the Non-essential column
    non_essential_col <- grep(paste0("^", study, "_Non-essential$"), all_cols, value = TRUE)

    if (length(non_essential_col) == 1) {
        col_name <- paste(study, "- Non-essential")
        plot_data[[col_name]] <- upset_data[[non_essential_col]]
        print(paste("Created:", col_name))
    } else {
        warning(paste("No Non-essential column found for study:", study))
    }
}

# Define Desired Order for Sets
ordered_sets <- c()

# 1. Add all Essential sets first
for (study in target_studies) {
    set_e <- paste(study, "- Essential")
    if (set_e %in% colnames(plot_data)) ordered_sets <- c(ordered_sets, set_e)
}

# 2. Add all Non-essential sets second
for (study in target_studies) {
    set_ne <- paste(study, "- Non-essential")
    if (set_ne %in% colnames(plot_data)) ordered_sets <- c(ordered_sets, set_ne)
}

# Inverse Order as requested
ordered_sets <- rev(ordered_sets)

print("Set Order:")
print(ordered_sets)

# Generate UpSet Plot
print("Generating UpSet plot...")
png(filename = output_file, width = 1600, height = 900, res = 150)

upset(plot_data,
    sets = ordered_sets,
    keep.order = TRUE,
    order.by = "freq",
    show.numbers = "yes",
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size",
    text.scale = c(1.5, 1.3, 1.5, 1.3, 1.3, 1),
    set_size.show = TRUE,
    set_size.numbers_size = 6,
    set_size.scale_max = max(colSums(plot_data)) + 500 # Dynamic scale max
)

dev.off()
print(paste("Plot saved to:", output_file))
