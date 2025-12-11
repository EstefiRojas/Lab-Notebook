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

# Generate UpSet Plot
print("Generating UpSet plot...")
png(filename = output_file, width = 1600, height = 900, res = 150)

upset(plot_data,
    nsets = ncol(plot_data),
    nintersects = NA,
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

# Generate Pairwise Plots
print("Generating Pairwise UpSet plots...")

pairs <- combn(target_studies, 2, simplify = FALSE)

for (pair in pairs) {
    study1 <- pair[1]
    study2 <- pair[2]

    print(paste("Processing pair:", study1, "vs", study2))

    # Select columns for this pair
    cols_to_plot <- c(
        paste(study1, "- Essential"),
        paste(study1, "- Non-essential"),
        paste(study2, "- Essential"),
        paste(study2, "- Non-essential")
    )

    # Check if columns exist
    missing_cols <- setdiff(cols_to_plot, colnames(plot_data))
    if (length(missing_cols) > 0) {
        warning(paste("Missing columns for pair", study1, "vs", study2, ":", paste(missing_cols, collapse = ", ")))
        next
    }

    pair_data <- plot_data[, cols_to_plot]

    # Output filename
    pair_output_file <- paste0(dirname(output_file), "/upset_pairwise_", study1, "_vs_", study2, ".png")

    png(filename = pair_output_file, width = 1200, height = 900, res = 150)

    upset(pair_data,
        nsets = 4,
        nintersects = NA,
        order.by = "freq",
        show.numbers = "yes",
        mainbar.y.label = "Intersection Size",
        sets.x.label = "Set Size",
        text.scale = c(1.5, 1.3, 1.5, 1.3, 1.3, 1),
        set_size.show = TRUE,
        set_size.numbers_size = 6,
        set_size.scale_max = max(colSums(pair_data)) + 500
    )

    dev.off()
    print(paste("Saved pairwise plot to:", pair_output_file))
}
