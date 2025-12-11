# Load necessary libraries
if (!require("pheatmap")) {
    install.packages("pheatmap", repos = "http://cran.us.r-project.org")
}
library(pheatmap)
library(grid)

# Define paths
input_file <- "/Volumes/ADATA HD710 PRO/Downloads/Estefi/Otago University/Lab-Notebook/GeneFunctionality/results/plots2/upset_binary_matrix.csv"
output_dir <- "/Volumes/ADATA HD710 PRO/Downloads/Estefi/Otago University/Lab-Notebook/GeneFunctionality/results/plots2/"

# Check if file exists
if (!file.exists(input_file)) {
    stop(paste("Input file not found:", input_file))
}

# Read data
print(paste("Reading matrix from:", input_file))
upset_data <- read.csv(input_file, row.names = 1, check.names = FALSE)

# Targeted studies
target_studies <- c("Huang", "Liang", "Montero")

# Function to get classification for a study
# Returns a vector: "Essential", "Non-essential", or NA (if neither) for each gene
get_study_classification <- function(df, study_name) {
    all_cols <- colnames(df)

    # Identify Essential Columns
    ess_cols <- grep(paste0("^", study_name, "_(Rare|Common|Core)$"), all_cols, value = TRUE)

    # Identify Non-essential Column
    non_ess_col <- grep(paste0("^", study_name, "_Non-essential$"), all_cols, value = TRUE)

    classification <- rep(NA, nrow(df))

    # Check Essential
    if (length(ess_cols) > 0) {
        is_essential <- rowSums(df[, ess_cols, drop = FALSE]) > 0
        classification[is_essential] <- "Essential"
    }

    # Check Non-essential
    if (length(non_ess_col) == 1) {
        # Note: A gene could technically be both if the input data allows overlap (though strictly it shouldn't)
        # If it's already marked Essential, we keep it as Essential (Priority) or flag as Ambiguous?
        # Assuming mutual exclusivity or Essential priority for this analysis:
        is_non_essential <- df[[non_ess_col]] > 0

        # Assign "Non-essential" only if not already "Essential"
        classification[is_non_essential & is.na(classification)] <- "Non-essential"
    }

    return(classification)
}

print("Running Pairwise Fisher's Exact Tests...")

pairs <- combn(target_studies, 2, simplify = FALSE)

for (pair in pairs) {
    study1 <- pair[1]
    study2 <- pair[2]

    print(paste("Processing:", study1, "vs", study2))

    # Get classifications
    class1 <- get_study_classification(upset_data, study1)
    class2 <- get_study_classification(upset_data, study2)

    # Create data frame for this pair
    pair_df <- data.frame(
        S1 = class1,
        S2 = class2,
        stringsAsFactors = FALSE
    )

    # Filter: Keep only genes that have a classification in BOTH studies
    # (i.e., ignore genes that are completely missing/not analyzed in one of the studies)
    valid_pair <- pair_df[!is.na(pair_df$S1) & !is.na(pair_df$S2), ]

    print(paste("  Genes matching criteria:", nrow(valid_pair)))

    if (nrow(valid_pair) == 0) {
        warning(paste("No overlapping genes found for", study1, "vs", study2))
        next
    }

    # Create Contingency Table (2x2)
    # Factor levels enforce the order: Essential, Non-essential
    tbl <- table(
        factor(valid_pair$S1, levels = c("Essential", "Non-essential")),
        factor(valid_pair$S2, levels = c("Essential", "Non-essential"))
    )

    # Run Fisher's Exact Test
    ft <- fisher.test(tbl)
    p_val <- ft$p.value
    or <- ft$estimate

    print(paste("  P-value:", p_val))
    print("  Table:")
    print(tbl)

    # Prepare Data for Heatmap (Counts)
    # pheatmap expects a matrix.
    plot_matrix <- as.matrix(tbl)

    # Proper Labels
    rownames(plot_matrix) <- paste(study1, rownames(plot_matrix))
    colnames(plot_matrix) <- paste(study2, colnames(plot_matrix))

    # Output Filename
    output_plot <- paste0(output_dir, "fisher_counts_", study1, "_vs_", study2, ".png")

    # Color palette (Light Blue to Dark Blue)
    colors_list <- colorRampPalette(c("#cfcfcf", "#08519C"))(50)

    # Title with Stats
    plot_title <- paste0(
        "Fisher's Exact Test: ", study1, " vs ", study2, "\n",
        "p-value = ", sprintf("%.3e", p_val),
        " | Odds Ratio = ", sprintf("%.2f", or)
    )

    png(filename = output_plot, width = 1000, height = 800, res = 150)

    pheatmap(plot_matrix,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        display_numbers = TRUE, # Show Counts
        number_format = "%.0f", # Integer format
        fontsize_number = 20,
        number_color = if (max(plot_matrix) > 50) "white" else "black", # dynamic readability
        main = plot_title,
        color = colors_list,
        legend = TRUE,
        angle_col = 0
    )

    dev.off()
    print(paste("  Saved plot to:", output_plot))
}

print("All pairwise comparisons completed.")
