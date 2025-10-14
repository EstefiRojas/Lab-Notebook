# ============================================================================
# RNA PCA Analysis: Diagnostic Follow-up Code
# ADJUSTED FOR PRE-COMPUTED PCA AND DATASET COLUMN
# ============================================================================

# Required Libraries
# ----------------------------------------------------------------------------
library(tidyverse)      # Data manipulation and visualization
library(ggplot2)        # Advanced plotting
library(corrplot)       # Correlation matrices
library(gridExtra)      # Multi-panel plots
library(RColorBrewer)   # Color palettes

# ============================================================================
# HELPER FUNCTIONS FOR DATASET COLUMN
# ============================================================================

get_rna_type <- function(dataset_values) {
  case_when(
    str_detect(dataset_values, "protein-coding") | str_detect(dataset_values, "protein-exon") ~ "mRNA",
    str_detect(dataset_values, "lncrna") ~ "lncRNA",
    str_detect(dataset_values, "short-ncrna") ~ "sncRNA",
    TRUE ~ NA_character_
  )
}

get_case_type <- function(dataset_values) {
  ifelse(str_detect(dataset_values, "negative-control"), "negative", "positive")
}

# ============================================================================
# 1. EXAMINE PC6 LOADINGS PRECISELY
# ============================================================================

examine_pc6_loadings <- function(pca_result, top_n = 10) {
  # Extract PC6 loadings (rotation matrix)
  pc6_loadings <- pca_result$rotation[, 6]
  
  # Sort by absolute value to find strongest contributors
  pc6_df <- data.frame(
    Feature = names(pc6_loadings),
    Loading = pc6_loadings,
    Abs_Loading = abs(pc6_loadings)
  ) %>%
    arrange(desc(Abs_Loading))
  
  # Print summary
  cat("\n", rep("=", 60), "\n")
  cat("PC6 LOADINGS ANALYSIS\n")
  cat(rep("=", 60), "\n")
  cat("\nTop", top_n, "features contributing to PC6:\n")
  print(pc6_df %>% head(top_n))
  
  # Visualize loadings
  p <- ggplot(pc6_df %>% head(top_n), 
              aes(x = reorder(Feature, Loading), y = Loading, fill = Loading > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "coral"),
                      labels = c("Positive", "Negative"),
                      name = "Direction") +
    labs(title = "PC6 Feature Loadings",
         subtitle = "Sign and magnitude of contributions",
         x = "Feature", y = "Loading Weight") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  
  # Check for cancellation patterns in key features
  key_features <- c("coding_potential", "RPKM_primary.cell", "RPKM_tissue")
  key_features_in_data <- intersect(key_features, names(pc6_loadings))
  
  if(length(key_features_in_data) > 0) {
    cat("\nKey feature loadings:\n")
    print(pc6_loadings[key_features_in_data])
    
    # Check if signs suggest cancellation
    signs <- sign(pc6_loadings[key_features_in_data])
    if(length(unique(signs)) > 1) {
      cat("\n⚠️  WARNING: Opposite signs detected - potential cancellation effect!\n")
    }
  }
  
  return(pc6_df)
}

# ============================================================================
# 2. PAIRWISE FEATURE CORRELATIONS STRATIFIED BY RNA TYPE
# ============================================================================

plot_correlation_by_rna_type <- function(df, features, dataset_col = "Dataset") {
  
  cat("\n", rep("=", 60), "\n")
  cat("CORRELATION ANALYSIS BY RNA TYPE\n")
  cat(rep("=", 60), "\n")
  
  # Add temporary RNA type column
  df_temp <- df %>%
    mutate(rna_type = get_rna_type(!!sym(dataset_col)))
  
  rna_types <- unique(df_temp$rna_type)
  rna_types <- rna_types[!is.na(rna_types)]
  
  # Create correlation comparison
  key_pairs <- data.frame(
    RNA_Type = character(),
    coding_RPKM_tissue = numeric(),
    coding_RPKM_cell = numeric(),
    RPKM_tissue_RPKM_cell = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(rna in rna_types) {
    cat("\n", rep("-", 50), "\n")
    cat("RNA Type:", rna, "\n")
    cat(rep("-", 50), "\n")
    
    # Subset data
    subset_data <- df_temp %>% 
      filter(rna_type == rna) %>%
      dplyr::select(all_of(features))
    
    # Remove NA
    subset_data <- subset_data[complete.cases(subset_data), ]
    
    if(nrow(subset_data) < 3) {
      cat("⚠️  Insufficient data for", rna, "\n")
      next
    }
    
    # Compute correlation matrix
    cor_matrix <- cor(subset_data, use = "complete.obs")
    
    # Create correlation plot
    corrplot(cor_matrix, 
             method = "color",
             type = "upper",
             tl.col = "black",
             tl.srt = 45,
             tl.cex = 0.7,
             title = paste(rna, "Correlation Matrix"),
             mar = c(0, 0, 2, 0),
             addCoef.col = "black",
             number.cex = 0.5)
    
    # Extract key correlations
    if(all(c("coding_potential", "RPKM_tissue", "RPKM_primary.cell") %in% colnames(cor_matrix))) {
      key_pairs <- rbind(key_pairs, data.frame(
        RNA_Type = rna,
        coding_RPKM_tissue = cor_matrix["coding_potential", "RPKM_tissue"],
        coding_RPKM_cell = cor_matrix["coding_potential", "RPKM_primary.cell"],
        RPKM_tissue_RPKM_cell = cor_matrix["RPKM_tissue", "RPKM_primary.cell"]
      ))
    }
  }
  
  cat("\n\nKey Feature Correlations Summary:\n")
  print(key_pairs)
  
  return(key_pairs)
}

# ============================================================================
# 3. VARIANCE EXPLAINED BY COMPONENTS
# ============================================================================

analyze_variance_explained <- function(pca_result, n_components = 10) {
  
  cat("\n", rep("=", 60), "\n")
  cat("VARIANCE EXPLAINED ANALYSIS\n")
  cat(rep("=", 60), "\n")
  
  # Calculate variance explained
  var_explained <- pca_result$sdev^2
  prop_var <- var_explained / sum(var_explained)
  cum_var <- cumsum(prop_var)
  
  # Create summary dataframe
  var_df <- data.frame(
    PC = paste0("PC", 1:length(prop_var)),
    Variance = var_explained,
    Proportion = prop_var * 100,
    Cumulative = cum_var * 100
  )
  
  # Print PC6 specifically
  cat("\nPC6 Variance Explained:\n")
  cat(sprintf("  Proportion: %.2f%%\n", var_df$Proportion[6]))
  cat(sprintf("  Cumulative through PC6: %.2f%%\n", var_df$Cumulative[6]))
  
  # Plot scree plot
  p1 <- ggplot(var_df %>% head(n_components), 
               aes(x = factor(PC, levels = PC), y = Proportion)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    geom_point(aes(y = Proportion), color = "darkblue", size = 3) +
    geom_line(aes(group = 1, y = Proportion), color = "darkblue") +
    geom_hline(yintercept = 5, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = "Scree Plot: Variance Explained by PC",
         subtitle = "Red line at 5% threshold",
         x = "Principal Component", y = "% Variance Explained") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot cumulative variance
  p2 <- ggplot(var_df %>% head(n_components), 
               aes(x = factor(PC, levels = PC), y = Cumulative, group = 1)) +
    geom_line(color = "darkgreen", size = 1) +
    geom_point(color = "darkgreen", size = 3) +
    geom_hline(yintercept = c(80, 90), linetype = "dashed", 
               color = c("orange", "red"), alpha = 0.5) +
    labs(title = "Cumulative Variance Explained",
         x = "Principal Component", y = "Cumulative % Variance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Display plots
  grid.arrange(p1, p2, ncol = 2)
  
  return(var_df)
}

# ============================================================================
# 4. VISUALIZE RAW FEATURE DISTRIBUTIONS
# ============================================================================

plot_feature_distributions <- function(df, dataset_col = "Dataset") {
  
  cat("\n", rep("=", 60), "\n")
  cat("FEATURE DISTRIBUTION ANALYSIS\n")
  cat(rep("=", 60), "\n")
  
  # Add temporary columns
  df_temp <- df %>%
    mutate(
      rna_type = get_rna_type(!!sym(dataset_col)),
      case_type = get_case_type(!!sym(dataset_col))
    )
  
  # Key feature pairs to examine
  feature_pairs <- list(
    c("coding_potential", "RPKM_tissue"),
    c("coding_potential", "RPKM_primary.cell"),
    c("RPKM_tissue", "RPKM_primary.cell")
  )
  
  plots <- list()
  idx <- 1
  
  for(pair in feature_pairs) {
    # Check if both features exist
    if(!all(pair %in% names(df_temp))) {
      cat("⚠️  Skipping pair", paste(pair, collapse = " vs "), "- features not found\n")
      next
    }
    
    for(rna in unique(df_temp$rna_type)) {
      if(is.na(rna)) next
      
      p <- df_temp %>%
        filter(rna_type == rna) %>%
        ggplot(aes(x = !!sym(pair[1]), y = !!sym(pair[2]), 
                   color = case_type)) +
        geom_point(alpha = 0.5, size = 2) +
        stat_ellipse(level = 0.95, size = 1) +
        scale_color_manual(values = c("positive" = "#E41A1C", 
                                      "negative" = "#377EB8")) +
        labs(title = paste(rna, ":", pair[1], "vs", pair[2]),
             color = "Case Type") +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      plots[[idx]] <- p
      idx <- idx + 1
    }
  }
  
  # Display all plots
  if(length(plots) > 0) {
    do.call(grid.arrange, c(plots, ncol = 3))
  }
  
  return(plots)
}


# ============================================================================
# 6. EXAMINE PC1-PC6 FOR SEPARATION
# ============================================================================

examine_pc_separation <- function(pca_result, df, pcs = 1:6, dataset_col = "Dataset") {
  
  cat("\n", rep("=", 60), "\n")
  cat("PC SEPARATION ANALYSIS\n")
  cat(rep("=", 60), "\n")
  
  # Add temporary columns
  df_temp <- df %>%
    mutate(
      rna_type = get_rna_type(!!sym(dataset_col)),
      case_type = get_case_type(!!sym(dataset_col))
    )
  
  # Add PC scores to dataframe
  pc_scores <- as.data.frame(pca_result$x[, pcs])
  df_with_pcs <- cbind(df_temp, pc_scores)
  
  # Calculate separation metrics for each PC and RNA type
  separation_stats <- expand.grid(
    PC = paste0("PC", pcs),
    RNA_Type = unique(df_temp$rna_type),
    stringsAsFactors = FALSE
  )
  
  separation_stats <- separation_stats %>% filter(!is.na(RNA_Type))
  
  separation_stats$KS_Statistic <- NA
  separation_stats$KS_PValue <- NA
  separation_stats$Mean_Diff <- NA
  
  for(i in 1:nrow(separation_stats)) {
    pc <- separation_stats$PC[i]
    rna <- separation_stats$RNA_Type[i]
    
    subset_data <- df_with_pcs %>% filter(rna_type == rna)
    
    pos_values <- subset_data %>% 
      filter(case_type == "positive") %>% 
      pull(!!sym(pc))
    
    neg_values <- subset_data %>% 
      filter(case_type == "negative") %>% 
      pull(!!sym(pc))
    
    if(length(pos_values) > 0 && length(neg_values) > 0) {
      # K-S test
      ks_test <- ks.test(pos_values, neg_values)
      separation_stats$KS_Statistic[i] <- ks_test$statistic
      separation_stats$KS_PValue[i] <- ks_test$p.value
      separation_stats$Mean_Diff[i] <- abs(mean(pos_values) - mean(neg_values))
    }
  }
  
  # Print results
  cat("\nSeparation Statistics for PC1-PC6:\n")
  print(separation_stats %>% arrange(desc(KS_Statistic)))
  
  # Visualize
  p1 <- ggplot(separation_stats, 
               aes(x = PC, y = KS_Statistic, fill = RNA_Type)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +
    labs(title = "K-S Statistics Across Principal Components",
         subtitle = "Red line at D = 0.4 threshold",
         y = "K-S Statistic (D)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p1)
  
  # Plot distributions for each PC
  for(pc_name in paste0("PC", pcs)) {
    p <- ggplot(df_with_pcs, aes(x = !!sym(pc_name), fill = case_type)) +
      geom_density(alpha = 0.5) +
      facet_wrap(~rna_type, ncol = 3, scales = "free") +
      scale_fill_manual(values = c("positive" = "#E41A1C", 
                                   "negative" = "#377EB8")) +
      labs(title = paste(pc_name, "Distribution by RNA Type"),
           x = pc_name, y = "Density") +
      theme_minimal()
    
    print(p)
  }
  
  return(separation_stats)
}


# ============================================================================
# 8. CLASS-CONDITIONAL COVARIANCE MATRICES
# ============================================================================

analyze_covariance_structure <- function(df, dataset_col = "Dataset") {
  
  cat("\n", rep("=", 60), "\n")
  cat("COVARIANCE STRUCTURE ANALYSIS\n")
  cat(rep("=", 60), "\n")
  
  # Add temporary columns
  df_temp <- df %>%
    mutate(
      rna_type = get_rna_type(!!sym(dataset_col)),
      case_type = get_case_type(!!sym(dataset_col))
    )
  
  key_features <- c("coding_potential", "RPKM_primary.cell", "RPKM_tissue")
  
  # Check which features exist
  key_features <- intersect(key_features, names(df_temp))
  
  if(length(key_features) < 2) {
    cat("⚠️  Insufficient key features found for covariance analysis\n")
    return(NULL)
  }
  
  results <- list()
  
  for(rna in unique(df_temp$rna_type)) {
    if(is.na(rna)) next
    
    cat("\n", rep("-", 50), "\n")
    cat("Covariance Analysis for", rna, "\n")
    cat(rep("-", 50), "\n")
    
    # Positive cases
    pos_data <- df_temp %>%
      filter(rna_type == rna, case_type == "positive") %>%
      dplyr::select(all_of(key_features)) %>%
      na.omit()
    
    # Negative cases
    neg_data <- df_temp %>%
      filter(rna_type == rna, case_type == "negative") %>%
      dplyr::select(all_of(key_features)) %>%
      na.omit()
    
    if(nrow(pos_data) < 3 || nrow(neg_data) < 3) {
      cat("⚠️  Insufficient data for covariance calculation\n")
      next
    }
    
    # Covariance matrices
    cov_pos <- cov(pos_data)
    cov_neg <- cov(neg_data)
    
    cat("\nPositive Cases Covariance Matrix:\n")
    print(round(cov_pos, 4))
    
    cat("\nNegative Cases Covariance Matrix:\n")
    print(round(cov_neg, 4))
    
    # Difference
    cov_diff <- cov_pos - cov_neg
    cat("\nDifference (Positive - Negative):\n")
    print(round(cov_diff, 4))
    
    # Frobenius norm of difference
    frob_norm <- sqrt(sum(cov_diff^2))
    cat("\nFrobenius Norm of Difference:", round(frob_norm, 4), "\n")
    
    # Visualize covariance matrices
    par(mfrow = c(1, 3))
    corrplot(cov2cor(cov_pos), method = "color", 
             title = paste(rna, "- Positive"),
             mar = c(0, 0, 2, 0), is.corr = FALSE)
    corrplot(cov2cor(cov_neg), method = "color", 
             title = paste(rna, "- Negative"),
             mar = c(0, 0, 2, 0), is.corr = FALSE)
    corrplot(cov2cor(abs(cov_diff)), method = "color", 
             title = paste(rna, "- |Difference|"),
             mar = c(0, 0, 2, 0), is.corr = FALSE)
    par(mfrow = c(1, 1))
    
    results[[rna]] <- list(
      cov_positive = cov_pos,
      cov_negative = cov_neg,
      difference = cov_diff,
      frobenius_norm = frob_norm
    )
  }
  
  return(results)
}

# ============================================================================
# 9. STATISTICAL TEST: CORRELATION WITH PC6 BY CASE TYPE
# ============================================================================

test_pc6_correlations <- function(pca_result, df, dataset_col = "Dataset") {
  
  cat("\n", rep("=", 60), "\n")
  cat("PC6 CORRELATION TEST\n")
  cat(rep("=", 60), "\n")
  
  # Add temporary columns
  df_temp <- df %>%
    mutate(
      rna_type = get_rna_type(!!sym(dataset_col)),
      case_type = get_case_type(!!sym(dataset_col))
    )
  
  # Add PC6 scores
  df_temp$PC6 <- pca_result$x[, 6]
  
  key_features <- c("coding_potential", "RPKM_primary.cell", "RPKM_tissue")
  key_features <- intersect(key_features, names(df_temp))
  
  results <- expand.grid(
    RNA_Type = unique(df_temp$rna_type),
    Case_Type = c("positive", "negative"),
    Feature = key_features,
    stringsAsFactors = FALSE
  )
  
  results <- results %>% filter(!is.na(RNA_Type))
  
  results$Correlation <- NA
  results$P_Value <- NA
  
  for(i in 1:nrow(results)) {
    rna <- results$RNA_Type[i]
    case <- results$Case_Type[i]
    feat <- results$Feature[i]
    
    # Get subset
    subset_data <- df_temp %>%
      filter(rna_type == rna, case_type == case) %>%
      dplyr::select(all_of(c(feat, "PC6"))) %>%
      na.omit()
    
    if(nrow(subset_data) > 3) {
      # Compute correlation
      cor_test <- cor.test(subset_data[[feat]], subset_data$PC6)
      results$Correlation[i] <- cor_test$estimate
      results$P_Value[i] <- cor_test$p.value
    }
  }
  
  # Print results
  cat("\nCorrelations of Key Features with PC6:\n")
  cat("(by RNA Type and Case Type)\n\n")
  print(results %>% 
          arrange(RNA_Type, Feature, Case_Type) %>%
          mutate(Correlation = round(Correlation, 3),
                 P_Value = format.pval(P_Value, digits = 3)))
  
  # Visualize
  p <- ggplot(results, aes(x = Feature, y = Correlation, 
                           fill = Case_Type)) +
    geom_col(position = "dodge") +
    facet_wrap(~RNA_Type, ncol = 3) +
    scale_fill_manual(values = c("positive" = "#E41A1C", 
                                 "negative" = "#377EB8")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Correlation of Key Features with PC6",
         subtitle = "Stratified by RNA Type and Case Type",
         y = "Pearson Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  # Test hypothesis
  cat("\n\nHypothesis Test Results:\n")
  for(rna in unique(df_temp$rna_type)) {
    if(is.na(rna)) next
    for(feat in key_features) {
      cor_pos <- results %>% 
        filter(RNA_Type == rna, Case_Type == "positive", Feature == feat) %>%
        pull(Correlation)
      
      cor_neg <- results %>% 
        filter(RNA_Type == rna, Case_Type == "negative", Feature == feat) %>%
        pull(Correlation)
      
      if(length(cor_pos) > 0 && length(cor_neg) > 0 && 
         !is.na(cor_pos) && !is.na(cor_neg)) {
        diff <- abs(cor_pos - cor_neg)
        cat(sprintf("%s - %s: |r_pos - r_neg| = %.3f\n", rna, feat, diff))
      }
    }
    cat("\n")
  }
  
  return(results)
}

# ============================================================================
# 10. COMPREHENSIVE SUMMARY REPORT
# ============================================================================

generate_summary_report <- function() {
  cat("\n")
  cat(rep("=", 70), "\n")
  cat("RNA PCA ANALYSIS - COMPREHENSIVE SUMMARY REPORT\n")
  cat(rep("=", 70), "\n\n")
  
  cat("ANALYSES COMPLETED:\n")
  cat("  ✓ PC6 loading analysis\n")
  cat("  ✓ Correlation structure by RNA type\n")
  cat("  ✓ Variance explained by components\n")
  cat("  ✓ Raw feature distribution visualization\n")
  cat("  ✓ PC1-PC6 separation analysis\n")
  cat("  ✓ Class-conditional covariance matrices\n")
  cat("  ✓ PC6 correlation tests\n\n")
  
  cat(rep("=", 70), "\n")
}

# ============================================================================
# COMPLETE WORKFLOW - RUN ALL ANALYSES
# ============================================================================

run_all_diagnostics <- function(df, pca_result, features = NULL) {
  
  cat("\n")
  cat(rep("=", 70), "\n")
  cat("RNA PCA DIAGNOSTIC ANALYSIS - STARTING\n")
  cat(rep("=", 70), "\n\n")
  
  # If features not provided, identify them automatically
  if(is.null(features)) {
    exclude_cols <- c("ID", "Functional", "Chromosome", "Start", "End", 
                      "Sequence", "GeneID", "Random", "Dataset")
    features <- setdiff(names(df), exclude_cols)
    
    # Convert to numeric
    for(col in features) {
      df[[col]] <- as.numeric(as.character(df[[col]]))
    }
    
    cat("Auto-detected", length(features), "features\n")
  }
  
  cat("Data:", nrow(df), "rows\n")
  cat("PCA components:", ncol(pca_result$x), "\n\n")
  
  # Run all diagnostic analyses
  cat("\n### 1. PC6 LOADINGS ANALYSIS ###\n")
  pc6_analysis <- examine_pc6_loadings(pca_result)
  
  cat("\n### 2. CORRELATION ANALYSIS ###\n")
  correlation_analysis <- plot_correlation_by_rna_type(df, features)
  
  cat("\n### 3. VARIANCE EXPLAINED ###\n")
  variance_summary <- analyze_variance_explained(pca_result)
  
  cat("\n### 4. FEATURE DISTRIBUTIONS ###\n")
  distribution_plots <- plot_feature_distributions(df)
  
  cat("\n### 6. PC SEPARATION ANALYSIS ###\n")
  pc_separation <- examine_pc_separation(pca_result, df)
  
  cat("\n### 8. COVARIANCE ANALYSIS ###\n")
  covariance_results <- analyze_covariance_structure(df)
  
  cat("\n### 9. PC6 CORRELATION TEST ###\n")
  pc6_correlation_test <- test_pc6_correlations(pca_result, df)
  
  cat("\n### 10. SUMMARY REPORT ###\n")
  generate_summary_report()
  
  # Compile all results
  all_results <- list(
    pc6_loadings = pc6_analysis,
    correlations = correlation_analysis,
    variance = variance_summary,
    pc_separation = pc_separation,
    covariance = covariance_results,
    pc6_correlations = pc6_correlation_test
  )
  
  # Save results
  saveRDS(all_results, "rna_pca_diagnostic_results.rds")
  cat("\n✓ All results saved to: rna_pca_diagnostic_results.rds\n")
  
  return(all_results)
}

to_numeric_safely <- function(x, na_vals = c("", "NA", "NaN", "nan", "NULL", "null", "N/A", "n/a"),
                              na_rate_threshold = 0.25) {
  if (is.numeric(x)) return(x)
  
  chr <- as.character(x)
  
  # Normalize whitespace and Unicode minus
  chr <- str_replace_all(chr, "\\s+", "")
  chr <- str_replace_all(chr, "\u2212", "-")  # Unicode minus → ASCII '-'
  
  # Early NA pass
  chr[chr %in% na_vals] <- NA_character_
  
  # Heuristics to detect decimal mark
  frac_dot   <- mean(str_detect(chr, "\\d+\\.\\d+"), na.rm = TRUE)
  frac_comma <- mean(str_detect(chr, "\\d+,\\d+"),    na.rm = TRUE)
  loc <- if (isTRUE(frac_comma > frac_dot)) {
    locale(decimal_mark = ",", grouping_mark = ".")
  } else {
    locale(decimal_mark = ".", grouping_mark = ",")
  }
  
  # Try strict parsing first (handles scientific notation correctly)
  num <- suppressWarnings(parse_double(chr, na = na_vals, locale = loc, trim_ws = TRUE))
  
  # If too many NAs, fall back to parse_number (tolerant: strips non-numeric chars)
  if (mean(is.na(num)) > na_rate_threshold) {
    num2 <- suppressWarnings(parse_number(chr, na = na_vals, locale = loc, trim_ws = TRUE))
    if (sum(!is.na(num2)) > sum(!is.na(num))) num <- num2
  }
  
  num
}


# ============================================================================
# RUN
# ============================================================================
# Load data
source("scripts/config.R")
source("scripts/load_gene_functionality_zscores.R")
source("scripts/load_gene_functionality_features.R")

zscores_all <- load_gene_functionality_zscores()
df <- load_gene_functionality_features()
# Creates a new numeric columns
df <- df %>%
  mutate(
    coding_potential_num = to_numeric_safely(coding_potential),
    GERP_91_mammals_max_num = to_numeric_safely(GERP_91_mammals_max),
    GERP_63_amniotes_max_num = to_numeric_safely(GERP_63_amniotes_max)
  )

# PCA
PCA_20_SELECT_FEATURES <- c("GC_percentage",
                            "CpG", "TA", "GA",
                            "phyloP_max_241w", "phyloP_max_100w",
                            "GERP_91_mammals_max_num", "GERP_63_amniotes_max_num",
                            "RPKM_tissue", "RPKM_primary.cell",
                            "copy_number", 
                            "coding_potential_num",
                            "Max_covariance", "MFE",
                            "methylome",
                            "Interaction_ave", 
                            "H3K9ac_MaxScaledSignal", "H3K79me2_MaxScaledSignal", "H3K79me1_MaxScaledSignal",
                            "chrm_acc_MaxScaledSignal")#, "repeat_distance")

PCA_20_SELECT_FEATURES_LABELS <- c("GC%",
                                   "CpG",  "GA", "TA",
                                   "PhyloP-mammals", "PhyloP-vertebrates",
                                   "GERP-mammals", "GERP-vertebrates",
                                   "Tissue RPKM", "Primary cell RPKM", 
                                   "Copies", "RNAcode",
                                   "Covariance", "MFE",
                                   "Methylome",
                                   "Interactions",
                                   "H3K9ac", "H3K79me2", "H3K79me1",
                                   "Chromatin")#, "Repeat free")
# Select all data
allrna_data_normalized <- zscores_all |>
  dplyr::select(all_of(PCA_20_SELECT_FEATURES), Dataset) |>
  na.omit()
pca_result <- prcomp(scale(allrna_data_normalized |> dplyr::select(-Dataset)))

# Run all analyses
results <- run_all_diagnostics(df |> dplyr::select(all_of(PCA_20_SELECT_FEATURES), Dataset) |> na.omit(), pca_result)

