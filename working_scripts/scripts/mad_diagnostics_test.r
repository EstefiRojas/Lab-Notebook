# ---- MAD(negative) diagnostics ----------------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(forcats)

# CONFIG: choose features to inspect (add more if you like)
features <- c("RPKM_primary.cell", "RPKM_tissue", "coding_potential_num")

# Helper: ensure RNA_type exists (fallback heuristic from Dataset)
ensure_type <- function(df) {
  if (!"RNA_type" %in% names(df)) {
    df <- df %>%
      mutate(
        RNA_type = case_when(
          str_detect(Dataset, regex("protein", ignore_case = TRUE)) ~ "mRNA",
          str_detect(Dataset, regex("lnc", ignore_case = TRUE)) ~ "lncRNA",
          str_detect(Dataset, regex("short[-_]?ncrna|sncrna|small", ignore_case = TRUE)) ~ "sncRNA",
          TRUE ~ "other"
        )
      )
  }
  df
}

# Helper: compute negative-only robust scale with fallback
neg_scale_stats <- function(x, tiny = 1e-6) {
  # R's mad() returns 1.4826 * median(|x - median|) by default
  mad_neg <- suppressWarnings(mad(x, na.rm = TRUE))
  iqr_neg <- IQR(x, na.rm = TRUE)
  fb      <- iqr_neg / 1.349
  use_fb  <- !is.finite(mad_neg) || mad_neg < tiny
  used    <- if (use_fb) fb else mad_neg
  list(
    mad_neg = mad_neg,
    iqr_neg = iqr_neg,
    fallback = fb,
    used_scale = used,
    used_scale_source = if (use_fb) "IQR/1.349" else "MAD",
    tiny_flag = use_fb
  )
}

# Main diagnostics function
mad_negative_diagnostics <- function(df, features,
                                     dataset_col = "Dataset",
                                     type_col = "RNA_type",
                                     neg_regex = "negative-control",
                                     tiny = 1e-6,
                                     make_plots = TRUE) {
  stopifnot(all(features %in% names(df)))
  df <- ensure_type(df)
  
  is_neg <- str_detect(df[[dataset_col]], regex(neg_regex, ignore_case = TRUE))
  df <- df %>% mutate(.is_neg = is_neg)
  
  # Long format for raw & log1p variants
  long_raw <- df %>%
    dplyr::select(all_of(c(dataset_col, type_col, ".is_neg", features))) %>%
    pivot_longer(cols = all_of(features), names_to = "feature", values_to = "value") %>%
    mutate(scale_variant = "raw")
  
  long_log <- long_raw %>%
    mutate(value = log1p(value), scale_variant = "log1p")
  
  long_all <- bind_rows(long_raw, long_log)
  
  # Compute negative-only scale stats per type x feature x variant
  diag_tbl <- long_all %>%
    group_by(.data[[type_col]], feature, scale_variant) %>%
    summarise(
      n_neg     = sum(.is_neg & is.finite(value)),
      n_pos     = sum((!.is_neg) & is.finite(value)),
      med_neg   = median(value[.is_neg], na.rm = TRUE),
      prop_zero_neg = mean(replace_na(value[.is_neg] == 0, FALSE)),
      # pack list column with robust scales
      stats = list(neg_scale_stats(value[.is_neg], tiny = tiny)),
      .groups = "drop"
    ) %>%
    mutate(
      mad_neg   = map_dbl(stats, "mad_neg"),
      iqr_neg   = map_dbl(stats, "iqr_neg"),
      fallback  = map_dbl(stats, "fallback"),
      used_scale = map_dbl(stats, "used_scale"),
      used_scale_source = map_chr(stats, "used_scale_source"),
      tiny_flag = map_lgl(stats, "tiny_flag")
    ) %>%
    dplyr::select(-stats) %>%
    arrange(scale_variant, .data[[type_col]], feature)
  
  # Add a rough standardized effect size (median_pos - median_neg)/used_scale
  eff_tbl <- long_all %>%
    group_by(.data[[type_col]], feature, scale_variant) %>%
    summarise(
      med_pos = median(value[!.is_neg], na.rm = TRUE),
      med_neg = median(value[ .is_neg], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(diag_tbl %>%
                dplyr::select(.data[[type_col]], feature, scale_variant, used_scale),
              by = c(type_col, "feature", "scale_variant")) %>%
    mutate(effect_size = (med_pos - med_neg) / (used_scale + 1e-12))
  
  diag_tbl <- diag_tbl %>%
    left_join(eff_tbl %>%
                dplyr::select(.data[[type_col]], feature, scale_variant, effect_size),
              by = c(type_col, "feature", "scale_variant"))
  
  if (make_plots) {
    # Plot 1: Used scale (MAD vs fallback) per feature x type
    p1 <- diag_tbl %>%
      mutate(feature = fct_reorder(feature, used_scale, .fun = max)) %>%
      ggplot(aes(x = feature, y = used_scale, fill = used_scale_source)) +
      geom_col(position = "dodge") +
      facet_grid(scale_variant ~ .data[[type_col]], scales = "free_y") +
      coord_flip() +
      labs(title = "Negative-only scale used (MAD vs IQR fallback)",
           y = "Scale (units of feature)", x = "Feature", fill = "Scale source")
    
    print(p1)
    
    # Plot 2: Proportion of zeros in negatives (heatmap)
    p2 <- diag_tbl %>%
      ggplot(aes(x = feature, y = .data[[type_col]], fill = prop_zero_neg)) +
      geom_tile() +
      facet_wrap(~ scale_variant) +
      coord_equal() +
      labs(title = "Proportion of exact zeros among negatives",
           x = "Feature", y = "RNA type", fill = "Prop zero") +
      scale_fill_continuous(limits = c(0,1)) +
      theme(
        axis.text.x = element_text(angle = 90)
      )
    print(p2)
    
    # Plot 3: Standardized effect size (median_pos - median_neg)/used_scale
    p3 <- diag_tbl %>%
      mutate(feature = fct_reorder(feature, effect_size, .fun = max)) %>%
      ggplot(aes(x = feature, y = effect_size)) +
      geom_col() +
      facet_grid(scale_variant ~ .data[[type_col]], scales = "free_y") +
      coord_flip() +
      labs(title = "Standardized median shift using negative-only scale",
           y = "Effect size (Δ median / used scale)", x = "Feature")
    print(p3)
  }
  
  # Brief console warnings for tiny/zero MAD
  tiny_hits <- diag_tbl %>% filter(tiny_flag)
  if (nrow(tiny_hits) > 0) {
    message("⚠️  Tiny/zero MAD detected; falling back to IQR/1.349 for:\n",
            paste0(" - ",
                   tiny_hits[[type_col]], " | ", tiny_hits$feature,
                   " [", tiny_hits$scale_variant, "]",
                   collapse = "\n"))
  }
  
  diag_tbl %>%
    arrange(scale_variant, .data[[type_col]], desc(used_scale))
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

# ---- Apply to your df ------------------------------------------------------
# Creates a new numeric column and shows a quick diagnostic of unparsed values
df <- df %>%
  mutate(
    coding_potential_num = to_numeric_safely(coding_potential),
    GERP_91_mammals_max_num = to_numeric_safely(GERP_91_mammals_max),
    GERP_63_amniotes_max_num = to_numeric_safely(GERP_63_amniotes_max)
  )

# Diagnostics: what didn't parse?
unparsed <- df %>%
  mutate(orig = as.character(coding_potential)) %>%
  filter(!is.na(orig) & is.na(coding_potential_num)) %>%
  count(orig, sort = TRUE)

if (nrow(unparsed) > 0) {
  message("⚠️ Some coding_potential entries could not be parsed. Top examples:")
  print(head(unparsed, 20))
}

# Quick sanity check
summary(df$coding_potential_num)


# ---- RUN IT ---------------------------------------------------------------
mad_diag <- mad_negative_diagnostics(
  df,
  features = PCA_20_SELECT_FEATURES,
  dataset_col = "Dataset",
  type_col = "RNA_type",     # if you don't have this, it will be created heuristically
  neg_regex = "negative-control",
  tiny = 1e-6,
  make_plots = TRUE
)

# Look at the table
print(mad_diag, n = 60)

