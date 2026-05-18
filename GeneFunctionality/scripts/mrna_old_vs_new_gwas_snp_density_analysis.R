# ===============================================
# GWAS SNP Density Analysis for mrna_old_vs_new (v7 vs v49)
# ===============================================
# Forked from mrna_gwas_snp_density_analysis.R for the 500-gene
# v7-deprecated vs v49-current paired mRNA dataset.

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

custom_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "red", "yellow", "green")

# Input: output of extract_gwas_data_mrna.sh
gwas_data <- read.csv("../results/mrna_old_vs_new_gwas_snp_density.csv", header = TRUE)

# Highest probability across exons
gwas_data <- gwas_data %>%
  mutate(
    highest_prob = pmax(ex2_prob_Yes, ex3_prob_Yes, na.rm = TRUE)
  )

cat("Dataset dimensions:", nrow(gwas_data), "genes x", ncol(gwas_data), "columns\n")
cat("Genes with transcript SNPs:", sum(gwas_data$transcript_SNP_count_ct > 0, na.rm = TRUE), "\n")
cat("Functional split:\n")
print(table(gwas_data$Functional))

#####################
# SNP Count Analysis #
#####################

filtered_snp_count_data <- gwas_data %>%
  mutate(
    selected_snp_count = ifelse(ex2_prob_Yes == ex3_prob_Yes,
                                ifelse(is.na(ex3_SNP_count_ct),
                                       ex2_SNP_count_ct,
                                       ifelse(!is.na(ex2_SNP_count_ct),
                                              ifelse(ex2_SNP_count_ct >= ex3_SNP_count_ct,
                                                     ex2_SNP_count_ct,
                                                     ex3_SNP_count_ct),
                                              ex3_SNP_count_ct)),
                                ifelse(ex2_prob_Yes > ex3_prob_Yes,
                                       ex2_SNP_count_ct,
                                       ex3_SNP_count_ct)),
    max_snp_count = transcript_SNP_count_ct,

    snp_count_category = case_when(
      transcript_SNP_count_ct < 1 ~ "0",
      transcript_SNP_count_ct >= 1 & transcript_SNP_count_ct <= 10 ~ "1-10",
      transcript_SNP_count_ct > 10 & transcript_SNP_count_ct <= 50 ~ "11-50",
      transcript_SNP_count_ct > 50 ~ ">50"
    ),

    snp_density_per_kb = transcript_SNP_count_ct / (End_Transcript - Start_Transcript) * 1000
  )


density_quintiles <- quantile(
  filtered_snp_count_data$snp_density_per_kb[filtered_snp_count_data$snp_density_per_kb > 0],
  probs = seq(0, 1, by = 0.25), na.rm = TRUE
)
cat("\n=== SNP Density Quartiles (non-zero only) ===\n")
print(density_quintiles)

q1 <- round(density_quintiles[2], 3)
q2 <- round(density_quintiles[3], 3)
q3 <- round(density_quintiles[4], 3)

filtered_snp_count_data <- filtered_snp_count_data %>%
  mutate(
    snp_density_category = case_when(
      is.na(snp_density_per_kb) ~ "NA",
      snp_density_per_kb == 0 ~ "0",
      snp_density_per_kb > 0 & snp_density_per_kb <= density_quintiles[2] ~ "Q1",
      snp_density_per_kb > density_quintiles[2] & snp_density_per_kb <= density_quintiles[3] ~ "Q2",
      snp_density_per_kb > density_quintiles[3] & snp_density_per_kb <= density_quintiles[4] ~ "Q3",
      snp_density_per_kb > density_quintiles[4] ~ "Q4"
    )
  )

count_quintiles <- quantile(
  filtered_snp_count_data$transcript_SNP_count_ct[filtered_snp_count_data$transcript_SNP_count_ct > 0],
  probs = seq(0, 1, by = 0.25), na.rm = TRUE
)
cat("\n=== SNP Count Quartiles (non-zero only) ===\n")
print(count_quintiles)

filtered_snp_count_data <- filtered_snp_count_data %>%
  mutate(
    snp_count_category = case_when(
      transcript_SNP_count_ct == 0 ~ "0",
      transcript_SNP_count_ct >= 1 & transcript_SNP_count_ct <= count_quintiles[2] ~ "Q1",
      transcript_SNP_count_ct > count_quintiles[2] & transcript_SNP_count_ct <= count_quintiles[3] ~ "Q2",
      transcript_SNP_count_ct > count_quintiles[3] & transcript_SNP_count_ct <= count_quintiles[4] ~ "Q3",
      transcript_SNP_count_ct > count_quintiles[4] ~ "Q4"
    )
  )

cat("\n=== SNP Count Summary ===\n")
print(summary(filtered_snp_count_data$transcript_SNP_count_ct))
cat("\n=== SNP Density per kb Summary ===\n")
print(summary(filtered_snp_count_data$snp_density_per_kb))
cat("\n=== SNP Count Categories ===\n")
print(table(filtered_snp_count_data$snp_count_category))
cat("\n=== SNP Density Categories ===\n")
print(table(filtered_snp_count_data$snp_density_category))

filtered_snp_count_data$snp_count_category <- factor(
  filtered_snp_count_data$snp_count_category,
  levels = c("0", "Q1", "Q2", "Q3", "Q4")
)

filtered_snp_count_data$snp_density_category <- factor(
  filtered_snp_count_data$snp_density_category,
  levels = c("0", "Q1", "Q2", "Q3", "Q4")
)

#####################
# SNP Count Plots   #
#####################

p_hist_log <- ggplot(filtered_snp_count_data %>% filter(transcript_SNP_count_ct > 0),
                     aes(x = transcript_SNP_count_ct)) +
  geom_histogram(bins = 30) +
  scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c(1, 10, 100, 1000)) +
  labs(
    title = "GWAS SNP Count per mRNA Gene (v7+v49)",
    x = "SNP Count (log scale)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave("../results/gwas/mrna_old_vs_new_SNP_count_histogram_log.png", p_hist_log, width = 10, height = 7)

p_hist_lin <- ggplot(filtered_snp_count_data, aes(x = transcript_SNP_count_ct)) +
  geom_histogram(bins = 30) +
  labs(
    title = "GWAS SNP Count per mRNA Gene (v7+v49)",
    x = "SNP Count (linear scale)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave("../results/gwas/mrna_old_vs_new_SNP_count_histogram_linear.png", p_hist_lin, width = 10, height = 7)

range_label <- function(cat, quintiles, digits = 3) {
  if (cat == "0") return("0")
  i  <- as.integer(sub("Q", "", cat))
  lo <- round(quintiles[i], digits)
  sprintf(">%s", lo)
}

# Bin highest_prob using fixed biological bins anchored on the 0.4 model threshold.
# Quantile-based binning fails here because the model output is heavily right-skewed
# (>75% of predictions cluster at prob_Yes == 1.000), so quartile cuts collapse.
prob_breaks <- c(0, 0.4, 0.7, 0.95, 1.0001)  # 1.0001 to include the inclusive upper bound
prob_labels <- c("Q1", "Q2", "Q3", "Q4")
# Bin boundaries shown to the user in plot labels
prob_quintiles <- c(0, 0.4, 0.7, 0.95, 1.0)
cat("\n=== Functional Probability Bins (fixed) ===\n")
print(setNames(prob_quintiles, c("min", "Q1/Q2", "Q2/Q3", "Q3/Q4", "max")))

filtered_snp_count_data <- filtered_snp_count_data %>%
  mutate(
    prob_category = cut(highest_prob, breaks = prob_breaks,
                        labels = prob_labels, include.lowest = TRUE, right = FALSE)
  )
filtered_snp_count_data$prob_category <- factor(
  filtered_snp_count_data$prob_category,
  levels = prob_labels
)
cat("\n=== Functional Probability Bin Sizes ===\n")
print(table(filtered_snp_count_data$prob_category, useNA = "ifany"))

legend_data <- filtered_snp_count_data %>%
  filter(!is.na(prob_category)) %>%
  count(prob_category)

# x-axis labels: half-open intervals matching cut(..., right=FALSE).
# The final bin is closed on both ends since prob == 1.0 is in scope.
prob_label <- function(cat, quintiles, digits = 2) {
  i  <- as.integer(sub("Q", "", cat))
  lo <- round(quintiles[i], digits)
  hi <- round(quintiles[i + 1], digits)
  closer <- if (i == length(quintiles) - 1) "]" else ")"
  sprintf("[%s, %s%s", lo, hi, closer)
}
new_labels <- setNames(
  vapply(as.character(legend_data$prob_category), function(cat) {
    n_val <- legend_data$n[as.character(legend_data$prob_category) == cat]
    paste0(prob_label(cat, prob_quintiles, digits = 2), "\n(n=", n_val, ")")
  }, character(1)),
  legend_data$prob_category
)

# KS tests: compare each higher probability quartile against Q1 (lowest)
reference_group <- "Q1"
comparison_groups <- c("Q2", "Q3", "Q4")

reference_data <- filtered_snp_count_data %>%
  filter(prob_category == reference_group) %>%
  pull(transcript_SNP_count_ct)

stats_list <- lapply(comparison_groups, function(group) {
  comparison_data <- filtered_snp_count_data %>%
    filter(prob_category == group) %>%
    pull(transcript_SNP_count_ct)

  if (length(comparison_data) < 2) {
    return(data.frame(prob_category = group, label = "n<2"))
  }

  ks_result <- ks.test(reference_data, comparison_data)
  data.frame(
    prob_category = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

stats_labels <- do.call(rbind, stats_list)
# Position KS labels just above the 95th-percentile of SNP count for headroom
count_ymax <- quantile(filtered_snp_count_data$transcript_SNP_count_ct, 0.95, na.rm = TRUE)
stats_labels$y_position <- count_ymax * 1.05

p_snp_count <- ggplot(filtered_snp_count_data %>% filter(!is.na(prob_category)),
                      aes(x = prob_category, y = transcript_SNP_count_ct, fill = prob_category)) +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  scale_fill_manual(values = custom_colors[1:length(levels(filtered_snp_count_data$prob_category))]) +
  scale_x_discrete(labels = new_labels) +
  geom_text(
    data = stats_labels,
    aes(x = prob_category, y = y_position, label = label),
    inherit.aes = FALSE,
    size = 5,
    color = "black",
    fontface = "bold"
  ) +
  coord_cartesian(ylim = c(0, count_ymax * 1.15), clip = "off") +
  labs(
    title = "GWAS SNP Count (v7+v49)",
    x = "Functional Probability",
    y = "Transcript SNP Count",
    caption = "*p-val < 5e-8; KS vs Q1"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28),
    axis.text.x = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    plot.caption = element_text(size = 16, hjust = -0.09, face = "bold.italic", color = "grey40"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave("../results/gwas/mrna_old_vs_new_SNP_count_boxPlot.png", p_snp_count, width = 12, height = 8)

#####################
# SNP Density Plots #
#####################

density_data <- filtered_snp_count_data %>% filter(!is.na(snp_density_category))

p_dens_hist_log <- ggplot(density_data %>% filter(snp_density_per_kb > 0),
                          aes(x = snp_density_per_kb)) +
  geom_histogram(bins = 30) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10), labels = c(0.01, 0.1, 1, 10)) +
  labs(
    title = "GWAS SNP Density per mRNA Gene (v7+v49)",
    x = "SNP Density per kb (log scale)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave("../results/gwas/mrna_old_vs_new_SNP_density_histogram_log.png", p_dens_hist_log, width = 10, height = 7)

p_dens_hist_lin <- ggplot(density_data, aes(x = snp_density_per_kb)) +
  geom_histogram(bins = 30) +
  labs(
    title = "GWAS SNP Density per mRNA Gene (v7+v49)",
    x = "SNP Density per kb (linear scale)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave("../results/gwas/mrna_old_vs_new_SNP_density_histogram_linear.png", p_dens_hist_lin, width = 10, height = 7)

# Reuse prob_category from above; legend/labels show probability quartile ranges + n
legend_data <- density_data %>%
  filter(!is.na(prob_category)) %>%
  count(prob_category)

new_labels <- setNames(
  vapply(as.character(legend_data$prob_category), function(cat) {
    n_val <- legend_data$n[as.character(legend_data$prob_category) == cat]
    paste0(prob_label(cat, prob_quintiles, digits = 2), "\n(n=", n_val, ")")
  }, character(1)),
  legend_data$prob_category
)

# KS tests: compare each higher-probability quartile's SNP density against Q1
reference_group <- "Q1"
comparison_groups <- c("Q2", "Q3", "Q4")

reference_data <- density_data %>%
  filter(prob_category == reference_group) %>%
  pull(snp_density_per_kb)

stats_list <- lapply(comparison_groups, function(group) {
  comparison_data <- density_data %>%
    filter(prob_category == group) %>%
    pull(snp_density_per_kb)

  if (length(comparison_data) < 2) {
    return(data.frame(prob_category = group, label = "n<2"))
  }

  ks_result <- ks.test(reference_data, comparison_data)
  data.frame(
    prob_category = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

stats_labels <- do.call(rbind, stats_list)
# Position KS labels above the 95th-percentile density for headroom
dens_ymax <- quantile(density_data$snp_density_per_kb, 0.95, na.rm = TRUE)
stats_labels$y_position <- dens_ymax * 1.05

p_snp_density <- ggplot(density_data %>% filter(!is.na(prob_category)),
                        aes(x = prob_category, y = snp_density_per_kb, fill = prob_category)) +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  scale_fill_manual(values = custom_colors[1:length(levels(density_data$prob_category))]) +
  scale_x_discrete(labels = new_labels) +
  geom_text(
    data = stats_labels,
    aes(x = prob_category, y = y_position, label = label),
    inherit.aes = FALSE,
    size = 5,
    color = "black",
    fontface = "bold"
  ) +
  coord_cartesian(ylim = c(0, dens_ymax * 1.15), clip = "off") +
  labs(
    title = "GWAS SNP Density per kb (v7+v49)",
    x = "Functional Probability",
    y = "SNP Density (per kb)",
    caption = "*p-val < 5e-8; KS vs Q1"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28),
    axis.text.x = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    plot.caption = element_text(size = 16, hjust = -0.09, face = "bold.italic", color = "grey40"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave("../results/gwas/mrna_old_vs_new_SNP_density_boxPlot.png", p_snp_density, width = 12, height = 8)

#####################
# Binary probability split: below vs above the 0.4 functional threshold
#####################

binary_data <- filtered_snp_count_data %>%
  filter(!is.na(highest_prob)) %>%
  mutate(prob_binary = factor(ifelse(highest_prob < 0.4, "[0, 0.4)", "[0.4, 1.0]"),
                              levels = c("[0, 0.4)", "[0.4, 1.0]")))

binary_legend <- binary_data %>% count(prob_binary)
binary_labels <- setNames(
  paste0(binary_legend$prob_binary, "\n(n=", binary_legend$n, ")"),
  as.character(binary_legend$prob_binary)
)

# KS tests across the binary split
ref_bin_count   <- binary_data %>% filter(prob_binary == "[0, 0.4)") %>% pull(transcript_SNP_count_ct)
cmp_bin_count   <- binary_data %>% filter(prob_binary == "[0.4, 1.0]") %>% pull(transcript_SNP_count_ct)
ks_bin_count    <- ks.test(ref_bin_count, cmp_bin_count)

bin_density_data <- binary_data %>% filter(!is.na(snp_density_per_kb))
ref_bin_dens   <- bin_density_data %>% filter(prob_binary == "[0, 0.4)") %>% pull(snp_density_per_kb)
cmp_bin_dens   <- bin_density_data %>% filter(prob_binary == "[0.4, 1.0]") %>% pull(snp_density_per_kb)
ks_bin_density <- ks.test(ref_bin_dens, cmp_bin_dens)

cat("\n=== Binary Probability Split (threshold 0.4) ===\n")
print(binary_data %>%
  group_by(prob_binary) %>%
  summarise(
    n = n(),
    median_count = median(transcript_SNP_count_ct, na.rm = TRUE),
    mean_count = round(mean(transcript_SNP_count_ct, na.rm = TRUE), 2),
    median_density = round(median(snp_density_per_kb, na.rm = TRUE), 3),
    mean_density = round(mean(snp_density_per_kb, na.rm = TRUE), 3),
    .groups = "drop"
  ))
cat(sprintf("KS (count):   D = %.3f  p = %.3g\n", ks_bin_count$statistic, ks_bin_count$p.value))
cat(sprintf("KS (density): D = %.3f  p = %.3g\n", ks_bin_density$statistic, ks_bin_density$p.value))

binary_colors <- c("[0, 0.4)" = "#FC8D62", "[0.4, 1.0]" = "#66C2A5")

# --- Count boxplot (binary) ---
count_ymax_bin <- quantile(binary_data$transcript_SNP_count_ct, 0.95, na.rm = TRUE)
p_count_bin <- ggplot(binary_data,
                      aes(x = prob_binary, y = transcript_SNP_count_ct, fill = prob_binary)) +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  scale_fill_manual(values = binary_colors) +
  scale_x_discrete(labels = binary_labels) +
  annotate("text", x = 1.5, y = count_ymax_bin * 1.05,
           label = sprintf("KS = %.2f   p = %.2g", ks_bin_count$statistic, ks_bin_count$p.value),
           size = 5, fontface = "bold") +
  coord_cartesian(ylim = c(0, count_ymax_bin * 1.15), clip = "off") +
  labs(
    title = "GWAS SNP Count by Functional Probability",
    x = "Functional Probability",
    y = "Transcript SNP Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 24),
    axis.text.x = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave("../results/gwas/mrna_old_vs_new_SNP_count_boxPlot_binary.png",
       p_count_bin, width = 10, height = 8)

# --- Density boxplots (binary): three versions — v7, v49, and v7+v49 pooled ---
density_binary_specs <- list(
  list(suffix = "_v7",  subset = "v7",     title = "GWAS SNP Density (v7)"),
  list(suffix = "_v49", subset = "v49",    title = "GWAS SNP Density (v49)"),
  list(suffix = "",     subset = "v7+v49", title = "GWAS SNP Density (v7+v49)")
)

for (dbs in density_binary_specs) {
  if (dbs$subset == "v7") {
    sub_df <- bin_density_data %>% filter(Functional == "No")
  } else if (dbs$subset == "v49") {
    sub_df <- bin_density_data %>% filter(Functional == "Yes")
  } else {
    sub_df <- bin_density_data
  }

  # Skip if a subset would be empty
  if (nrow(sub_df) < 2) next

  # KS within this subset
  ref_d <- sub_df %>% filter(prob_binary == "[0, 0.4)")  %>% pull(snp_density_per_kb)
  cmp_d <- sub_df %>% filter(prob_binary == "[0.4, 1.0]") %>% pull(snp_density_per_kb)
  if (length(ref_d) >= 2 && length(cmp_d) >= 2) {
    ks_d <- ks.test(ref_d, cmp_d)
    ks_label <- sprintf("KS = %.2f   p = %.2g", ks_d$statistic, ks_d$p.value)
  } else {
    ks_label <- "n<2"
  }

  sub_legend <- sub_df %>% count(prob_binary)
  sub_labels <- setNames(
    paste0(sub_legend$prob_binary, "\n(n=", sub_legend$n, ")"),
    as.character(sub_legend$prob_binary)
  )

  dens_ymax <- quantile(sub_df$snp_density_per_kb, 0.95, na.rm = TRUE)

  p <- ggplot(sub_df, aes(x = prob_binary, y = snp_density_per_kb, fill = prob_binary)) +
    geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
    scale_fill_manual(values = binary_colors) +
    scale_x_discrete(labels = sub_labels) +
    annotate("text", x = 1.5, y = dens_ymax * 1.05, label = ks_label,
             size = 5, fontface = "bold") +
    coord_cartesian(ylim = c(0, dens_ymax * 1.15), clip = "off") +
    labs(
      title = dbs$title,
      x = "Functional Probability",
      y = "SNP Density (per kb)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 24),
      axis.text.x = element_text(hjust = 0.5),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )

  ggsave(sprintf("../results/gwas/mrna_old_vs_new_SNP_density_boxPlot_binary%s.png", dbs$suffix),
         p, width = 10, height = 8)
}

#####################
# SNP density vs sequence length scatter
#####################

# Transcript length is the natural "sequence length" here since SNP density is
# normalized by it. Drop rows without resolved transcript coords or zero span.
scatter_data <- filtered_snp_count_data %>%
  filter(!is.na(End_Transcript), !is.na(Start_Transcript)) %>%
  mutate(transcript_length = End_Transcript - Start_Transcript) %>%
  filter(transcript_length > 0, !is.na(snp_density_per_kb))

n_total <- nrow(scatter_data)
n_nonzero <- sum(scatter_data$snp_density_per_kb > 0)

# Correlations (linear scale and log-log)
pearson_lin <- cor(scatter_data$transcript_length, scatter_data$snp_density_per_kb,
                   method = "pearson", use = "pairwise.complete.obs")
spearman    <- cor(scatter_data$transcript_length, scatter_data$snp_density_per_kb,
                   method = "spearman", use = "pairwise.complete.obs")

# log-log Pearson restricted to nonzero density (log(0) undefined)
nz <- scatter_data %>% filter(snp_density_per_kb > 0)
pearson_log <- cor(log10(nz$transcript_length), log10(nz$snp_density_per_kb),
                   method = "pearson")

cat("\n=== SNP Density vs Transcript Length ===\n")
cat(sprintf("n = %d (%d with density > 0)\n", n_total, n_nonzero))
cat(sprintf("Pearson  (linear):   r = %.3f\n", pearson_lin))
cat(sprintf("Spearman (rank):     rho = %.3f\n", spearman))
cat(sprintf("Pearson  (log-log, density>0): r = %.3f\n", pearson_log))

# Linear-axes scatter
p_scatter_lin <- ggplot(scatter_data,
                        aes(x = transcript_length, y = snp_density_per_kb, color = Functional)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "grey20", linewidth = 0.7) +
  scale_color_manual(values = c("No" = "#FC8D62", "Yes" = "#66C2A5"),
                     labels = c("No" = "v7", "Yes" = "v49")) +
  labs(
    title = "SNP density vs transcript length",
    subtitle = sprintf("Pearson r = %.3f   Spearman rho = %.3f   n = %d",
                       pearson_lin, spearman, n_total),
    x = "Transcript length (bp)",
    y = "GWAS SNP density (per kb)",
    color = "Source"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    text = element_text(size = 18),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave("../results/gwas/mrna_old_vs_new_density_vs_length_linear.png",
       p_scatter_lin, width = 11, height = 7)

# Log-log scatter (drops density == 0 points which can't be log-transformed)
p_scatter_log <- ggplot(nz,
                        aes(x = transcript_length, y = snp_density_per_kb, color = Functional)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "grey20", linewidth = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("No" = "#FC8D62", "Yes" = "#66C2A5"),
                     labels = c("No" = "v7", "Yes" = "v49")) +
  labs(
    title = "SNP density vs transcript length (log-log)",
    subtitle = sprintf("Pearson r = %.3f (log-log, density>0)   n = %d", pearson_log, nrow(nz)),
    x = "Transcript length (bp, log scale)",
    y = "GWAS SNP density (per kb, log scale)",
    color = "Source"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    text = element_text(size = 18),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave("../results/gwas/mrna_old_vs_new_density_vs_length_loglog.png",
       p_scatter_log, width = 11, height = 7)

#####################
# Per-GENCODE-version boxplots: separate files for v7 and v49.
# Reuses the probability-bin x-axis from the combined boxplots, but
# subsets each plot to a single GENCODE version (v7 = deprecated genes,
# v49 = current genes).
#####################

version_specs <- list(
  list(suffix = "v7",  func_value = "No",  label = "v7",  fill = "#FC8D62"),
  list(suffix = "v49", func_value = "Yes", label = "v49", fill = "#66C2A5")
)

cat("\n=== Per-GENCODE-version per-probability-bin breakdown ===\n")

for (spec in version_specs) {
  sub_count   <- filtered_snp_count_data %>%
    filter(Functional == spec$func_value, !is.na(prob_category))
  sub_density <- sub_count %>% filter(!is.na(snp_density_per_kb))

  cat(sprintf("\n-- %s --\n", spec$label))
  print(sub_count %>%
    group_by(prob_category) %>%
    summarise(
      n = n(),
      median_count = median(transcript_SNP_count_ct, na.rm = TRUE),
      median_density = round(median(snp_density_per_kb, na.rm = TRUE), 3),
      .groups = "drop"
    ))

  # --- KS tests for SNP count, each prob bin vs Q1 (within this version) ---
  ref_count <- sub_count %>% filter(prob_category == "Q1") %>% pull(transcript_SNP_count_ct)
  stats_count <- lapply(c("Q2", "Q3", "Q4"), function(g) {
    d <- sub_count %>% filter(prob_category == g) %>% pull(transcript_SNP_count_ct)
    if (length(d) < 2 || length(ref_count) < 2) return(data.frame(prob_category = g, label = "n<2"))
    ks <- ks.test(ref_count, d)
    data.frame(prob_category = g, label = paste0("KS=", round(ks$statistic, 2)))
  })
  stats_count <- do.call(rbind, stats_count)
  count_ymax <- quantile(sub_count$transcript_SNP_count_ct, 0.95, na.rm = TRUE)
  stats_count$y_position <- count_ymax * 1.05

  # x-axis labels per bin: range + n for this subset
  count_legend <- sub_count %>% count(prob_category)
  count_labels <- setNames(
    vapply(as.character(count_legend$prob_category), function(cat) {
      n_val <- count_legend$n[as.character(count_legend$prob_category) == cat]
      paste0(prob_label(cat, prob_quintiles, digits = 2), "\n(n=", n_val, ")")
    }, character(1)),
    count_legend$prob_category
  )

  p_count <- ggplot(sub_count,
                    aes(x = prob_category, y = transcript_SNP_count_ct)) +
    geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black",
                 fill = spec$fill, staplewidth = 0.5) +
    scale_x_discrete(labels = count_labels) +
    geom_text(
      data = stats_count,
      aes(x = prob_category, y = y_position, label = label),
      inherit.aes = FALSE,
      size = 5, color = "black", fontface = "bold"
    ) +
    coord_cartesian(ylim = c(0, count_ymax * 1.15), clip = "off") +
    labs(
      title = sprintf("GWAS SNP Count (%s)", spec$label),
      x = "Functional Probability",
      y = "Transcript SNP Count",
      caption = "*p-val < 5e-8; KS vs Q1"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 24),
      axis.text.x = element_text(hjust = 0.5),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      plot.caption = element_text(size = 14, hjust = -0.09, face = "bold.italic", color = "grey40"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )

  ggsave(sprintf("../results/gwas/mrna_old_vs_new_SNP_count_boxPlot_%s.png", spec$suffix),
         p_count, width = 12, height = 8)

  # --- KS tests for SNP density, each prob bin vs Q1 ---
  ref_dens <- sub_density %>% filter(prob_category == "Q1") %>% pull(snp_density_per_kb)
  stats_dens <- lapply(c("Q2", "Q3", "Q4"), function(g) {
    d <- sub_density %>% filter(prob_category == g) %>% pull(snp_density_per_kb)
    if (length(d) < 2 || length(ref_dens) < 2) return(data.frame(prob_category = g, label = "n<2"))
    ks <- ks.test(ref_dens, d)
    data.frame(prob_category = g, label = paste0("KS=", round(ks$statistic, 2)))
  })
  stats_dens <- do.call(rbind, stats_dens)
  dens_ymax <- quantile(sub_density$snp_density_per_kb, 0.95, na.rm = TRUE)
  stats_dens$y_position <- dens_ymax * 1.05

  dens_legend <- sub_density %>% count(prob_category)
  dens_labels <- setNames(
    vapply(as.character(dens_legend$prob_category), function(cat) {
      n_val <- dens_legend$n[as.character(dens_legend$prob_category) == cat]
      paste0(prob_label(cat, prob_quintiles, digits = 2), "\n(n=", n_val, ")")
    }, character(1)),
    dens_legend$prob_category
  )

  p_dens <- ggplot(sub_density,
                   aes(x = prob_category, y = snp_density_per_kb)) +
    geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black",
                 fill = spec$fill, staplewidth = 0.5) +
    scale_x_discrete(labels = dens_labels) +
    geom_text(
      data = stats_dens,
      aes(x = prob_category, y = y_position, label = label),
      inherit.aes = FALSE,
      size = 5, color = "black", fontface = "bold"
    ) +
    coord_cartesian(ylim = c(0, dens_ymax * 1.15), clip = "off") +
    labs(
      title = sprintf("GWAS SNP Density (%s)", spec$label),
      x = "Functional Probability",
      y = "SNP Density (per kb)",
      caption = "*p-val < 5e-8; KS vs Q1"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 24),
      axis.text.x = element_text(hjust = 0.5),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      plot.caption = element_text(size = 14, hjust = -0.09, face = "bold.italic", color = "grey40"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )

  ggsave(sprintf("../results/gwas/mrna_old_vs_new_SNP_density_boxPlot_%s.png", spec$suffix),
         p_dens, width = 12, height = 8)
}

cat("\n=== All plots saved to results/gwas/ ===\n")
cat("Files generated:\n")
cat("  - mrna_old_vs_new_SNP_count_histogram_log.png\n")
cat("  - mrna_old_vs_new_SNP_count_histogram_linear.png\n")
cat("  - mrna_old_vs_new_SNP_count_boxPlot.png\n")
cat("  - mrna_old_vs_new_SNP_density_histogram_log.png\n")
cat("  - mrna_old_vs_new_SNP_density_histogram_linear.png\n")
cat("  - mrna_old_vs_new_SNP_density_boxPlot.png\n")
cat("  - mrna_old_vs_new_density_vs_length_linear.png\n")
cat("  - mrna_old_vs_new_density_vs_length_loglog.png\n")
cat("  - mrna_old_vs_new_SNP_count_boxPlot_v7.png\n")
cat("  - mrna_old_vs_new_SNP_count_boxPlot_v49.png\n")
cat("  - mrna_old_vs_new_SNP_density_boxPlot_v7.png\n")
cat("  - mrna_old_vs_new_SNP_density_boxPlot_v49.png\n")
cat("  - mrna_old_vs_new_SNP_count_boxPlot_binary.png\n")
cat("  - mrna_old_vs_new_SNP_density_boxPlot_binary.png\n")
cat("  - mrna_old_vs_new_SNP_density_boxPlot_binary_v7.png\n")
cat("  - mrna_old_vs_new_SNP_density_boxPlot_binary_v49.png\n")
