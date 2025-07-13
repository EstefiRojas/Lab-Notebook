# Load libraries
library(dplyr)
library(ggplot2)
library(Hmisc)
library(ggsignif)
library(scales)  # For percentage formatting

# Load gwas association counts and minimum p-values file
lncrna_gwas_data <- read.csv("../data/datasets/gwas_pval_feature/latest_gencode-lncrna-ranking-5.3.25.csv", header = TRUE)

# Filter out NAs
lncrna_gwas_data_w <- lncrna_gwas_data

# Replace 0 with NA at number of gwas associations when the min p-value is NA for both transcripts
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  mutate(tl_gwas_assoc = case_when(
    is.na(tl_min_p_val) ~ NA,
    !is.na(tl_min_p_val) ~ tl_gwas_assoc
  ))

lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  mutate(tr_gwas_assoc = case_when(
    is.na(tr_min_p_val) ~ NA,
    !is.na(tr_min_p_val) ~ tr_gwas_assoc
  ))

# Get the maximum number of associations between transcripts
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  mutate(max_gwas_assoc = pmax(tl_gwas_assoc, tr_gwas_assoc, na.rm = TRUE))


# Compute transcript lengths
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  mutate(tl_length = End_Transcript_Left - Start_Transcript_Left,
         tr_length = End_Transcript_Right - Start_Transcript_Right)

# Compute SNP density per kilo base
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  mutate(tl_snp_density_kb = tl_gwas_assoc / tl_length * 1000,
         tr_snp_density_kb = tr_gwas_assoc / tr_length * 1000
         )

# Get the maximum SNP density per kilo base between transcripts
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  mutate(max_snp_density_kb = pmax(tl_snp_density_kb, tr_snp_density_kb, na.rm = TRUE))

# Check new fields created so far
summary(lncrna_gwas_data_w$tl_length)
summary(lncrna_gwas_data_w$tr_length)
summary(lncrna_gwas_data_w$max_snp_density_kb)
summary(lncrna_gwas_data_w$max_gwas_assoc)
summary(lncrna_gwas_data_w$tl_gwas_assoc)
summary(lncrna_gwas_data_w$tr_gwas_assoc)


# Group by SNPDensity
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  mutate(snp_density_category = case_when(
    is.na(max_snp_density_kb) ~ "NA",
    max_snp_density_kb > 0 & max_snp_density_kb <= 1 ~ ">0 - 1",
    max_snp_density_kb > 1 & max_snp_density_kb <= 5 ~ ">1 - 5",
    max_snp_density_kb > 5 ~ ">5"
  ))
# Count data points and update category labels
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  group_by(snp_density_category) %>%
  mutate(snp_density_label = paste0(snp_density_category, " (", n(), ")")) %>%
  ungroup()

unique(lncrna_gwas_data_w$snp_density_label)
# Ensure the legend labels are ordered as required
lncrna_gwas_data_w$snp_density_label <- factor(lncrna_gwas_data_w$snp_density_label, levels = c(
  "NA (6850)",">0 - 1 (6519)", ">1 - 5 (4339)", ">5 (394)")
)
#help("stat_ecdf")
#Plot cumulative distribution
ggplot(lncrna_gwas_data_w, aes(x = highest_prob, color = snp_density_label)) +
  stat_ecdf(geom = "step", linewidth = 1.2, ) +  # Increased line thickness
  scale_color_manual(
    values = c(
      "#1f77b4",  # Strong blue
      "#ffd700",  # Gold
      "#2ca02c",  # Strong green
      "#d62728",  # Bright red
      "#9467bd",  # Strong purple
      "#e377c2",  # Bright pink
      "#000000",  # Black
      "#bcbd22",  # Olive
      "#00ff00"  # Lime green
    )
  ) +
  coord_flip() +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +  # Format as percentages
  scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Format as percentages
  labs(
    title = "Cumulative Distribution by SNP Density",
    x = "Functional Probability",
    y = "Fraction of Genes",
    color = "GWAS SNP density per kb (lncRNAs)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    panel.grid.major = element_line(color = "gray90"),  # Lighter grid lines
    panel.grid.minor = element_line(color = "gray95")   # Lighter grid lines
  )

##########################
# GWAS Association counts
corr_obj_max_assoc <- rcorr(lncrna_gwas_data$highest_prob,
                       lncrna_gwas_data$max_gwas_assoc, 
                       type = "spearman")

# Group by number of associations
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  mutate(association_category = case_when(
    is.na(max_gwas_assoc) ~ "NA",
    max_gwas_assoc == 0 ~ "0",
    max_gwas_assoc >= 1 & max_gwas_assoc <= 10 ~ "1-10",
    max_gwas_assoc >= 11 & max_gwas_assoc <= 100 ~ "11-100",
    max_gwas_assoc >= 101 & max_gwas_assoc <= 500 ~ "101-500",
    max_gwas_assoc > 500 ~ ">500"
  ))
unique(lncrna_gwas_data_w$association_category)
# Count data points and update category labels
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  group_by(association_category) %>%
  mutate(label_with_count = paste0(association_category, " (", n(), ")")) %>%
  ungroup()

unique(lncrna_gwas_data_w$label_with_count)
# Ensure the legend labels are ordered as required
lncrna_gwas_data_w$label_with_count <- factor(lncrna_gwas_data_w$label_with_count, levels = c(
  "NA (6845)", "1-10 (7154)", "11-100 (3422)",
  "101-500 (633)", ">500 (48)")
  )

ggplot(lncrna_gwas_data_w, aes(x = highest_prob, color = label_with_count)) +
  stat_ecdf(geom = "step") +
  scale_color_manual(
    values = c("blue", "yellow4", "green", "red", "purple", "brown", "lightblue", "pink", "grey", "black", "lightgreen", "magenta", "yellow4")  # Custom colors
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Format as percentages
  labs(
    title = "Cumulative Distribution by GWAS Association Counts",
    x = "Functional Probability",
    y = "Fraction of Genes",
    color = "GWAS SNPs (Number of Genes)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    #legend.title = element_blank()
    #legend.position = "none"
  )



ggplot(lncrna_gwas_data_w, aes(x = highest_prob, color = label_with_count)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +  # Increased line thickness
  scale_color_manual(
    values = c(
      "#1f77b4",  # Strong blue
      "#ffd700",  # Gold
      "#2ca02c",  # Strong green
      "#d62728",  # Bright red
      "#9467bd",  # Strong purple
      
      "#e377c2",  # Bright pink
      "#000000",  # Black
      "#bcbd22",  # Olive
      "#00ff00"  # Lime green
      
      
    )
  ) +
  coord_flip() +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +  # Format as percentages
  scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Format as percentages
  labs(
    title = "Cumulative Distribution by GWAS Association Counts",
    x = "Functional Probability",
    y = "Fraction of Genes",
    color = "GWAS SNPs (lncRNAs)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    panel.grid.major = element_line(color = "gray90"),  # Lighter grid lines
    panel.grid.minor = element_line(color = "gray95")   # Lighter grid lines
  )




ggplot(lncrna_gwas_data_w, aes(y = highest_prob, x = max_gwas_assoc, color = label_with_count)) +
  geom_point() +
  scale_color_manual(
    values = c("blue", "darkorange2", "green", "red", "purple", "brown", "lightblue", "pink", "grey", "black", "lightgreen", "magenta", "yellow4")  # Custom colors
  ) +
  #scale_x_continuous(limits = c(0.25, 1)) +
  #scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "GWAS Association Counts",
    y = "Functional Probability",
    x = "GWAS SNPs",
    color = "GWAS SNPs (Number of Genes)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    #legend.title = element_blank()
    #legend.position = "top"
  )


# Minimum p-value
# Get the minimum p-value between both transcripts analysed
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  mutate(min_p_value = pmin(tl_min_p_val, tr_min_p_val, na.rm = TRUE))


corr_obj_minp <- rcorr(lncrna_gwas_data_w$highest_prob,
                  lncrna_gwas_data$min_p_value, 
                  type = "spearman")

# Group by 
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  mutate(pval_association_category = case_when(
    is.na(min_p_value) ~ "NA",
    min_p_value <= 1e-05 & min_p_value >= 1e-10 ~ "1e-05 - 1e-10",
    min_p_value < 1e-10 & min_p_value >= 1e-100 ~ "1e-11 - 1e-100",
    min_p_value < 1e-100 & min_p_value >= 1e-300 ~ "1e-101 - 1e-300",
    min_p_value < 1e-300 ~ "<1e-300"
  ))

unique(lncrna_gwas_data_w$pval_association_category)
# Count data points and update category labels
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  group_by(pval_association_category) %>%
  mutate(pval_label_with_count = paste0(pval_association_category, " (", n(), ")")) %>%
  ungroup()

unique(lncrna_gwas_data_w$pval_label_with_count)
# Ensure the legend labels are ordered as required
lncrna_gwas_data_w$pval_label_with_count <- factor(lncrna_gwas_data_w$pval_label_with_count, levels = c(
  "NA (6845)", "1e-05 - 1e-10 (7115)", "1e-11 - 1e-100 (3751)", "1e-101 - 1e-300 (150)",
  "<1e-300 (241)")
)

lncrna_gwas_data_w$pval_association_category <- factor(lncrna_gwas_data_w$pval_association_category, levels = c(
  "NA", "1e-05 - 1e-10", "1e-11 - 1e-100", "1e-101 - 1e-300",
  "<1e-300")
)

ggplot(lncrna_gwas_data_w, aes(x = highest_prob, color = pval_label_with_count)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  scale_color_manual(
    values = c(
      "#1f77b4",  # Strong blue
      "#ffd700",  # Gold
      "#2ca02c",  # Strong green
      "#d62728",  # Bright red
      "#9467bd",  # Strong purple
      "#e377c2",  # Bright pink
      "#000000",  # Black
      "#bcbd22",  # Olive
      "#00ff00"  # Lime green
    )
  ) +
  #scale_x_continuous(limits = c(0.25, 1)) +
  #scale_y_continuous(limits = c(0, 1)) +
  coord_flip() +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +  # Format as percentages
  scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Format as percentages
  labs(
    title = "Cumulative Distribution by GWAS min p-value",
    x = "Functional Probability",
    y = "Fraction of Genes",
    color = "GWAS min p-value (lncRNAs)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    #legend.title = element_blank()
    #legend.position = "top"
  )


ggplot(lncrna_gwas_data_w, aes(x = -log10(min_p_value), colour = label_with_count)) +
  stat_bin(alpha = 0.8, binwidth = 0.5, na.rm = TRUE) +
  scale_x_continuous(breaks = seq(0, 300, by = 20)) +
  scale_color_manual(
    values = c("blue", "darkorange2", "green", "red", "purple")
  ) +
  labs(
    title = "Distribution of GWAS P-values by Association Count",
    x = "-log10(P-value)",
    y = "Count",
    color = "GWAS SNPs (Number of Genes)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20)
  )

ggplot(lncrna_gwas_data_w, aes(x = highest_prob, y = min_p_value, color = label_with_count)) +
  geom_step() +
  scale_color_manual(
    values = c("blue", "darkorange2", "green", "red", "purple", "brown", "lightblue", "pink", "grey", "black", "lightgreen", "magenta", "yellow4")  # Custom colors
  ) +
  #scale_x_continuous(limits = c(0.25, 1)) +
  #scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Cumulative Distribution by GWAS minimum P-value of SNP associations",
    x = "Functional Probability",
    y = "Fraction of Genes",
    color = "GWAS SNPs (Number of Genes)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    #legend.title = element_blank()
    #legend.position = "top"
  )


# Plot density distribution
lncrna_gwas_data_w |> 
  mutate(Type = "lncRNA") |>
  filter(!is.na(min_p_value)) |>
  ggplot(aes(x = -log10(min_p_value), fill = Type, colour = Type)) +
  geom_histogram(alpha = 0.6, binwidth = 0.5) +
  geom_vline(xintercept = median(-log10(lncrna_gwas_data_w$min_p_value), na.rm = TRUE)) +
  labs(title = "GWAS catalog - min P-value lncRNA") + 
  #theme_minimal() +
  scale_fill_manual(values = c("#e37b88FF")) +
  scale_x_continuous(breaks = seq(0, 300, by = 20)) +
  #scale_x_reverse() +
  #scale_color_manual(values = c("#53a4f5FF")) +
  theme(axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32),  # Increase y-axis text size
        axis.title.x = element_text(size = 44),  # Increase x-axis title size
        axis.title.y = element_text(size = 44),  # Increase y-axis title size
        legend.position = "none",
        plot.title = element_text(size = 64, face = "italic", hjust = 0.5)  # center the title
  )
summary(-log10(lncrna_gwas_data_w$min_p_value))
# Plot density distribution
lncrna_gwas_data_w |> 
  mutate(Type = "lncRNA") |>
  filter(min_p_value>0) |>
  ggplot(aes(x = -log10(min_p_value), y = highest_prob, fill = Type, colour = Type)) +
  geom_point() +
  geom_vline(xintercept = median(-log10(lncrna_gwas_data_w$min_p_value), na.rm = TRUE)) +
  labs(title = "GWAS catalog - min P-value",
       x="p-value (log scale)",
       y="Fuctional Probability") + 
  theme_minimal() +
  scale_fill_manual(values = c("#e37b88FF")) +
  scale_x_continuous(limits = c(0, 310), breaks = seq(0, 300, by = 20)) +
  #coord_cartesian(xlim = c(0, 290)) +
  #scale_x_reverse() +
  #scale_color_manual(values = c("#53a4f5FF")) +
  theme(axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32),  # Increase y-axis text size
        axis.title.x = element_text(size = 44),  # Increase x-axis title size
        axis.title.y = element_text(size = 44),  # Increase y-axis title size
        legend.position = "none",
        plot.title = element_text(size = 64, face = "italic", hjust = 0.5)  # center the title
  )

# --- Data Preparation ---
# 1. Filter out non-positive p-values (important for log10)
# 2. Calculate the -log10(p-value)
# 3. Ensure 'pval_label_with_count' is treated as a factor if needed
plot_data_prep <- lncrna_gwas_data_w |>
  filter(min_p_value > 0) |>
  mutate(
    neg_log10_p = -log10(min_p_value),
    # Ensure the grouping variable is a factor
    pval_label_with_count = factor(pval_label_with_count)
  )

# Calculate the median from the *filtered* data's transformed p-values
median_neg_log10_p <- median(plot_data_prep$neg_log10_p, na.rm = TRUE)

# --- Calculate Spearman Correlations and Create Legend Labels ---
# Group by the category, calculate correlation
correlation_labels <- plot_data_prep |>
  group_by(pval_label_with_count) |>
  # Calculate Spearman correlation; requires at least 2 complete pairs per group
  summarise(
    correlation = cor(neg_log10_p, highest_prob, method = "spearman", use = "pairwise.complete.obs"),
    .groups = 'drop' # Drop grouping after summarise
  ) |>
  # Create the combined label string for the legend
  mutate(
    legend_label = case_when(
      is.na(correlation) ~ sprintf("%s (NA)", pval_label_with_count), # Handle NA case
      TRUE ~ sprintf("%s (ρ = %.2f)", pval_label_with_count, correlation) # Format if valid
    )
  ) |>
  # Select only the original category name and the new combined label
  select(pval_label_with_count, legend_label)

# --- Join Legend Labels back to Plot Data ---
# Join the created labels back to the main data used for plotting
plot_data <- plot_data_prep |>
  left_join(correlation_labels, by = "pval_label_with_count") |>
  # Create a factor from the legend labels for ggplot mapping.
  # Order levels based on the original factor to maintain consistency (optional but good practice)
  mutate(legend_label_factor = factor(legend_label, levels = unique(legend_label[order(pval_label_with_count)])))

# --- Plotting ---
# Map 'colour' aesthetic to the new 'legend_label_factor' variable
ggplot(plot_data, aes(x = neg_log10_p, y = highest_prob, colour = legend_label_factor)) +
  # geom_point no longer needs direct color setting
  geom_point(size = 2) +
  # Add vertical line using the pre-calculated median of the plotted data
  geom_vline(xintercept = median_neg_log10_p, linetype = "dashed", colour = "grey50") +
  # Add a color scale. Brewer palettes are good for categories.
  # The legend will now use the values from 'legend_label_factor'
  scale_colour_brewer(palette = "Set1") +
  labs(
    title = "GWAS catalog - min P-value",
    x = "-log10(p-value)", # More precise axis label
    y = "Functional Probability",
    colour = "Group (Spearman ρ)" # Update legend title
  ) +
  theme_minimal() +
  # Set x-axis limits and breaks (adjust if needed based on your data range)
  scale_x_continuous(limits = c(0, 310), breaks = seq(0, 300, by = 40)) +
  # Removed specific y-axis limits, ggplot will set defaults
  # scale_y_continuous(limits = c(NA, 1.05)) + # Re-add if needed for spacing
  # Theme adjustments - Legend is shown by default. Customize if needed.
  theme(
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),  # Increase y-axis text size
    axis.title.x = element_text(size = 44),  # Increase x-axis title size
    axis.title.y = element_text(size = 44),  # Increase y-axis title size
    plot.title = element_text(size = 64, face = "italic", hjust = 0.5),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 26),
    legend.position = "right" # Explicitly position legend
  )




# Custom function to perform K-S test and format the D-statistic and p-value stars
ks_test_custom <- function(x, y) {
  test <- ks.test(x, y)
  
  # Convert p-value to significance stars
  # 'cut' is a clean way to assign stars based on p-value ranges
  stars <- cut(test$p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
               labels = c("***", "**", "*", "ns"))
  
  # Format the label as: D-value (stars)
  label <- paste0(format(test$statistic, digits = 2), " (", format(test$p.value, digits = 2), ")")
  
  # Return a list with the custom label for ggsignif
  return(list(p.value = label))
}

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("NA","1e-05 - 1e-10"),
                       c("NA","1e-11 - 1e-100"),
                       c("NA","1e-101 - 1e-300"),
                       c("NA","<1e-300"))

# Create violin plots
violin_plot <- ggplot(lncrna_gwas_data_w, aes(x = pval_association_category, y = highest_prob, fill = pval_label_with_count)) +
  geom_violin(scale="area", na.rm = TRUE, trim = FALSE, bounds = c(0,1)) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Functional probability distribution by GWAS minimum p-value",
    x = "GWAS SNP minimum p-value",
    y = "Highest Probability",
    fill = "Minimum p-value" # Update legend title to match 'color'
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 28),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


violin_plot_with_pvals <- violin_plot +
  geom_signif(
    comparisons = my_comparisons,
    test = "ks_test_custom",
    step_increase = 0.2,
    textsize = 6.5,
    tip_length = 0.01,
    y_position = 1.2
  ) +
  # This prevents nonsensical axis ticks (e.g., > 1) while keeping room for annotations.
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 2.1), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations

#violin_plot +
#  annotate("segment", x = 1, xend = 2, y = 1.2, yend = 1.2) +  # Line for NA to 1e-05 - 1e-10
#  annotate("text", x = 1.5, y = 1.25, label = "0.032 (0.001)", size = 5) +   # Text for first p-value
#  annotate("segment", x = 1, xend = 3, y = 1.3, yend = 1.3) +  # Line for NA to 1e-11 - 1e-100
#  annotate("text", x = 2, y = 1.35, label = "0.129 (2.79e-35)", size = 5) +   # Text for second p-value
#  annotate("segment", x = 1, xend = 4, y = 1.4, yend = 1.4) +  # Line for NA to 1e-101 - 1e-300
#  annotate("text", x = 2.5, y = 1.45, label = "0.099 (0.110)", size = 5) +   # Text for second p-value
#  annotate("segment", x = 1, xend = 5, y = 1.5, yend = 1.5) +  # Line for NA to <1e-300
#  annotate("text", x = 3, y = 1.55, label = "0.164 (7.03e-06)", size = 5) +   # Text for second p-value
#  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
#  coord_cartesian(ylim = c(0, 1.7), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations

violin_plot_with_pvals


##################################
# Associations KS stat computation
##################################
# Compute KS between first and last bin
assoc_max_assoc_bin1 <- lncrna_gwas_data_w |>
  filter(association_category=="NA")

assoc_max_assoc_bin2 <- lncrna_gwas_data_w |>
  filter(association_category=="1-10")

assoc_max_assoc_bin3 <- lncrna_gwas_data_w |>
  filter(association_category=="11-100")

assoc_max_assoc_bin4 <- lncrna_gwas_data_w |>
  filter(association_category=="101-500")

assoc_max_assoc_bin5 <- lncrna_gwas_data_w |>
  filter(association_category==">500")

colnames(assoc_max_assoc_bin1)

assoc_ks_results_bin1_5 <- run_ks_tests(assoc_max_assoc_bin1, 
                                        assoc_max_assoc_bin5, 
                                        c("highest_prob"))
assoc_ks_results_bin1_4 <- run_ks_tests(assoc_max_assoc_bin1, 
                                        assoc_max_assoc_bin4, 
                                        c("highest_prob"))
assoc_ks_results_bin1_3 <- run_ks_tests(assoc_max_assoc_bin1, 
                                        assoc_max_assoc_bin3, 
                                        c("highest_prob"))
assoc_ks_results_bin1_2 <- run_ks_tests(assoc_max_assoc_bin1, 
                                        assoc_max_assoc_bin2, 
                                        c("highest_prob"))

assoc_comparisons <- list(
  c("NA", "1-10"), 
  c("NA", "11-100"),
  c("NA", "101-500"),
  c("NA", ">500")
)
assoc_ks_results_bin1_5
assoc_p_values <- c(2.780165e-06, 1.557839e-12, 3.421503e-08, 0.05397672)

# Create violin plots
assoc_violin_plot <- ggplot(lncrna_gwas_data_w, aes(x = label_with_count, y = highest_prob, fill = label_with_count)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Functional probability distribution by GWAS associations",
    x = "GWAS SNP count",
    y = "Highest Probability"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

assoc_violin_plot +
  annotate("segment", x = 1, xend = 2, y = 1.2, yend = 1.2) +  # Line for NA to 1e-05 - 1e-10
  annotate("text", x = 1.5, y = 1.25, label = "0.044 (****)") +   # Text for first p-value
  annotate("segment", x = 1, xend = 3, y = 1.3, yend = 1.3) +  # Line for NA to 1e-11 - 1e-100
  annotate("text", x = 2, y = 1.35, label = "0.078 (****)") +   # Text for second p-value
  annotate("segment", x = 1, xend = 4, y = 1.4, yend = 1.4) +  # Line for NA to 1e-101 - 1e-300
  annotate("text", x = 2.5, y = 1.45, label = "0.124 (****)") +   # Text for second p-value
  annotate("segment", x = 1, xend = 5, y = 1.5, yend = 1.5) +  # Line for NA to <1e-300
  annotate("text", x = 3, y = 1.55, label = "0.195 (-)") +   # Text for second p-value
  # ... similar lines for other comparisons
  ylim(-0.5, 1.6)

##################################
# SNP density KS stat computation
##################################
# Compute KS between first and last bin
snp_den_max_assoc_bin1 <- lncrna_gwas_data_w |>
  filter(snp_density_category=="NA")

snp_den_max_assoc_bin2 <- lncrna_gwas_data_w |>
  filter(snp_density_category==">0 - 1")

snp_den_max_assoc_bin3 <- lncrna_gwas_data_w |>
  filter(snp_density_category==">1 - 5")

snp_den_max_assoc_bin4 <- lncrna_gwas_data_w |>
  filter(snp_density_category==">5")

colnames(snp_den_max_assoc_bin1)


snp_den_ks_results_bin1_4 <- run_ks_tests(snp_den_max_assoc_bin1, 
                                        snp_den_max_assoc_bin4, 
                                        c("highest_prob"))
snp_den_ks_results_bin1_3 <- run_ks_tests(snp_den_max_assoc_bin1, 
                                          snp_den_max_assoc_bin3, 
                                        c("highest_prob"))
snp_den_ks_results_bin1_2 <- run_ks_tests(snp_den_max_assoc_bin1, 
                                          snp_den_max_assoc_bin2, 
                                        c("highest_prob"))

snp_den_comparisons <- list(
  c("NA", ">0 - 1"), 
  c("NA", ">1 - 5"),
  c("NA", ">5")
)
snp_den_ks_results_bin1_4
snp_den_p_values <- c(1.407950e-07, 3.294273e-29, 6.971807e-12)

# Create violin plots
snp_den_violin_plot <- ggplot(lncrna_gwas_data_w, aes(x = snp_density_label, y = highest_prob, fill = snp_density_label)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Functional probability distribution by GWAS SNP density per kb",
    x = "GWAS SNP density per kb",
    y = "Highest Probability"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

snp_den_violin_plot +
  annotate("segment", x = 1, xend = 2, y = 1.2, yend = 1.2) +  # Line for NA to 1e-05 - 1e-10
  annotate("text", x = 1.5, y = 1.25, label = "0.049 (****)") +   # Text for first p-value
  annotate("segment", x = 1, xend = 3, y = 1.3, yend = 1.3) +  # Line for NA to 1e-11 - 1e-100
  annotate("text", x = 2, y = 1.35, label = "0.112 (****)") +   # Text for second p-value
  annotate("segment", x = 1, xend = 4, y = 1.4, yend = 1.4) +  # Line for NA to 1e-101 - 1e-300
  annotate("text", x = 2.5, y = 1.45, label = "-0.188 (****)") +   # Text for second p-value

  ylim(-0.25, 1.6)




##############################################
# Analyze rows classified as functional only #
##############################################
lncrna_gwas_data_top_prob <- lncrna_gwas_data_w %>%
  filter(highest_prob >= 0.5)

# Group by SNPDensity
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  mutate(highes_prob_category = case_when(
    is.na(highest_prob) ~ "NA",
    highest_prob >= 0 & highest_prob <= 0.1 ~ "0% - 10%",
    highest_prob > 0.1 & highest_prob <= 0.2 ~ "10% - 20%",
    highest_prob > 0.2 & highest_prob <= 0.4 ~ "20% - 40%",
    highest_prob > 0.4 & highest_prob <= 0.7 ~ "40% - 70%",
    highest_prob > 0.7 ~ "70% - 100%"
  ))
# Count data points and update category labels
lncrna_gwas_data_w <- lncrna_gwas_data_w %>%
  group_by(highes_prob_category) %>%
  mutate(highes_prob_label = paste0(highes_prob_category, " (", n(), ")")) %>%
  ungroup()

unique(lncrna_gwas_data_w$highes_prob_label)
# Ensure the legend labels are ordered as required
lncrna_gwas_data_w$highes_prob_label <- factor(lncrna_gwas_data_w$highes_prob_label, levels = c(
  "0% - 10% (3827)",
  "10% - 20% (3574)",
  "20% - 40% (3553)",
  "40% - 70% (3938)",
  "70% - 100% (3206)")
)

#Plot cumulative distribution
ggplot(lncrna_gwas_data_w, aes(x = tl_gwas_assoc, color = highes_prob_label)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +  # Increased line thickness
  scale_color_manual(
    values = c(
      "#1f77b4",  # Strong blue
      "#2ca02c",  # Strong green
      "#ffd700",  # Gold
      "#d62728",  # Bright red
      "#9467bd",  # Strong purple
      "#e377c2",  # Bright pink
      "#bcbd22",  # Olive
      "#00ff00",  # Lime green
      "#ff7f0e",  # Bright orange
      "#17becf",  # Cyan/teal
      "#8c564b",  # Brown
      "#e6194B",  # Crimson
      "#3cb44b",  # Emerald green
      "#ffe119",  # Bright yellow
      "#4363d8",  # Royal blue
      "#f58231",  # Burnt orange
      "#911eb4",  # Deep purple
      "#42d4f4",  # Sky blue
      "#f032e6"   # Magenta
    )
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Format as percentages
  scale_x_log10() +
  coord_flip() +
  labs(
    title = "Cumulative Distribution of GWAS SNPs by Functional Probability",
    x = "GWAS SNPs",
    y = "Fraction of Genes",
    color = "Functional Probability (lncRNAs)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    panel.grid.major = element_line(color = "gray90"),  # Lighter grid lines
    panel.grid.minor = element_line(color = "gray95")   # Lighter grid lines
  )


ggplot(lncrna_gwas_data_w, aes(x = tl_snp_density_kb, color = highes_prob_label)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +  # Increased line thickness
  scale_color_manual(
    values = c(
      "#1f77b4",  # Strong blue
      "#2ca02c",  # Strong green
      "#ffd700",  # Gold
      "#d62728",  # Bright red
      "#9467bd",  # Strong purple
      "#e377c2",  # Bright pink
      "#bcbd22",  # Olive
      "#00ff00",  # Lime green
      "#ff7f0e",  # Bright orange
      "#17becf",  # Cyan/teal
      "#8c564b",  # Brown
      "#e6194B",  # Crimson
      "#3cb44b",  # Emerald green
      "#ffe119",  # Bright yellow
      "#4363d8",  # Royal blue
      "#f58231",  # Burnt orange
      "#911eb4",  # Deep purple
      "#42d4f4",  # Sky blue
      "#f032e6"   # Magenta
    )
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Format as percentages
  scale_x_log10() +
  coord_flip() +
  labs(
    title = "Cumulative Distribution of GWAS SNP density by Functional Probability assigned",
    x = "GWAS SNP density",
    y = "Fraction of Genes",
    color = "Functional Probability (lncRNAs)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    panel.grid.major = element_line(color = "gray90"),  # Lighter grid lines
    panel.grid.minor = element_line(color = "gray95")   # Lighter grid lines
  )

  # Create violin plots
  ggplot(lncrna_gwas_data_w, aes(x = highes_prob_label, y = tl_gwas_assoc, fill = highes_prob_label)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
    scale_color_manual(
      values = c(
        "#1f77b4",  # Strong blue
        "#2ca02c",  # Strong green
        "#ffd700",  # Gold
        "#d62728",  # Bright red
        "#9467bd",  # Strong purple
        "#e377c2",  # Bright pink
        "#bcbd22",  # Olive
        "#00ff00",  # Lime green
        "#ff7f0e",  # Bright orange
        "#17becf",  # Cyan/teal
        "#8c564b",  # Brown
        "#e6194B",  # Crimson
        "#3cb44b",  # Emerald green
        "#ffe119",  # Bright yellow
        "#4363d8",  # Royal blue
        "#f58231",  # Burnt orange
        "#911eb4",  # Deep purple
        "#42d4f4",  # Sky blue
        "#f032e6"   # Magenta
      )
    ) +
    scale_y_continuous(trans = "pseudo_log") +
    labs(
      title = "GWAS SNP counts distribution by Functional probability",
      x = "Highest Probability",
      y = "GWAS SNP count"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 20),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  