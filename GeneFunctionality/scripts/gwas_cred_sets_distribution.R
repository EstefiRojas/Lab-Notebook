# Load libraries
library(dplyr)
library(tidyr)
library(stringr) # For str_detect
library(ggplot2)
library(ggsignif)

# Load csv file with PIP = 1 from query
#gwas_var_pip_data <- read.csv("../data/gwas/gwas_cred_sets/bq-results-20250813-220807-1755122994396.csv", 
#                              header = TRUE)
# Load csv file with PIP = 0.95 from query
gwas_var_pip_data <- read.csv("../data/gwas/gwas_cred_sets/bq-results-20250814-021949-1755138038271.csv", 
                              header = TRUE)

# Extract chromosome and position from variantId
gwas_var_pip_data <- gwas_var_pip_data %>%
  separate(variantId, 
           into = c("chromosome", "position", "nucl1", "nucl2"), 
           sep = "_",
           remove = FALSE)
#help("write.csv")
write.csv(gwas_var_pip_data, file="../data/gwas/gwas_cred_sets/gwas_cred_sets_095.csv", sep = ",", na = "NA", row.names = FALSE, quote = FALSE)

# RUN get_pip_sum.sh script, which uses bedtools intercept to match with lncrna data,
# then continue with the rest if this script.
# That is, run the following on the terminal
#./get_pip_sum.sh ../data/gwas/lncrna-gwas-ranked-5-aug-2025.csv ../data/gwas/gwas_cred_sets/gwas_cred_sets.csv ../results/ranked_lncrna_gwas_cred_sets_pips.csv

# Load lncRNA 1 data:
lncrna_prob_data <- read.csv("../results/ranked_lncrna_gwas_cred_sets_pips.csv", header = TRUE) %>%
  separate(Probability_Functional, 
           into = c("tl_prob", "tr_prob"), 
           sep = "\\|",
           remove = FALSE) %>%
  mutate(tl_prob = as.numeric(tl_prob),
         tr_prob = as.numeric(tr_prob))

#####################
# Sum PIP analysis #
filterd_pip_data <- lncrna_prob_data %>% 
  mutate(
    tl_pip_per_kb = 1000 * (tl_sum_pip / (End_Transcript_Left - Start_Transcript_Left)),
    tr_pip_per_kb = 1000 * (tr_sum_pip / (End_Transcript_Right - Start_Transcript_Right)),
    selected_pip = ifelse(!is.na(tr_prob), 
                               ifelse(tl_prob == tr_prob, 
                                      ifelse(is.na(tr_pip_per_kb), 
                                             tl_pip_per_kb,
                                             ifelse(is.na(tl_pip_per_kb),
                                                    tr_pip_per_kb,
                                                    ifelse(tl_pip_per_kb >= tr_pip_per_kb,
                                                           tl_pip_per_kb,
                                                           tr_pip_per_kb
                                                    )
                                             )
                                      ), 
                                      ifelse(tl_prob > tr_prob, 
                                             tl_pip_per_kb, 
                                             tr_pip_per_kb
                                      )
                               ), 
                               tl_pip_per_kb
    ),
    # Then classify based on the selected beta
    pip_category = case_when(
      selected_pip <= 0.014277408 ~ ">0",
      selected_pip >= 0.014277408 & selected_pip < 0.020334379 ~ ">0.01",
      selected_pip > 0.020334379 & selected_pip < 0.026895388 ~ ">0.02",
      
      selected_pip >= 0.026895388 & selected_pip < 0.034763746 ~ ">0.025",
      selected_pip > 0.034763746 & selected_pip < 0.045502116 ~ ">0.035",
      
      selected_pip >= 0.045502116 & selected_pip < 0.060745600 ~ ">0.045",
      selected_pip > 0.060745600 & selected_pip < 0.084308511 ~ ">0.06",
      
      selected_pip >= 0.084308511 & selected_pip < 0.129309238 ~ ">0.085",
      selected_pip > 0.129309238 & selected_pip < 0.252665625 ~ ">0.13",
      selected_pip >= 0.252665625 ~ "≥0.25"
    ),
    tl_cred_sets_per_kb = tl_cred_set_count / (End_Transcript_Left - Start_Transcript_Left) * 1000,
    tr_cred_sets_per_kb = tr_cred_set_count / (End_Transcript_Right - Start_Transcript_Right) * 1000,
    selected_cred_set_per_kb = ifelse(!is.na(tr_prob), 
                          ifelse(tl_prob == tr_prob, 
                                 ifelse(is.na(tr_cred_sets_per_kb), 
                                        tl_cred_sets_per_kb,
                                        ifelse(is.na(tl_cred_sets_per_kb),
                                               tr_cred_sets_per_kb,
                                               ifelse(tl_cred_sets_per_kb >= tr_cred_sets_per_kb,
                                                      tl_cred_sets_per_kb,
                                                      tr_cred_sets_per_kb
                                               )
                                        )
                                 ), 
                                 ifelse(tl_prob > tr_prob, 
                                        tl_cred_sets_per_kb, 
                                        tr_cred_sets_per_kb
                                 )
                          ), 
                          tl_cred_sets_per_kb
    ),
    # Then classify based on the selected count of cred_sets
    cred_set_category = case_when(
      is.na(selected_cred_set_per_kb) ~ "NA",
      selected_cred_set_per_kb <= 0.014354115 ~ ">0",
      selected_cred_set_per_kb >= 0.014354115 & selected_cred_set_per_kb < 0.020465774 ~ ">0.015",
      selected_cred_set_per_kb > 0.020465774 & selected_cred_set_per_kb < 0.027204503 ~ ">0.02",
      
      selected_cred_set_per_kb >= 0.027204503 & selected_cred_set_per_kb < 0.035138456 ~ ">0.03",
      selected_cred_set_per_kb > 0.035138456 & selected_cred_set_per_kb < 0.045810619 ~ ">0.035",
      
      selected_cred_set_per_kb >= 0.045810619 & selected_cred_set_per_kb < 0.061348999 ~ ">0.045",
      selected_cred_set_per_kb > 0.061348999 & selected_cred_set_per_kb < 0.085045975 ~ ">0.06",
      
      selected_cred_set_per_kb >= 0.085045975 & selected_cred_set_per_kb < 0.129782510 ~ ">0.085",
      selected_cred_set_per_kb > 0.129782510 & selected_cred_set_per_kb < 0.254465971 ~ ">0.13",
      
      selected_cred_set_per_kb >= 0.254465971 ~ "≥0.25"
    )
  ) %>%
  select(highest_prob,tl_prob,tr_prob,
         tl_sum_pip,tr_sum_pip,selected_pip,pip_category,
         tl_cred_sets_per_kb,tr_cred_sets_per_kb,selected_cred_set_per_kb,cred_set_category) #%>%
  #filter(!is.na(selected_pip))
summary(filterd_pip_data$selected_pip)
summary(filterd_pip_data$selected_cred_set_per_kb)

table(filterd_pip_data$pip_category)
table(filterd_pip_data$cred_set_category)

# Get the decile probabilities
decile_probs <- seq(0, 1, by = 0.1)

# 3. Calculate the decile values
decile_values <- quantile(filterd_pip_data$selected_pip, probs = decile_probs, na.rm = TRUE)
decile_values <- quantile(filterd_pip_data$selected_cred_set_per_kb, probs = decile_probs, na.rm = TRUE)

# Print the results
print(decile_values)

filterd_pip_data$pip_category <- factor(filterd_pip_data$pip_category,
                                                    levels = c(">0",">0.01",">0.02",">0.025",">0.035",
                                                               ">0.045",">0.06",">0.085",">0.13","≥0.25"))
filterd_pip_data$cred_set_category <- factor(filterd_pip_data$cred_set_category,
                                        levels = c("NA",">0",">0.015",">0.02",">0.03",">0.035",">0.045",
                                                   ">0.06",">0.085",">0.13","≥0.25"))

# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filterd_pip_data %>%
  count(pip_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$pip_category, "\n(n=", legend_data$n, ")"),
  legend_data$pip_category
)

# Histogram in log scale
ggplot(filterd_pip_data, aes(x = selected_pip)) +
  geom_histogram(bins = 100) +
  scale_x_log10(breaks=c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100)) +
  labs(
    title = "GWAS PIP",
    x = "Sum PIP (log scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Histogram normal scale
ggplot(filterd_pip_data, aes(x = selected_pip)) +
  geom_histogram(bins = 100) +
  #coord_cartesian(c(0,6)) +
  labs(
    title = "GWAS PIP",
    x = "Sum PIP (linear scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c(">0",">0.01"),
                       c(">0",">0.02"),
                       c(">0",">0.025"),
                       c(">0",">0.035"),
                       c(">0",">0.045"),
                       c(">0",">0.06"),
                       c(">0",">0.085"),
                       c(">0",">0.13"),
                       c(">0","≥0.25"))

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


# Compute KS stats

# --- Step 1: Pre-calculate statistics for labels ---
reference_group <- ">0"
comparison_groups <- c(">0.01",">0.02",">0.025",">0.035",
                       ">0.045",">0.06",">0.085",">0.13","≥0.25")

# Extract the data for the reference group
reference_data <- filterd_pip_data %>%
  filter(pip_category == reference_group) %>%
  pull(highest_prob)

# Calculate KS statistic for each comparison group
stats_list <- lapply(comparison_groups, function(group) {
  # Extract data for the current comparison group
  comparison_data <- filterd_pip_data %>%
    filter(pip_category == group) %>%
    pull(highest_prob)
  
  # Perform the KS test
  ks_result <- ks.test(reference_data, comparison_data)
  
  # Return a data frame with the necessary info for plotting
  data.frame(
    pip_category = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

# Combine the list of data frames into a single data frame
stats_labels <- do.call(rbind, stats_list)

# Define the y-position for the labels (adjust as needed)
stats_labels$y_position <- 1.1


plot_modified <- ggplot(filterd_pip_data, aes(x = pip_category, y = highest_prob, fill = pip_category)) +
  #geom_violin(scale = "area") +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_manual(values = paletteer_d("ggsci::planetexpress_futurama")) +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS Sum PIP (kb)",
    #x = "Sum PIP per kb",
    y = "lncRNA Probability",
    caption = "*p-val<5e-8"
  ) +
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    text = element_text(size = 36),
    plot.title = element_text(size = 46, hjust = 0.5),
    axis.text.x = element_text(hjust = 0.5, size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(size = 44),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
    plot.caption = element_text(size = 20, hjust = -0.09, face = "bold.italic", color = "grey40")
  )

# Add the statistical comparison layer
plot_with_text_stats <- plot_modified +
  
  # Add the pre-calculated stats as text
  geom_text(
    data = stats_labels,
    aes(x = pip_category, y = y_position, label = label),
    inherit.aes = FALSE, 
    size = 6.5,
    color = "black",
    fontface = "bold"
  ) +
  
  # Adjust y-axis limits to ensure text is visible
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.15), clip = "off")
#plot_with_stats <- plot_modified +
#  geom_signif(
#    comparisons = my_comparisons,
#    test = "ks_test_custom",
#    step_increase = 0.2,
#    textsize = 9.5,
#    tip_length = 0.01,
#    y_position = 1.2
#  ) +
#  # KEY CHANGE: Set explicit breaks for the y-axis and use coord_cartesian to set the visual range.
#  # This prevents nonsensical axis ticks (e.g., > 1) while keeping room for annotations.
#  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
#  coord_cartesian(ylim = c(0, 1.7), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations


# Display the final plot
print(plot_with_text_stats)
#print(plot_with_stats)



########################################
# - Number of credible sets analysis - #
# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filterd_pip_data %>%
  count(cred_set_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$cred_set_category, "\n(n=", legend_data$n, ")"),
  legend_data$cred_set_category
)

# Histogram in log scale
ggplot(filterd_pip_data, aes(x = selected_cred_set_per_kb)) +
  geom_histogram(bins = 100) +
  scale_x_log10(breaks=c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100)) +
  labs(
    title = "GWAS Credible Sets",
    x = "Number of credible sets (log scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Histogram normal scale
ggplot(filterd_pip_data, aes(x = selected_cred_set_per_kb)) +
  geom_histogram(bins = 100) +
  #coord_cartesian(c(0,6)) +
  labs(
    title = "GWAS Credible Sets",
    x = "Number of credible sets (linear scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("NA",">0"),
                       c("NA",">0.015"),
                       c("NA",">0.02"),
                       c("NA",">0.035"),
                       c("NA",">0.045"),
                       c("NA",">0.06"),
                       c("NA",">0.085"),
                       c("NA",">0.13"),
                       c("NA","≥0.25"))

# Compute KS stats

# --- Step 1: Pre-calculate statistics for labels ---
reference_group <- "NA"
comparison_groups <- c(">0",">0.015",">0.02",">0.03",">0.035",">0.045",
                       ">0.06",">0.085",">0.13","≥0.25")

# Extract the data for the reference group
reference_data <- filterd_pip_data %>%
  filter(cred_set_category == reference_group) %>%
  pull(highest_prob)

# Calculate KS statistic for each comparison group
stats_list <- lapply(comparison_groups, function(group) {
  # Extract data for the current comparison group
  comparison_data <- filterd_pip_data %>%
    filter(cred_set_category == group) %>%
    pull(highest_prob)
  
  # Perform the KS test
  ks_result <- ks.test(reference_data, comparison_data)
  
  # Return a data frame with the necessary info for plotting
  data.frame(
    cred_set_category = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

# Combine the list of data frames into a single data frame
stats_labels <- do.call(rbind, stats_list)

# Define the y-position for the labels (adjust as needed)
stats_labels$y_position <- 1.1


plot_modified <- ggplot(filterd_pip_data, aes(x = cred_set_category, y = highest_prob, fill = cred_set_category)) +
  #geom_violin(scale = "area") +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_manual(values = paletteer_d("ggsci::planetexpress_futurama")) +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS: Causal SNPs Open Targets",
    #x = "Number of Credible sets per kb",
    y = "lncRNA Probability",
    caption = "*p-val<5e-8"
  ) +
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    text = element_text(size = 36),
    plot.title = element_text(size = 46, hjust = 0.5),
    axis.text.x = element_text(hjust = 0.5, size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(size = 44),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
    plot.caption = element_text(size = 20, hjust = -0.09, face = "bold.italic", color = "grey40")
  )

# Add the statistical comparison layer
plot_with_text_stats <- plot_modified +
  
  # Add the pre-calculated stats as text
  geom_text(
    data = stats_labels,
    aes(x = cred_set_category, y = y_position, label = label),
    inherit.aes = FALSE, 
    size = 6.5,
    color = "black",
    fontface = "bold"
  ) +
  
  # Adjust y-axis limits to ensure text is visible
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.15), clip = "off")
#plot_with_stats <- plot_modified +
#  geom_signif(
#    comparisons = my_comparisons,
#    test = "ks_test_custom",
#    step_increase = 0.2,
#    textsize = 9.5,
#    tip_length = 0.01,
#    y_position = 1.2
#  ) +
  # KEY CHANGE: Set explicit breaks for the y-axis and use coord_cartesian to set the visual range.
  # This prevents nonsensical axis ticks (e.g., > 1) while keeping room for annotations.
#  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
#  coord_cartesian(ylim = c(0, 1.9), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations


# Display the final plot
print(plot_with_text_stats)
#print(plot_with_stats)






############################
#--------------------------#
# Threshold 0.95 cred sets #
#####################
# Load lncRNA 0.95 data:
lncrna_prob_data <- read.csv("../results/ranked_lncrna_gwas_cred_sets_pips_095.csv", header = TRUE) %>%
  separate(Probability_Functional, 
           into = c("tl_prob", "tr_prob"), 
           sep = "\\|",
           remove = FALSE) %>%
  mutate(tl_prob = as.numeric(tl_prob),
         tr_prob = as.numeric(tr_prob))

# Sum PIP analysis #
filterd_pip_data <- lncrna_prob_data %>% 
  mutate(
    tl_pip_per_kb = tl_sum_pip / (End_Transcript_Left - Start_Transcript_Left) * 1000,
    tr_pip_per_kb = tr_sum_pip / (End_Transcript_Right - Start_Transcript_Right) * 1000,
    selected_pip = ifelse(!is.na(tr_prob), 
                          ifelse(tl_prob == tr_prob, 
                                 ifelse(is.na(tr_pip_per_kb), 
                                        tl_pip_per_kb,
                                        ifelse(is.na(tl_pip_per_kb),
                                               tr_pip_per_kb,
                                               ifelse(tl_pip_per_kb >= tr_pip_per_kb,
                                                      tl_pip_per_kb,
                                                      tr_pip_per_kb
                                               )
                                        )
                                 ), 
                                 ifelse(tl_prob > tr_prob, 
                                        tl_pip_per_kb, 
                                        tr_pip_per_kb
                                 )
                          ), 
                          tl_pip_per_kb
    ),
    # Then classify based on the selected beta
    pip_category = case_when(
      selected_pip <= 0.065 ~ "≤0.07",
      #selected_pip >= 0.032 & selected_pip < 0.065 ~ "0.03-0.07",
      #selected_pip > 0.031 & selected_pip <= 0.065 ~ ">0.03",
      selected_pip > 0.065 & selected_pip <= 0.148 ~ ">0.07",
      selected_pip >= 0.148 ~ "≥0.15"
    ),
    tl_cred_sets_per_kb = tl_cred_set_count / (End_Transcript_Left - Start_Transcript_Left) * 1000,
    tr_cred_sets_per_kb = tr_cred_set_count / (End_Transcript_Right - Start_Transcript_Right) * 1000,
    selected_cred_set_per_kb = ifelse(!is.na(tr_prob), 
                                      ifelse(tl_prob == tr_prob, 
                                             ifelse(is.na(tr_cred_sets_per_kb), 
                                                    tl_cred_sets_per_kb,
                                                    ifelse(is.na(tl_cred_sets_per_kb),
                                                           tr_cred_sets_per_kb,
                                                           ifelse(tl_cred_sets_per_kb >= tr_cred_sets_per_kb,
                                                                  tl_cred_sets_per_kb,
                                                                  tr_cred_sets_per_kb
                                                           )
                                                    )
                                             ), 
                                             ifelse(tl_prob > tr_prob, 
                                                    tl_cred_sets_per_kb, 
                                                    tr_cred_sets_per_kb
                                             )
                                      ), 
                                      tl_cred_sets_per_kb
    ),
    # Then classify based on the selected count of cred_sets
    cred_set_category = case_when(
      is.na(selected_cred_set_per_kb) ~ "NA",
      selected_cred_set_per_kb <= 0.059062 ~ "≤0.06",
      selected_cred_set_per_kb > 0.059062 & selected_cred_set_per_kb <= 0.122556 ~ ">0.06",
      selected_cred_set_per_kb > 0.122556 & selected_cred_set_per_kb <= 0.272276 ~ ">0.12",
      selected_cred_set_per_kb >= 0.272276 ~ "≥0.3"
    )
  ) %>%
  select(highest_prob,tl_prob,tr_prob,
         tl_sum_pip,tr_sum_pip,selected_pip,pip_category,
         tl_cred_sets_per_kb,tr_cred_sets_per_kb,selected_cred_set_per_kb,cred_set_category) #%>%
  #filter(!is.na(selected_pip))
summary(filterd_pip_data$selected_pip)
summary(filterd_pip_data$selected_cred_set_per_kb)

table(filterd_pip_data$pip_category)
table(filterd_pip_data$cred_set_category)
filterd_pip_data$pip_category <- factor(filterd_pip_data$pip_category,
                                        levels = c("≤0.07",">0.07","≥0.15"))
filterd_pip_data$cred_set_category <- factor(filterd_pip_data$cred_set_category,
                                             levels = c("NA","≤0.06",">0.06",">0.12","≥0.3"))

# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filterd_pip_data %>%
  count(pip_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$pip_category, "\n(n=", legend_data$n, ")"),
  legend_data$pip_category
)

# Histogram in log scale
ggplot(filterd_pip_data, aes(x = selected_pip)) +
  geom_histogram(bins = 100) +
  scale_x_log10(breaks=c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100)) +
  labs(
    title = "GWAS PIP",
    x = "Sum PIP (log scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Histogram normal scale
ggplot(filterd_pip_data, aes(x = selected_pip)) +
  geom_histogram(bins = 100) +
  #coord_cartesian(c(0,6)) +
  labs(
    title = "GWAS PIP",
    x = "Sum PIP (linear scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("≤0.07",">0.07"),
                       c("≤0.07","≥0.15"))

# Compute KS stats

# --- Step 1: Pre-calculate statistics for labels ---
reference_group <- "≤0.07"
comparison_groups <- c(">0.07", "≥0.15")

# Extract the data for the reference group
reference_data <- filterd_pip_data %>%
  filter(pip_category == reference_group) %>%
  pull(highest_prob)

# Calculate KS statistic for each comparison group
stats_list <- lapply(comparison_groups, function(group) {
  # Extract data for the current comparison group
  comparison_data <- filterd_pip_data %>%
    filter(pip_category == group) %>%
    pull(highest_prob)
  
  # Perform the KS test
  ks_result <- ks.test(reference_data, comparison_data)
  
  # Return a data frame with the necessary info for plotting
  data.frame(
    pip_category = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

# Combine the list of data frames into a single data frame
stats_labels <- do.call(rbind, stats_list)

# Define the y-position for the labels (adjust as needed)
stats_labels$y_position <- 1.1

plot_modified <- ggplot(filterd_pip_data, aes(x = pip_category, y = highest_prob, fill = pip_category)) +
  #geom_violin(scale = "area") +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS Sum PIP (kb)",
    x = "Sum PIP per kb",
    y = "lncRNA Probability",
    caption = "*p-val<5e-8"
  ) +
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    text = element_text(size = 36),
    plot.title = element_text(size = 46, hjust = 0.5),
    axis.text.x = element_text(hjust = 0.5, size = 30),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 44),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
    plot.caption = element_text(size = 20, hjust = -0.12, face = "bold.italic", color = "grey40")
  )

# Add the statistical comparison layer
plot_with_text_stats <- plot_modified +
  
  # Add the pre-calculated stats as text
  geom_text(
    data = stats_labels,
    aes(x = pip_category, y = y_position, label = label),
    inherit.aes = FALSE, 
    size = 9.5,
    color = "black",
    fontface = "bold"
  ) +
  
  # Adjust y-axis limits to ensure text is visible
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.15), clip = "off")
#plot_with_stats <- plot_modified +
#  geom_signif(
#    comparisons = my_comparisons,
#    test = "ks_test_custom",
#    step_increase = 0.2,
#    textsize = 9.5,
#    tip_length = 0.01,
#    y_position = 1.2
#  ) +
  # KEY CHANGE: Set explicit breaks for the y-axis and use coord_cartesian to set the visual range.
  # This prevents nonsensical axis ticks (e.g., > 1) while keeping room for annotations.
#  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
#  coord_cartesian(ylim = c(0, 1.7), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations


# Display the final plot
print(plot_with_text_stats)
#print(plot_with_stats)



########################################
# - Number of credible sets analysis - #
# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filterd_pip_data %>%
  count(cred_set_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$cred_set_category, "\n(n=", legend_data$n, ")"),
  legend_data$cred_set_category
)

# Histogram in log scale
ggplot(filterd_pip_data, aes(x = selected_cred_set_per_kb)) +
  geom_histogram(bins = 100) +
  scale_x_log10(breaks=c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100)) +
  labs(
    title = "GWAS Credible Sets",
    x = "Number of credible sets (log scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Histogram normal scale
ggplot(filterd_pip_data, aes(x = selected_cred_set_per_kb)) +
  geom_histogram(bins = 100) +
  #coord_cartesian(c(0,6)) +
  labs(
    title = "GWAS Credible Sets",
    x = "Number of credible sets (linear scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("NA","≤0.06"),
                       c("NA",">0.06"),
                       c("NA",">0.12"),
                       c("NA","≥0.3"))

# Compute KS stats

# --- Step 1: Pre-calculate statistics for labels ---
reference_group <- "NA"
comparison_groups <- c("≤0.06",">0.06", ">0.12", "≥0.3")

# Extract the data for the reference group
reference_data <- filterd_pip_data %>%
  filter(cred_set_category == reference_group) %>%
  pull(highest_prob)

# Calculate KS statistic for each comparison group
stats_list <- lapply(comparison_groups, function(group) {
  # Extract data for the current comparison group
  comparison_data <- filterd_pip_data %>%
    filter(cred_set_category == group) %>%
    pull(highest_prob)
  
  # Perform the KS test
  ks_result <- ks.test(reference_data, comparison_data)
  
  # Return a data frame with the necessary info for plotting
  data.frame(
    cred_set_category = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

# Combine the list of data frames into a single data frame
stats_labels <- do.call(rbind, stats_list)

# Define the y-position for the labels (adjust as needed)
stats_labels$y_position <- 1.1




plot_modified <- ggplot(filterd_pip_data, aes(x = cred_set_category, y = highest_prob, fill = cred_set_category)) +
  #geom_violin(scale = "area") +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS Credible Sets (kb)",
    #x = "Number of Credible sets per kb",
    y = "lncRNA Probability",
    caption = "*p-val<5-e8"
  ) +
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    text = element_text(size = 36),
    plot.title = element_text(size = 46, hjust = 0.5),
    axis.text.x = element_text(hjust = 0.5, size = 30),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 44),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
    plot.caption = element_text(size = 20, hjust = -0.12, face = "bold.italic", color = "grey40")
  )

# Add the statistical comparison layer
plot_with_text_stats <- plot_modified +
  
  # Add the pre-calculated stats as text
  geom_text(
    data = stats_labels,
    aes(x = cred_set_category, y = y_position, label = label),
    inherit.aes = FALSE, 
    size = 9.5,
    color = "black",
    fontface = "bold"
  ) +
  
  # Adjust y-axis limits to ensure text is visible
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.15), clip = "off")
#plot_with_stats <- plot_modified +
#  geom_signif(
#    comparisons = my_comparisons,
#    test = "ks_test_custom",
#    step_increase = 0.2,
#    textsize = 9.5,
#    tip_length = 0.01,
#    y_position = 1.2
#  ) +
#  # KEY CHANGE: Set explicit breaks for the y-axis and use coord_cartesian to set the visual range.
#  # This prevents nonsensical axis ticks (e.g., > 1) while keeping room for annotations.
#  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
#  coord_cartesian(ylim = c(0, 1.9), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations


# Display the final plot
print(plot_with_text_stats)
#print(plot_with_stats)


