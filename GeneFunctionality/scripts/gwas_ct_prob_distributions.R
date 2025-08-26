# Load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggsignif)

# Load csv file
gwas_data <- read.csv("../results/all_lncrna_gwas_dataset_11_aug_2025.csv", header = TRUE) %>%
  separate(Probability_Functional, 
           into = c("tl_prob", "tr_prob"), 
           sep = "\\|",
           remove = FALSE) %>%
  mutate(tl_prob = as.numeric(tl_prob),
         tr_prob = as.numeric(tr_prob))

#####################
# Max Beta analysis #
filterd_beta_data <- gwas_data %>% 
  mutate(# First determine which beta to use based on probability
         selected_abs_beta = ifelse(!is.na(tr_prob), 
                                    ifelse(tl_prob == tr_prob, 
                                           ifelse(is.na(tr_max_abs_beta_ct), 
                                                  tl_max_abs_beta_ct,
                                                  ifelse(is.na(tl_max_abs_beta_ct),
                                                         tr_max_abs_beta_ct,
                                                         ifelse(tl_max_abs_beta_ct >= tr_max_abs_beta_ct,
                                                                tl_max_abs_beta_ct,
                                                                tr_max_abs_beta_ct
                                                                )
                                                         )
                                                  ), 
                                           ifelse(tl_prob > tr_prob, 
                                                  tl_max_abs_beta_ct, 
                                                  tr_max_abs_beta_ct
                                                  )
                                           ), 
                                    tl_max_abs_beta_ct
                                    ),
         # Then classify based on the selected beta
         abs_beta_ratio_category = case_when(
           selected_abs_beta <= 0.031 ~ "≤0.03",
           selected_abs_beta > 0.031 & selected_abs_beta <= 0.152 ~ ">0.03",
           selected_abs_beta > 0.152 & selected_abs_beta <= 0.67 ~ ">0.15",
           selected_abs_beta >= 0.67 ~ "≥0.7"
         )
  ) %>%
  select(highest_prob,tl_prob,tr_prob,tl_max_abs_beta_ct,tr_max_abs_beta_ct,selected_abs_beta,abs_beta_ratio_category) %>%
  filter(!is.na(selected_abs_beta))
summary(filterd_beta_data$selected_abs_beta)
summary(gwas_data$tl_max_abs_beta_ct)
summary(gwas_data)
table(filterd_beta_data$abs_beta_ratio_category)
filterd_beta_data$abs_beta_ratio_category <- factor(filterd_beta_data$abs_beta_ratio_category,
                                              levels = c("≤0.03",">0.03",">0.15","≥0.7"))

# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filterd_beta_data %>%
  count(abs_beta_ratio_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$abs_beta_ratio_category, "\n(n=", legend_data$n, ")"),
  legend_data$abs_beta_ratio_category
)

# Histogram in log scale
ggplot(filterd_beta_data, aes(x = selected_abs_beta)) +
  geom_histogram(bins = 100) +
  scale_x_log10(breaks=c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100)) +
  labs(
    title = "GWAS Beta",
    x = "Max |Beta| (log scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Histogram normal scale
ggplot(filterd_beta_data, aes(x = selected_abs_beta)) +
  geom_histogram(bins = 50) +
  #coord_cartesian(c(0,100)) +
  labs(
    title = "GWAS Beta",
    x = "Max |Beta| (linear scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("≤0.03",">0.03"),
                       c("≤0.03",">0.15"),
                       c("≤0.03","≥0.7"))

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

# Define your reference group and comparison groups
reference_group <- "≤0.03"
comparison_groups <- c(">0.03", ">0.15", "≥0.7")

# Extract the data for the reference group
reference_data <- filterd_beta_data %>%
  filter(abs_beta_ratio_category == reference_group) %>%
  pull(highest_prob)

# Calculate KS statistic for each comparison group
stats_list <- lapply(comparison_groups, function(group) {
  # Extract data for the current comparison group
  comparison_data <- filterd_beta_data %>%
    filter(abs_beta_ratio_category == group) %>%
    pull(highest_prob)
  
  # Perform the KS test
  ks_result <- ks.test(reference_data, comparison_data)
  
  # Return a data frame with the necessary info for plotting
  data.frame(
    abs_beta_ratio_category = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

# Combine the list of data frames into a single data frame
stats_labels <- do.call(rbind, stats_list)

# Define the y-position for the labels (adjust as needed)
stats_labels$y_position <- 1.1 


plot_modified <- ggplot(filterd_beta_data, aes(x = abs_beta_ratio_category, y = highest_prob, fill = abs_beta_ratio_category)) +
  #geom_violin(scale = "width") +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS Max |Beta|",
    #x = "Max |Beta|",
    y = "lncRNA Probability",
    caption = "*p-val < 5e-8"
  ) +
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 34), # Adjusted size for better readability
    axis.text.x = element_text(hjust = 0.5), # Angle x-axis labels if they overlap
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(), # Hiding minor grid lines for a cleaner look
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
    aes(x = abs_beta_ratio_category, y = y_position, label = label),
    inherit.aes = FALSE, 
    size = 9.5,
    color = "black",
    fontface = "bold"
  ) +
  
  # Adjust y-axis limits to ensure text is visible
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.15), clip = "off") # Increased ylim slightly

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
  
  
  
########################
# Min p-value analysis #
filterd_pval_data <- gwas_data %>% 
  mutate(#tl_prob = as.numeric(tl_prob),
         #tr_prob = as.numeric(tr_prob),
         # First determine which p-val to use based on probability
         selected_pval = ifelse(!is.na(tr_prob), 
                                    ifelse(tl_prob == tr_prob, 
                                           ifelse(is.na(tr_min_p_value_ct), 
                                                  tl_min_p_value_ct, 
                                                  ifelse(is.na(tl_min_p_value_ct),
                                                         tr_min_p_value_ct,
                                                         ifelse(tl_min_p_value_ct <= tr_min_p_value_ct,
                                                                tl_min_p_value_ct,
                                                                tr_min_p_value_ct
                                                         )
                                                  )
                                           ), 
                                           ifelse(tl_prob >= tr_prob, 
                                                  tl_min_p_value_ct, 
                                                  tr_min_p_value_ct
                                           )
                                    ), 
                                    tl_min_p_value_ct
         ),
         # Then classify based on the selected beta
         pval_category = case_when(
           selected_pval <= 1e-7 & selected_pval > 1e-11 ~ "≤1e-7",
           selected_pval <= 1e-11 & selected_pval > 1e-17 ~ "<1e-11",
           selected_pval <= 1e-17 & selected_pval > 1e-38 ~ "<1e-17",
           selected_pval <= 1e-38 ~ "≤1e-38"
         )
  ) %>%
  select(highest_prob,tl_prob,tr_prob,tl_min_p_value_ct,tr_min_p_value_ct,selected_pval,pval_category) %>%
  filter(!is.na(selected_pval))
summary(log10(filterd_pval_data$selected_pval))
table(filterd_pval_data$pval_category)
filterd_pval_data$pval_category <- factor(filterd_pval_data$pval_category,
                                                    levels = c("≤1e-7","<1e-11","<1e-17","≤1e-38"))

ggplot(filterd_pval_data, aes(x=selected_pval))+
  geom_histogram(bins = 100) +
  scale_x_log10() +
  labs(
    title = "GWAS Min p-value",
    x = "Min p-value (log scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

ggplot(filterd_pval_data, aes(x=selected_pval))+
  geom_histogram(bins = 100) +
  #scale_x_log10() +
  labs(
    title = "GWAS Min p-value",
    x = "Min p-value (linear scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filterd_pval_data %>%
  count(pval_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$pval_category, "\n(n=", legend_data$n, ")"),
  legend_data$pval_category
)

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("≤1e-7","<1e-11"),
                       c("≤1e-7","<1e-17"),
                       c("≤1e-7","≤1e-38"))

# Compute KS stats

# --- Step 1: Pre-calculate statistics for labels ---

# Define your reference group and comparison groups
reference_group <- "≤1e-7"
comparison_groups <- c("<1e-11", "<1e-17", "≤1e-38")

# Extract the data for the reference group
reference_data <- filterd_pval_data %>%
  filter(pval_category == reference_group) %>%
  pull(highest_prob)

# Calculate KS statistic for each comparison group
stats_list <- lapply(comparison_groups, function(group) {
  # Extract data for the current comparison group
  comparison_data <- filterd_pval_data %>%
    filter(pval_category == group) %>%
    pull(highest_prob)
  
  # Perform the KS test
  ks_result <- ks.test(reference_data, comparison_data)
  
  # Return a data frame with the necessary info for plotting
  data.frame(
    pval_category = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

# Combine the list of data frames into a single data frame
stats_labels <- do.call(rbind, stats_list)

# Define the y-position for the labels (adjust as needed)
stats_labels$y_position <- 1.1 

plot_modified <- ggplot(filterd_pval_data, aes(x = pval_category, y = highest_prob, fill = pval_category)) +
  #geom_violin(scale = "area") +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS Min p-value",
    #x = "Min p-value",
    y = "lncRNA Probability",
    caption = "*p-val < 5e-8"
  ) +
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 34), # Adjusted size for better readability
    axis.text.x = element_text(hjust = 0.5), # Angle x-axis labels if they overlap
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(), # Hiding minor grid lines for a cleaner look
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.caption = element_text(size = 20, hjust = -0.12, face = "bold.italic", color = "grey40"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3)
  )

# Add the statistical comparison layer
plot_with_text_stats <- plot_modified +
  
  # Add the pre-calculated stats as text
  geom_text(
    data = stats_labels,
    aes(x = pval_category, y = y_position, label = label),
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
#  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
#  coord_cartesian(ylim = c(0, 1.7), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations


# Display the final plot
print(plot_with_text_stats)
#print(plot_with_stats)


######################
# SNP count analysis #
summary(gwas_data$tl_SNP_count_ct)
summary(gwas_data$tl_max_beta)
summary(gwas_data$tl_prob)

summary(gwas_data$tr_gwas_assoc)
summary(gwas_data$tr_max_beta)
summary(gwas_data$tr_prob)
filtered_snp_count_data <- gwas_data %>%
  mutate(#tl_prob = as.numeric(tl_prob),
         #tr_prob = as.numeric(tr_prob),
         selected_snp_count = ifelse(!is.na(tr_prob), 
                                ifelse(tl_prob == tr_prob, 
                                       ifelse(is.na(tr_SNP_count_ct), 
                                              tl_SNP_count_ct, 
                                              ifelse(!is.na(tl_SNP_count_ct),
                                                     ifelse(tl_SNP_count_ct >= tr_SNP_count_ct,
                                                            tl_SNP_count_ct,
                                                            tr_SNP_count_ct
                                                            ),
                                                     tr_SNP_count_ct
                                              )
                                       ), 
                                       ifelse(tl_prob > tr_prob, 
                                              tl_SNP_count_ct, 
                                              tr_SNP_count_ct
                                       )
                                ), 
                                tl_SNP_count_ct
                                ),
         max_snp_count = pmax(tl_SNP_count_ct,tr_SNP_count_ct,na.rm = TRUE),
         # Then classify based on the selected beta
         snp_count_category = case_when(
           selected_snp_count < 1 ~ "NA",
           selected_snp_count >= 1 & selected_snp_count <= 10 ~ "≥1",
           selected_snp_count > 10 & selected_snp_count <= 100 ~ ">10",
           #selected_snp_count > 100 & selected_snp_count <= 500 ~ ">100",
           selected_snp_count >= 100 ~ "≥100"
         ),
         tl_snp_density_per_kb = tl_SNP_count_ct / (End_Transcript_Left - Start_Transcript_Left) * 1000,
         tr_snp_density_per_kb = tr_SNP_count_ct / (End_Transcript_Right - Start_Transcript_Right) * 1000,
         selected_snp_density = ifelse(!is.na(tr_prob), 
                                       ifelse(tl_prob == tr_prob, 
                                              ifelse(is.na(tr_snp_density_per_kb), 
                                                     tl_snp_density_per_kb, 
                                                     ifelse(!is.na(tl_snp_density_per_kb),
                                                            ifelse(tl_snp_density_per_kb >= tr_snp_density_per_kb,
                                                                   tl_snp_density_per_kb,
                                                                   tr_snp_density_per_kb
                                                            ),
                                                            tr_snp_density_per_kb
                                                     )
                                              ), 
                                              ifelse(tl_prob > tr_prob, 
                                                     tl_snp_density_per_kb, 
                                                     tr_snp_density_per_kb
                                              )
                                       ), 
                                       tl_snp_density_per_kb
         ),
         max_snp_density = pmax(tl_snp_density_per_kb,tr_snp_density_per_kb,na.rm = TRUE),
         # Then classify based on the selected beta
         snp_density_category = case_when(
           is.na(selected_snp_density) ~ "NA",
           selected_snp_density == 0 ~ "NA",
           selected_snp_density > 0 & selected_snp_density <= 0.076 ~ ">0",
           selected_snp_density > 0.076 & selected_snp_density <= 0.2122 ~ ">0.08",
           selected_snp_density > 0.2122 & selected_snp_density <= 0.6169 ~ ">0.2",
           selected_snp_density >= 0.6169 ~ "≥0.6"
         ),
  ) %>%
  select(highest_prob,tl_prob,tr_prob,
         tl_SNP_count_ct,tr_SNP_count_ct,selected_snp_count,snp_count_category,max_snp_count,
         tl_snp_density_per_kb,tr_snp_density_per_kb,selected_snp_density,snp_density_category,max_snp_density)# %>%
  #filter(selected_snp_density > 0)

summary(filtered_snp_count_data$selected_snp_density)
table(filtered_snp_count_data$snp_density_category)
summary(filtered_snp_count_data$max_snp_count)
table(filtered_snp_count_data$snp_count_category)

filtered_snp_count_data$snp_count_category <- factor(filtered_snp_count_data$snp_count_category,
                                          levels = c("NA","≥1",">10","≥100"))
filtered_snp_count_data$snp_density_category <- factor(filtered_snp_count_data$snp_density_category,
                                                     levels = c("NA",">0",">0.08",">0.2","≥0.6"))

ggplot(filtered_snp_count_data, 
       aes(x=selected_snp_count)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  labs(
    title = "GWAS SNP count",
    x = "SNP count (log scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

ggplot(filtered_snp_count_data, 
       aes(x=selected_snp_count)) +
  geom_histogram(bins = 80) +
  #scale_x_log10() +
  coord_cartesian(c(0,500)) +
  labs(
    title = "GWAS SNP count",
    x = "SNP count (linear scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )


# Calculate the number of records for each category.
legend_data <- filtered_snp_count_data %>%
  count(snp_count_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$snp_count_category, "\n(", legend_data$n, ")"),
  legend_data$snp_count_category
)

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("NA","≥1"),
                       c("NA",">10"),
                       #c("NA",">100"),
                       c("NA","≥100")
                       )

# Compute KS stats

# --- Step 1: Pre-calculate statistics for labels ---
reference_group <- "NA"
comparison_groups <- c("≥1",">10", "≥100")

# Extract the data for the reference group
reference_data <- filtered_snp_count_data %>%
  filter(snp_count_category == reference_group) %>%
  pull(highest_prob)

# Calculate KS statistic for each comparison group
stats_list <- lapply(comparison_groups, function(group) {
  # Extract data for the current comparison group
  comparison_data <- filtered_snp_count_data %>%
    filter(snp_count_category == group) %>%
    pull(highest_prob)
  
  # Perform the KS test
  ks_result <- ks.test(reference_data, comparison_data)
  
  # Return a data frame with the necessary info for plotting
  data.frame(
    snp_count_category = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

# Combine the list of data frames into a single data frame
stats_labels <- do.call(rbind, stats_list)

# Define the y-position for the labels (adjust as needed)
stats_labels$y_position <- 1.1 

plot_modified <- ggplot(filtered_snp_count_data, aes(x = snp_count_category, y = highest_prob, fill = snp_count_category)) +
  #geom_violin(scale = "area") +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS SNP count",
    #x = "SNP count",
    y = "lncRNA Probability",
    caption = "*p-val<5e-8"
  ) +
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 34), # Adjusted size for better readability
    axis.text.x = element_text(hjust = 0.5), # Angle x-axis labels if they overlap
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(), # Hiding minor grid lines for a cleaner look
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.caption = element_text(size = 20, hjust = -0.12, face = "bold.italic", color = "grey40"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3)
  )

# Add the statistical comparison layer
plot_with_text_stats <- plot_modified +
  
  # Add the pre-calculated stats as text
  geom_text(
    data = stats_labels,
    aes(x = snp_count_category, y = y_position, label = label),
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
#  coord_cartesian(ylim = c(0, 2), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations


# Display the final plot
print(plot_with_text_stats)
#print(plot_with_stats)





###############
# SNP density #
# Calculate the number of records for each Essential_Status category.
legend_data <- filtered_snp_count_data %>%
  count(snp_density_category)

# Create a named vector for the new labels.
new_labels <- setNames(
  paste0(legend_data$snp_density_category, "\n(n=", legend_data$n, ")"),
  legend_data$snp_density_category
)

ggplot(filtered_snp_count_data, 
       aes(x=selected_snp_density)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  labs(
    title = "GWAS SNP density",
    x = "SNP density (log scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

ggplot(filtered_snp_count_data, 
       aes(x=selected_snp_density)) +
  geom_histogram(bins = 500) +
  coord_cartesian(c(0,100)) +
  labs(
    title = "GWAS SNP density",
    x = "SNP density (linear scale)"
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
                       c("NA",">0.08"),
                       c("NA",">0.2"),
                       c("NA","≥0.6")
)

# Compute KS stats

# --- Step 1: Pre-calculate statistics for labels ---
reference_group <- "NA"
comparison_groups <- c(">0",">0.08", ">0.2", "≥0.6")

# Extract the data for the reference group
reference_data <- filtered_snp_count_data %>%
  filter(snp_density_category == reference_group) %>%
  pull(highest_prob)

# Calculate KS statistic for each comparison group
stats_list <- lapply(comparison_groups, function(group) {
  # Extract data for the current comparison group
  comparison_data <- filtered_snp_count_data %>%
    filter(snp_density_category == group) %>%
    pull(highest_prob)
  
  # Perform the KS test
  ks_result <- ks.test(reference_data, comparison_data)
  
  # Return a data frame with the necessary info for plotting
  data.frame(
    snp_density_category = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

# Combine the list of data frames into a single data frame
stats_labels <- do.call(rbind, stats_list)

# Define the y-position for the labels (adjust as needed)
stats_labels$y_position <- 1.1 



plot_modified <- ggplot(filtered_snp_count_data, aes(x = snp_density_category, y = highest_prob, fill = snp_density_category)) +
  #geom_violin(scale = "area") +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = new_labels) +
  labs(
    title = "GWAS SNP density (kb)",
    #x = "SNP density (kb)",
    y = "lncRNA Probability",
    caption = "*p-val<5e-8"
  ) +
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
    aes(x = snp_density_category, y = y_position, label = label),
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
#  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
#  coord_cartesian(ylim = c(0, 2), clip = "off")


# Display the final plot
print(plot_with_text_stats)
#print(plot_with_stats)




#########################
# Sum Abs beta analysis #
filterd_sumbeta_data <- gwas_data %>% 
  mutate(tl_sumbeta_per_kb = (tl_sum_beta_ct / (End_Transcript_Left - Start_Transcript_Left)) * 1000,
         tr_sumbeta_per_kb = (tr_sum_beta_ct / (End_Transcript_Right - Start_Transcript_Right)) * 1000,
         selected_sumbeta = ifelse(!is.na(tr_prob), 
                           ifelse(tl_prob == tr_prob, 
                                  ifelse(is.na(tr_sumbeta_per_kb), 
                                         tl_sumbeta_per_kb, 
                                         ifelse(is.na(tl_sumbeta_per_kb),
                                                tr_sumbeta_per_kb,
                                                ifelse(tl_sumbeta_per_kb <= tr_sumbeta_per_kb,
                                                       tl_sumbeta_per_kb,
                                                       tr_sumbeta_per_kb
                                                )
                                         )
                                  ), 
                                  ifelse(tl_prob >= tr_prob, 
                                         tl_sumbeta_per_kb, 
                                         tr_sumbeta_per_kb
                                  )
                           ), 
                           tl_sumbeta_per_kb
    ),
    # Then classify based on the selected beta
    sumbeta_category = case_when(
      selected_sumbeta <= 0.00385 ~ "≤4e-3",
      selected_sumbeta > 0.00385 & selected_sumbeta < 0.01656 ~ ">4e-3",
      selected_sumbeta > 0.01656 & selected_sumbeta < 0.06478 ~ ">2e-2",
      selected_sumbeta >= 0.06478 ~ "≥6e-2"
    )
  ) %>%
  select(highest_prob,tl_prob,tr_prob,tl_sum_beta_ct,tr_sum_beta_ct,
         tl_sumbeta_per_kb,tr_sumbeta_per_kb,selected_sumbeta,sumbeta_category) %>%
  filter(!is.na(selected_sumbeta))

filterd_sumbeta_data$sumbeta_category <- factor(filterd_sumbeta_data$sumbeta_category,
                                          levels = c("≤4e-3",">4e-3",">2e-2","≥6e-2"))
summary(filterd_sumbeta_data$selected_sumbeta)
table(filterd_sumbeta_data$sumbeta_category)


ggplot(filterd_sumbeta_data, aes(x=selected_sumbeta))+
  geom_histogram(bins = 100) +
  scale_x_log10() +
  labs(
    title = "GWAS Beta",
    x = "Sum(|Beta|) per kb (log scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

ggplot(filterd_sumbeta_data, aes(x=selected_sumbeta))+
  geom_histogram(bins = 200) +
  coord_cartesian(c(0,50)) +
  labs(
    title = "GWAS Beta",
    x = "Sum(|Beta|) per kb (linear scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filterd_sumbeta_data %>%
  count(sumbeta_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$sumbeta_category, "\n(n=", legend_data$n, ")"),
  legend_data$sumbeta_category
)

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("≤4e-3",">4e-3"),
                       c("≤4e-3",">2e-2"),
                       c("≤4e-3","≥6e-2"))

# Compute KS stats

# --- Step 1: Pre-calculate statistics for labels ---
reference_group <- "≤4e-3"
comparison_groups <- c(">4e-3",">2e-2", "≥6e-2")

# Extract the data for the reference group
reference_data <- filterd_sumbeta_data %>%
  filter(sumbeta_category == reference_group) %>%
  pull(highest_prob)

# Calculate KS statistic for each comparison group
stats_list <- lapply(comparison_groups, function(group) {
  # Extract data for the current comparison group
  comparison_data <- filterd_sumbeta_data %>%
    filter(sumbeta_category == group) %>%
    pull(highest_prob)
  
  # Perform the KS test
  ks_result <- ks.test(reference_data, comparison_data)
  
  # Return a data frame with the necessary info for plotting
  data.frame(
    sumbeta_category = group,
    label = paste0("KS=", round(ks_result$statistic, 2))
  )
})

# Combine the list of data frames into a single data frame
stats_labels <- do.call(rbind, stats_list)

# Define the y-position for the labels (adjust as needed)
stats_labels$y_position <- 1.1 



plot_modified <- ggplot(filterd_sumbeta_data, aes(x = sumbeta_category, y = highest_prob, fill = sumbeta_category)) +
  #geom_violin(scale = "area") +
  geom_boxplot(linewidth = 0.9, na.rm = TRUE, outlier.shape = NA, color = "black", staplewidth = 0.5) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS Sum(|Beta|) (kb)",
    #x = "Sum(|Beta|) per kb",
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
    aes(x = sumbeta_category, y = y_position, label = label),
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

