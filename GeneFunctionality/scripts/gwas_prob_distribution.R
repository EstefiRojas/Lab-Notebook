# Load libraries
library(dplyr)
library(tidyr)
library(stringr) # For str_detect
library(ggplot2)
library(ggsignif)

# Load csv file
gwas_data <- read.csv("../data/gwas/lncrna-gwas-ranked-5-aug-2025.csv", header = TRUE) %>%
  separate(Probability_Functional, 
           into = c("tl_prob", "tr_prob"), 
           sep = "\\|",
           remove = FALSE) %>%
  mutate(tl_prob = as.numeric(tl_prob),
         tr_prob = as.numeric(tr_prob))

#####################
# Max Beta analysis #
filterd_beta_data <- gwas_data %>% 
  mutate(tl_prob = as.numeric(tl_prob),
         tr_prob = as.numeric(tr_prob),
         tl_abs_max_beta = abs(tl_max_beta),
         tr_abs_max_beta = abs(tr_max_beta),
         # First determine which beta to use based on probability
         selected_abs_beta = ifelse(!is.na(tr_prob), 
                                    ifelse(tl_prob == tr_prob, 
                                           ifelse(is.na(tr_abs_max_beta), 
                                                  tl_abs_max_beta, 
                                                  ifelse(tl_abs_max_beta >= tr_abs_max_beta,
                                                         tl_abs_max_beta,
                                                         tr_abs_max_beta
                                                         )
                                                  ), 
                                           ifelse(tl_prob >= tr_prob, 
                                                  tl_abs_max_beta, 
                                                  tr_abs_max_beta
                                                  )
                                           ), 
                                    tl_abs_max_beta
                                    ),
         # Then classify based on the selected beta
         abs_beta_ratio_category = case_when(
           selected_abs_beta < 0.3 ~ "<0.3",
           selected_abs_beta >= 0.3 & selected_abs_beta < 0.7 ~ "0.3-0.7",
           selected_abs_beta >= 0.7 ~ "≥0.7"
         )
  ) %>%
  select(highest_prob,tl_prob,tr_prob,tl_abs_max_beta,tr_abs_max_beta,selected_abs_beta,abs_beta_ratio_category) %>%
  filter(!is.na(selected_abs_beta))
summary(filterd_beta_data$selected_abs_beta)
filterd_beta_data$abs_beta_ratio_category <- factor(filterd_beta_data$abs_beta_ratio_category,
                                              levels = c("<0.3","0.3-0.7","≥0.7"))

# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filterd_beta_data %>%
  count(abs_beta_ratio_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$abs_beta_ratio_category, " (", legend_data$n, ")"),
  legend_data$abs_beta_ratio_category
)

# Histogram in log scale
ggplot(filterd_beta_data, aes(x = selected_abs_beta)) +
  geom_histogram(bins = 100) +
  scale_x_log10(breaks=c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100)) +
  labs(
    title = "GWAS Beta",
    x = "Absolute Max Beta (log scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Histogram normal scale
ggplot(filterd_beta_data, aes(x = selected_abs_beta)) +
  geom_histogram(bins = 10000) +
  coord_cartesian(c(0,6)) +
  labs(
    title = "GWAS Beta",
    x = "Absolute Max Beta (linear scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    legend.position = "none"
  )

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("<0.3","0.3-0.7"),
                       c("<0.3","≥0.7"))

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


plot_modified <- ggplot(filterd_beta_data, aes(x = abs_beta_ratio_category, y = highest_prob, fill = abs_beta_ratio_category)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.3, outliers=TRUE, na.rm = TRUE, position = position_dodge(width = 0.9), width=0.2) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS Beta",
    x = "Absolute Max Beta",
    y = "lncRNA Probability"
  ) +
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(), # Hiding minor grid lines for a cleaner look
    legend.position = "none"
  )
  
# Add the statistical comparison layer
plot_with_stats <- plot_modified +
  geom_signif(
    comparisons = my_comparisons,
    test = "ks_test_custom",
    step_increase = 0.2,
    textsize = 9.5,
    tip_length = 0.01,
    y_position = 1.2
  ) +
  # KEY CHANGE: Set explicit breaks for the y-axis and use coord_cartesian to set the visual range.
  # This prevents nonsensical axis ticks (e.g., > 1) while keeping room for annotations.
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.7), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations


# Display the final plot
print(plot_with_stats)
  
  
  
########################
# Min p-value analysis #
filterd_pval_data <- gwas_data %>% 
  mutate(tl_prob = as.numeric(tl_prob),
         tr_prob = as.numeric(tr_prob),
         # First determine which beta to use based on probability
         selected_pval = ifelse(!is.na(tr_prob), 
                                    ifelse(tl_prob == tr_prob, 
                                           ifelse(is.na(tr_min_p_val.1), 
                                                  tl_min_p_val.1, 
                                                  ifelse(tl_min_p_val.1 <= tr_min_p_val.1,
                                                         tl_min_p_val.1,
                                                         tr_min_p_val.1
                                                  )
                                           ), 
                                           ifelse(tl_prob >= tr_prob, 
                                                  tl_min_p_val.1, 
                                                  tr_min_p_val.1
                                           )
                                    ), 
                                    tl_min_p_val.1
         ),
         # Then classify based on the selected beta
         pval_category = case_when(
           selected_pval <= 5e-8 & selected_pval > 1e-15 ~ "<5e-8",
           selected_pval <= 1e-15 & selected_pval > 1e-25 ~ "<1e-15",
           selected_pval <= 1e-25 & selected_pval > 1e-100 ~ "<1e-25",
           selected_pval <= 1e-100 ~ "<1e-100"
         )
  ) %>%
  select(highest_prob,tl_prob,tr_prob,tl_min_p_val.1,tr_min_p_val.1,selected_pval,pval_category) %>%
  filter(!is.na(selected_pval))

filterd_pval_data$pval_category <- factor(filterd_pval_data$pval_category,
                                                    levels = c("<5e-8","<1e-15","<1e-25","<1e-100"))

ggplot(filterd_pval_data, aes(x=selected_pval))+
  geom_histogram(bins = 1000) +
  scale_x_log10()


# --- Prepare Legend Labels with Record Counts ---

# Calculate the number of records for each Essential_Status category.
legend_data <- filterd_pval_data %>%
  count(pval_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$pval_category, " (", legend_data$n, ")"),
  legend_data$pval_category
)

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("<5e-8","<1e-15"),
                       c("<5e-8","<1e-25"),
                       c("<5e-8","<1e-100"))

plot_modified <- ggplot(filterd_pval_data, aes(x = pval_category, y = highest_prob, fill = pval_category)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.3, outliers=TRUE, na.rm = TRUE, position = position_dodge(width = 0.9), width=0.2) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS p-value",
    x = "Min p-value",
    y = "lncRNA Probability"
  ) +
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(), # Hiding minor grid lines for a cleaner look
    legend.position = "none"
  )

# Add the statistical comparison layer
plot_with_stats <- plot_modified +
  geom_signif(
    comparisons = my_comparisons,
    test = "ks_test_custom",
    step_increase = 0.2,
    textsize = 9.5,
    tip_length = 0.01,
    y_position = 1.2
  ) +
  # KEY CHANGE: Set explicit breaks for the y-axis and use coord_cartesian to set the visual range.
  # This prevents nonsensical axis ticks (e.g., > 1) while keeping room for annotations.
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.7), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations


# Display the final plot
print(plot_with_stats)


######################
# SNP count analysis #
summary(gwas_data$tl_gwas_assoc)
summary(gwas_data$tl_max_beta)
summary(gwas_data$tl_prob)

summary(gwas_data$tr_gwas_assoc)
summary(gwas_data$tr_max_beta)
summary(gwas_data$tr_prob)
filtered_snp_count_data <- gwas_data %>%
  mutate(tl_prob = as.numeric(tl_prob),
         tr_prob = as.numeric(tr_prob),
         selected_snp_count = ifelse(!is.na(tr_prob), 
                                ifelse(tl_prob == tr_prob, 
                                       ifelse(is.na(tr_gwas_assoc), 
                                              tl_gwas_assoc, 
                                              ifelse(!is.na(tl_gwas_assoc),
                                                     ifelse(tl_gwas_assoc >= tr_gwas_assoc,
                                                            tl_gwas_assoc,
                                                            tr_gwas_assoc
                                                            ),
                                                     tr_gwas_assoc
                                              )
                                       ), 
                                       ifelse(tl_prob > tr_prob, 
                                              tl_gwas_assoc, 
                                              tr_gwas_assoc
                                       )
                                ), 
                                tl_gwas_assoc
                                ),
         max_snp_count = pmax(tl_gwas_assoc,tr_gwas_assoc,na.rm = TRUE),
         # Then classify based on the selected beta
         snp_count_category = case_when(
           selected_snp_count <= 10 ~ "≥1",
           selected_snp_count > 10 & selected_snp_count <= 100 ~ ">10",
           selected_snp_count > 100 & selected_snp_count <= 500 ~ ">100",
           #selected_snp_count >= 500 ~ ">=500",
           TRUE ~ "NA"
         )
  ) %>%
  select(highest_prob,tl_prob,tr_prob,tl_gwas_assoc,tr_gwas_assoc,selected_snp_count,snp_count_category,max_snp_count)# %>%
  #filter(!is.na(selected_snp_count))

summary(filtered_snp_count_data$selected_snp_count)
table(filtered_snp_count_data$snp_count_category)
summary(filtered_snp_count_data$max_snp_count)
table(filtered_snp_count_data$max_snp_count)

filtered_snp_count_data$snp_count_category <- factor(filtered_snp_count_data$snp_count_category,
                                          levels = c("NA","≥1",">10",">100"))

ggplot(filtered_snp_count_data, 
       aes(x=selected_snp_count)) +
  geom_histogram(bins = 100)

# Calculate the number of records for each Essential_Status category.
legend_data <- filtered_snp_count_data %>%
  count(snp_count_category)

# Create a named vector for the new labels.
# The names of the vector are the original categories (e.g., "Essential")
# The values are the new labels with counts (e.g., "Essential (n=1)")
new_labels <- setNames(
  paste0(legend_data$snp_count_category, " (", legend_data$n, ")"),
  legend_data$snp_count_category
)

# Define the pairwise comparisons to be performed
#my_comparisons <- combn(essentiality_labels, 2, simplify = FALSE)
my_comparisons <- list(c("NA","≥1"),
                       c("NA",">10"),
                       c("NA",">100"))

plot_modified <- ggplot(filtered_snp_count_data, aes(x = snp_count_category, y = highest_prob, fill = snp_count_category)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha=0.3, outliers=TRUE, na.rm = TRUE, position = position_dodge(width = 0.9), width=0.2) +
  # Use scale_fill_manual since we mapped the 'fill' aesthetic
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = new_labels) +
  # Update the labels for the new plot layout
  labs(
    title = "GWAS SNPs",
    x = "SNP count",
    y = "lncRNA Probability"
  ) +
  # Your custom theme remains the same
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 28), # Adjusted size for better readability
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle x-axis labels if they overlap
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(), # Hiding minor grid lines for a cleaner look
    legend.position = "none"
  )

# Add the statistical comparison layer
plot_with_stats <- plot_modified +
  geom_signif(
    comparisons = my_comparisons,
    test = "ks_test_custom",
    step_increase = 0.2,
    textsize = 9.5,
    tip_length = 0.01,
    y_position = 1.2
  ) +
  # KEY CHANGE: Set explicit breaks for the y-axis and use coord_cartesian to set the visual range.
  # This prevents nonsensical axis ticks (e.g., > 1) while keeping room for annotations.
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.7), clip = "off") # Use coord_cartesian to "zoom" without clipping annotations


# Display the final plot
print(plot_with_stats)
