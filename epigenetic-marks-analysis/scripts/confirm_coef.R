# Script to verify correlations reported in previous tables

colnames(lncrna_positive_feature_matrix)


# PhyloP Mammals vs Max Covariance
ggplot(protein_positive_feature_matrix,aes(x=phyloP_max_241w, y=Max_covariance)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +  # This adds the linear regression line
  labs(title = "Scatter Plot with Regression Line mRNA(+)",
       x = "phyloP_max_241w",
       y = "Covariance")

ggplot(protein_negative_feature_matrix,aes(x=phyloP_max_241w, y=Max_covariance)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # This adds the linear regression line
  labs(title = "Scatter Plot with Regression Line mRNA(-)",
       x = "phyloP_max_241w",
       y = "Covariance")
rcorr(protein_positive_feature_matrix$phyloP_max_241w,protein_positive_feature_matrix$Max_covariance,type = "spearman")
rcorr(protein_negative_feature_matrix$phyloP_max_241w,protein_negative_feature_matrix$Max_covariance,type = "spearman")




# CpG vs Methylome
ggplot(lncrna_positive_feature_matrix,aes(x=CpG_minmax, y=methylome_minmax)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") + 
  labs(title = "Scatter Plot with Regression Line (Min-Max Normalized)",
       x = "Normalized CpG (Min-Max)",
       y = "Normalized Methylome (Min-Max)")

ggplot(lncrna_positive_feature_matrix,aes(x=CpG, y=methylome)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") + 
  labs(title = "Scatter Plot with Regression Line",
       x = "CpG zscores",
       y = "Methylome zscores")

rcorr(lncrna_positive_feature_matrix$CpG,lncrna_positive_feature_matrix$methylome,type = "spearman")




mRNA_for_plot <- protein_positive_feature_matrix %>% mutate(group="mRNA(+)", Dataset="mRNA(+)")
mRNA_neg_for_plot <- protein_negative_feature_matrix %>% mutate(group="mRNA(-)", Dataset="mRNA(-)")
combined_df <- rbind(mRNA_for_plot, mRNA_neg_for_plot)
mRNA_joined_for_plot <- combined_df %>% mutate(group="mRNA")
df_for_plot <- rbind(mRNA_for_plot,mRNA_neg_for_plot,mRNA_joined_for_plot)

rcorr(protein_positive_feature_matrix$phyloP_max_241w,protein_positive_feature_matrix$Max_covariance,type = "spearman")
rcorr(protein_negative_feature_matrix$phyloP_max_241w,protein_negative_feature_matrix$Max_covariance,type = "spearman")
rcorr(combined_df$phyloP_max_241w,combined_df$Max_covariance,type = "spearman")

annotation_df <- data.frame(
  x_values = c(2, 2, 2), # X-coordinates for text in each facet
  y_values = c(325, 325, 325),  # Y-coordinates for text in each facet
  group = factor(c("mRNA(+)", "mRNA(-)", "mRNA"), levels = c("mRNA(+)", "mRNA(-)", "mRNA")), # Faceting variable
  custom_label = c("rho = 0.09", "rho = 0.8", "rho = 0.77") # Text labels
)

ggplot(df_for_plot, aes(x = phyloP_max_241w, y = Max_covariance, color = Dataset)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ group) + # Facet by the 'group' column
  geom_smooth(method = "lm", se = TRUE, color = "red") + 
  scale_color_manual(values = c("mRNA(+)" = "darkred", "mRNA(-)" = "darkblue")) +
  labs(title = "Scatter Plot mRNA") +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 26),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "right",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend text size
    plot.title = element_text(size = 30, hjust = 0.5), # Increase plot title size and center it
    strip.text = element_text(face="bold", size = 18) # Added size = 12
  ) + 
  geom_text(
    data = annotation_df, # Use the new data frame for these text elements
    aes(x = x_values, y = y_values, label = custom_label),
    color = "black",
    size = 6,
    fontface = "bold",
    show.legend = FALSE # Usually don't want separate legend for text annotations
  )



###################################



lncRNA_for_plot <- lncrna_positive_feature_matrix %>% mutate(group="lncRNA(+)", Dataset="lncRNA(+)")
lncRNA_neg_for_plot <- lncrna_negative_feature_matrix %>% mutate(group="lncRNA(-)", Dataset="lncRNA(-)")
combined_df <- rbind(lncRNA_for_plot, lncRNA_neg_for_plot)
lncRNA_joined_for_plot <- combined_df %>% mutate(group="lncRNA")
df_for_plot <- rbind(lncRNA_for_plot,lncRNA_neg_for_plot,lncRNA_joined_for_plot)

rcorr(lncrna_positive_feature_matrix$CpG,lncrna_positive_feature_matrix$methylome,type = "spearman")
rcorr(lncrna_negative_feature_matrix$CpG,lncrna_negative_feature_matrix$methylome,type = "spearman")
rcorr(combined_df$CpG,combined_df$methylome,type = "spearman")

annotation_df <- data.frame(
  x_values = c(-3, -3, -3), # X-coordinates for text in each facet
  y_values = c(110, 110, 110),  # Y-coordinates for text in each facet
  group = factor(c("lncRNA(+)", "lncRNA(-)", "lncRNA"), levels = c("lncRNA(+)", "lncRNA(-)", "lncRNA")), # Faceting variable
  custom_label = c("rho = -0.04", "rho = 0.45", "rho = 0.41") # Text labels
)

ggplot(df_for_plot, aes(x = CpG, y = methylome, color = Dataset)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ group) + # Facet by the 'group' column
  geom_smooth(method = "lm", se = TRUE, color = "red") + 
  scale_color_manual(values = c("lncRNA(+)" = "darkred", "lncRNA(-)" = "darkblue")) +
  labs(title = "Scatter Plot lncRNA") +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 26),  # Increase x-axis title size
    axis.title.y = element_text(size = 26),  # Increase y-axis title size
    legend.position = "right",
    legend.title = element_text(size = 20),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend text size
    plot.title = element_text(size = 30, hjust = 0.5), # Increase plot title size and center it
    strip.text = element_text(face="bold", size = 18) # Added size = 12
  ) + 
  geom_text(
    data = annotation_df, # Use the new data frame for these text elements
    aes(x = x_values, y = y_values, label = custom_label),
    color = "black",
    size = 6,
    fontface = "bold",
    show.legend = FALSE # Usually don't want separate legend for text annotations
  )
