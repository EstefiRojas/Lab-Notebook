# Install and load necessary libraries
install.packages("gwasrapidd")
library(gwasrapidd)
library(dplyr)
install.packages("tidyverse")
library(tidyverse)

# Load lncRNA data with coordinates
lncrna_coords_data <- as_tibble(read.csv("../data/datasets/ensembl_feature/lncRNAs_functionality_data_feature.csv", header = TRUE))

# Add empty column for min_pval
lncrna_coords_data |> mutate(min_pval = NA)

# Using gwasrapidd functions, query variant associations from GWAS catalog with coords from data row-wise, then get the minimum p-value for each lncRNA in the data
lncrna_coords_data <- lncrna_coords_data |>
  rowwise() |>
  mutate(min_pval = 
           tryCatch( # If no associations or variants found, return NA
             {
               min( # Obtain minimum value
                 get_associations( # Get associations
                 variant_id = get_variants( # Get variants
                   genomic_range = list(chromosome = chromosome,
                                        start = start, 
                                        end = end))
                 @variants$variant_id) # From all the fields returned by the variants query, only keep variant_id list
                 @associations$pvalue # From all the fields returned by the associations query, only keep pvalue list
                 )
             },
             error = function(e) NA_real_
           )
  )

# Save matrix with min_pval column into a csv 
write.csv(lncrna_coords_data, "../data/datasets/ensembl_feature/lncrna_GWAS_min_p_value.csv", row.names = FALSE)

# Plot density distribution
lncrna_coords_data |> 
  mutate(Type = "lncRNA") |>
  ggplot(aes(x = -log10(min_pval), fill = Type, colour = Type)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = median(-log10(lncrna_coords_data$min_pval), na.rm = TRUE)) +
  labs(title = "GWAS catalog - min P-value lncRNA") + 
  #theme_minimal() +
  scale_fill_manual(values = c("#e37b88FF")) +
  #scale_color_manual(values = c("#53a4f5FF")) +
  theme(axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32),  # Increase y-axis text size
        axis.title.x = element_text(size = 44),  # Increase x-axis title size
        axis.title.y = element_text(size = 44),  # Increase y-axis title size
        legend.position = "none",
        plot.title = element_text(size = 64, face = "italic", hjust = 0.5)  # center the title
  )
ggsave(file.path("../results/latest1000all/GWAS_min_pval", paste0("lncRNA", "_", "GWAS_min_pval", "_density.png")), scale = 3, width = 10, height = 6, bg = "white", dpi = 300)

summary(-log10(lncrna_coords_data$min_pval), na.rm = TRUE)

###########
# Extra code to analyse problematic rows:
# Investigate row 5: This coords do not return any variant.
# 2 rows have this problem
variants_data <- get_variants(genomic_range=list(
                                 chromosome=18, 
                                 start=5238059, 
                                 end=5246508)
                                 )@variants$variant_id

# Investigate row 15: this returns only one association (rs112588865) with no data in the catalog for it
# 38 rows have this problem
variants_data <- get_variants(genomic_range=list(
                                chromosome=12, 
                                start=25936683, 
                                end=25959442))@variants$variant_id

associations_data <- get_associations(variant_id = variants_data)
