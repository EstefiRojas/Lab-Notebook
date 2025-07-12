# Utility functions for data manipulation and plotting

melt_zscore_data <- function(dataframe, selected_features) {
  plot_list <- list()
  # Create an empty list to store data frames
  long_data_frames <- list()

  for (col in selected_features) {
    # Melt z-scores for violin plot
    df_with_zscores <- dataframe %>%
      pivot_longer(cols = col, names_to = "feature", values_to = "z_score") %>%
      dplyr::select(feature, z_score)

    # Store the result in the list
    long_data_frames[[col]] <- df_with_zscores
  }

  # Combine all melted data frames using bind_rows
  df_long <- dplyr::bind_rows(long_data_frames)
  return(df_long)
}

create_jitter_plot <- function(combined_df, plot_file_path) {
  p <- ggplot(combined_df, aes(x = feature, y = z_score, fill = Dataset)) +
    geom_jitter(aes(color = Dataset), width = 0.2, alpha = 0.5) +
    scale_color_manual(values = c(
      "mRNA (+)" = "blue", "lncRNA (+)" = "red", "sncRNA (+)" = "green",
      "mRNA (-)" = "lightblue", "lncRNA (-)" = "pink", "sncRNA (-)" = "lightgreen"
    )) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(
      title = "Z-scores of Epigenetic Features",
      x = "Epigenetic Feature",
      y = "Robust Z-score"
    )

  ggsave(plot_file_path, p, width = 15, height = 10)
}

create_violin_plot <- function(df_long_format, ylimit, title, violin_adjust = 1, violin_scale = "width", plot_path) {
  p <- ggplot(df_long_format, aes(x = feature, y = z_score, fill = `Gene type`)) +
    #geom_jitter() + 
    geom_violin(scale = violin_scale, na.rm = TRUE, adjust = violin_adjust, trim = FALSE) +
    geom_boxplot(alpha=0.0, outliers=TRUE, na.rm = TRUE, position = position_dodge(width = 0.9), width=0.2) +
    facet_wrap(~ feature, scales = "free") + 
    labs(title = "Conservation", x = "Feature", y = "Robust Z-score") +
    theme_minimal(base_size = 34) +
    theme(
      axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
      axis.text.y = element_text(size = 20),  # Increase y-axis text size
      axis.title.x = element_text(size = 0),  # Increase x-axis title size
      axis.title.y = element_text(size = 26),  # Increase y-axis title size
      legend.position = "none",
      legend.title = element_text(size = 20),  # Increase legend title size
      legend.text = element_text(size = 18),  # Increase legend text size
      plot.title = element_text(size = 0, hjust = 0.5), # Increase plot title size and center it
      plot.subtitle = element_blank(),
      axis.line.x = element_blank(),          # Remove x-axis line
      axis.ticks.x = element_blank(),         # Remove x-axis ticks
      panel.grid.major.x = element_blank(),   # Remove major grid lines along x-axis
      panel.grid.minor.x = element_blank(),    # Remove minor grid lines along x-axis
      strip.text = element_text(margin = ggplot2::margin(0,0,40,0))
    ) +
    scale_fill_manual(values = c("#F4A582FF","#c9e3f6FF","#D6604DFF","#56bdfcFF","#e37b88FF","#53a4f5FF")) + # colors
    coord_cartesian(ylim = ylimit) # Zoom into an specific area without removing data points
  
  ggsave(plot_path, p, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)
  return(p)
}