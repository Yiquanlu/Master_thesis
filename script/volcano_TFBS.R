library(ggplot2)
library(dplyr)
library(ggrepel)

# List of input files, you can just change it to interneuron TF enriched binding motifs
input_files <- c('astro_GPC_chromvar_with_TF_names.tsv', 'astro_opc_chromvar_with_TF_names.tsv', 'GPC_RG_chromvar_with_TF_names.tsv', 'opc_GPC_chromvar_with_TF_names.tsv')

create_volcano_plot <- function(input_file, p_val_cutoff = 0.05, avg_diff_cutoff = 1, output_file) {
  data <- tryCatch({
    read.table(input_file, header = TRUE, sep = "\t")
  }, error = function(e) {
    message(paste("Error reading", input_file, ":", e$message))
    return(NULL)
  })
  if (is.null(data)) return()

  # Add a small constant to p-values that are 0
  pseudo_count <- min(data$p_val[data$p_val > 0]) / 10
  data$p_val[data$p_val == 0] <- pseudo_count

  data$PointColor <- with(data, ifelse(p_val < p_val_cutoff & avg_diff > avg_diff_cutoff, "red",
                                       ifelse(p_val < p_val_cutoff & avg_diff < -avg_diff_cutoff, "blue", "gray")))

  # Filter out rows with NA in p_val, avg_diff, or TF_name before plotting
  significant_data <- subset(data, p_val < p_val_cutoff & abs(avg_diff) > avg_diff_cutoff & !is.na(TF_name))

  volcano_plot <- ggplot(data, aes(x = avg_diff, y = -log10(p_val))) +
    geom_point(aes(color = PointColor), alpha = 0.5, size = 1) +
    scale_color_identity() +
    geom_text_repel(data = significant_data, aes(label = TF_name),
                    size = 2, color = "black", max.overlaps = Inf) +
    labs(title = paste("Volcano Plot for", gsub("\\.tsv$", "", basename(input_file))), x = "Average Difference", y = "-log10(p-value)") +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_vline(xintercept = c(-avg_diff_cutoff, avg_diff_cutoff), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(p_val_cutoff), col = "black", linetype = "dashed")

  ggsave(output_file, plot = volcano_plot, width = 10, height = 6, dpi = 300)
}

# Iterate over input files and create Volcano Plots
for (input_file in input_files) {
  # Generate output file name
  output_file <- gsub("with_TF_names.tsv", "volcano_plot.png", input_file)
  
  # Create Volcano Plot and save figure
  create_volcano_plot(input_file, output_file = output_file)
}
