# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)

# Read the input files
tpm_counts <- read_delim('/Users/dalsasso/Desktop/Posdoc/CAU/People/Leon_Hofmann/Zpa796/tpm_counts_Zpa796_invitro_inplanta_DEseq2_final.csv', delim = "\t")
network_node_info <- read_delim('/Users/dalsasso/Desktop/Posdoc/CAU/annotations/Zpa796_secretome_metadata.tsv', delim = "\t")

# Preprocess names
tpm_counts <- tpm_counts %>%
  mutate(Gene_ID = gsub("\\.t1$", "", Gene_ID))

network_node_info <- network_node_info %>%
  mutate(`Protein ID` = gsub("\\.t1$", "", `Protein ID`))

# Filter proteins in a structural cluster
network_node_info <- network_node_info %>%
  filter(str_starts(`Structural subgraph`, "G\\.")) %>%
  mutate(`Structural subgraph` = as.character(`Structural subgraph`))


subgraph_ids <- unique(network_node_info$`Structural subgraph`)

peak_time_points <- data.frame(Subgraph = character(), Peak_Time = character(), stringsAsFactors = FALSE)

# Function to create plots and track peak time points
create_plots_for_subgraph <- function(subgraph_id) {
  filtered_network_node_info <- network_node_info %>%
    filter(`Structural subgraph` == subgraph_id)
  
  # Merge network info with TPM counts
  merged_data <- filtered_network_node_info %>%
    inner_join(tpm_counts, by = c("Protein ID" = "Gene_ID"))
  
  if (nrow(merged_data) == 0) {
    message("No data found for subgraph ", subgraph_id)
    return(NULL)
  }
  
  long_data <- merged_data %>%
    pivot_longer(cols = starts_with("invitro") | starts_with("dpi"),
                 names_to = "condition", values_to = "count") %>%
    mutate(count = as.numeric(count))
  
  # Calculate the average counts per time point for each gene
  gene_average_counts <- long_data %>%
    group_by(`Protein ID`, Time_Point = case_when(
      grepl("invitro", condition) ~ "axenic",
      grepl("dpi_10", condition) ~ "10",
      grepl("dpi_4", condition) ~ "4",
      grepl("dpi_7", condition) ~ "7"
    )) %>%
    summarise(average_count = mean(count, na.rm = TRUE), .groups = 'drop')
  
  # Calculate the overall average for the subgraph
  subgraph_average_counts <- gene_average_counts %>%
    group_by(Time_Point) %>%
    summarise(average_count = mean(average_count, na.rm = TRUE), .groups = 'drop')
  
  # Identify the time point with the highest expression
  peak_time <- subgraph_average_counts %>%
    filter(average_count == max(average_count, na.rm = TRUE)) %>%
    pull(Time_Point)
  

  peak_time_points <<- rbind(peak_time_points, data.frame(Subgraph = subgraph_id, Peak_Time = peak_time, stringsAsFactors = FALSE))
  
  gene_average_counts <- gene_average_counts %>%
    mutate(Time_Point = factor(Time_Point, levels = c("axenic", "4", "7", "10")))
  
  subgraph_average_counts <- subgraph_average_counts %>%
    mutate(Time_Point = factor(Time_Point, levels = c("axenic", "4", "7", "10")))
  
  # Log2 transformation
  gene_average_counts <- gene_average_counts %>%
    mutate(log2_average_count = log2(average_count + 1)) # Avoid log2(0)
  
  subgraph_average_counts <- subgraph_average_counts %>%
    mutate(log2_average_count = log2(average_count + 1)) # Avoid log2(0)
  
  # Plot
  log2meanTPM_plot <- ggplot() +
    # Plot individual gene averages as grey lines
    geom_line(data = gene_average_counts, aes(x = Time_Point, y = log2_average_count, group = `Protein ID`), color = "grey") +
    geom_point(data = gene_average_counts, aes(x = Time_Point, y = log2_average_count, group = `Protein ID`), color = "grey") +
    # Plot overall subgraph average as a red line
    geom_line(data = subgraph_average_counts, aes(x = Time_Point, y = log2_average_count, group = 1), color = "red", size = 2) +
    geom_point(data = subgraph_average_counts, aes(x = Time_Point, y = log2_average_count), color = "red", size = 2.5) +
    labs(title = paste("Subgraph", subgraph_id), x = "Time Point", y = "Log2 Mean TPM") +
    scale_x_discrete(labels = c("axenic" = "Axenic", "4" = "4 dpi", "7" = "7 dpi", "10" = "10 dpi")) +
    ylim(0, 15) +
    theme_classic() +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 22, face = "bold"), 
      axis.title = element_text(size = 20), 
      axis.text = element_text(size = 18),   
      legend.text = element_text(size = 22)
    )
  
  ggsave(filename = paste0("subgraph_", subgraph_id, "_normalizedTPM_plot_log2.png"), plot = log2meanTPM_plot, width = 4, height = 4)
}


for (subgraph_id in subgraph_ids) {
  create_plots_for_subgraph(subgraph_id)
}

peak_summary <- peak_time_points %>%
  group_by(Peak_Time) %>%
  summarise(Count = n(), Subgraphs = paste(Subgraph, collapse = ", "), .groups = "drop")

print(peak_summary)
write_csv(peak_summary, "peak_expression_summary.csv")
