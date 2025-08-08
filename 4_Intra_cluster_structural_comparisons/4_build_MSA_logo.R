library(ggmsa)
library(ggplot2)


protein_sequences <- "./G.12_mature_mafft-linsi_trimmed-gappyout.fasta"

# Define specific start and end positions 
start_pos <- 1  
end_pos<- 125


# Generate the plot
msa_plot <- ggmsa(protein_sequences, start = start_pos, end = end_pos, char_width = 0.5, seq_name = TRUE, border = "white") +
  geom_seqlogo() +
  geom_msaBar()
msa_plot

output_file_name <- paste0("msa_plot_residues_", start_pos, "-", end_pos, ".pdf")

ggsave(filename = output_file_name, plot = msa_plot, width = 15, height = 4, device = "pdf")
