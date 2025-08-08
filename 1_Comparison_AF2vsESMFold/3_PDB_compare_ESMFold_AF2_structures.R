library(dplyr)
library(ggplot2)
library(tidyr)
library(colorspace)
library(rstatix)
library(effsize)
library(multcompView)
library(seqinr)
library(extrafont)


font_import()  
loadfonts(device = "pdf") 



get_significance <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}


#########################################
# Load the data and get protein IDs

af2_vs_esmfold <- read.delim("Protein_structures/Protein_Data_Bank/af2_vs_esmfold_TM-align.tsv")
pdb_vs_af2 <- read.delim("Protein_structures/Protein_Data_Bank/pdb_vs_af2_TM-align.tsv")
pdb_vs_esmfold <- read.delim("Protein_structures/Protein_Data_Bank/pdb_vs_esmfold_TM-align.tsv")


af2_vs_esmfold$Protein.ID <- sub("_af2_vs_esmfold", "", af2_vs_esmfold$Protein_Name)
pdb_vs_af2$Protein.ID <- sub("_pdb_vs_af2", "", pdb_vs_af2$Protein_Name)
pdb_vs_esmfold$Protein.ID <- sub("_pdb_vs_esmfold", "", pdb_vs_esmfold$Protein_Name)


tm_data <- bind_rows(
  af2_vs_esmfold %>% select(Protein.ID, TM_Score, RMSD) %>% mutate(Method = "AF2_vs_ESMFold"),
  pdb_vs_af2 %>% select(Protein.ID, TM_Score, RMSD) %>% mutate(Method = "PDB_vs_AF2"),
  pdb_vs_esmfold %>% select(Protein.ID, TM_Score, RMSD) %>% mutate(Method = "PDB_vs_ESMFold")
)


tm_data$TM_Score <- as.numeric(as.character(tm_data$TM_Score))
tm_data$RMSD <- as.numeric(as.character(tm_data$RMSD))


tm_data$Method <- factor(tm_data$Method,
                         levels = c("AF2_vs_ESMFold", "PDB_vs_ESMFold", "PDB_vs_AF2"),
                         labels = c("AlphaFold2 vs ESMFold", "PDB vs ESMFold", "PDB vs AlphaFold2"))
colors <- c("AlphaFold2 vs ESMFold" = "grey50", "PDB vs ESMFold" = "#2F79B5", "PDB vs AlphaFold2" = "#C13639")


#########################################
#Create plots for mean TM-score and mean RMSD for all proteins

#TM-score plot
plot_mean_TMscore <- ggplot(tm_data, aes(x = Method, y = TM_Score, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.4, aes(color = Method)) +
  geom_boxplot(width = 0.1, fill = colors) +
  scale_fill_manual(values = colors) + 
  scale_color_manual(values = colors) + 
  labs(title = "",
       x = "", y = "TM-score")+
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2),  
         color = guide_legend(nrow = 2))  
  

#RMSD plot
plot_mean_RMSD <- ggplot(tm_data, aes(x = Method, y = RMSD, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.4, aes(color = Method)) +
  geom_boxplot(width = 0.1, fill = colors) +
  scale_fill_manual(values = colors) + 
  scale_color_manual(values = colors) + 
  labs(title = "",
       x = "", y = "RMSD")+
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2),  
         color = guide_legend(nrow = 2)) 

plot_mean_TMscore
plot_mean_RMSD


#########################################
#Test for normality to decide upon the statistical test

# Normality test
tm_normality <- tm_data %>%
  group_by(Method) %>%
  shapiro_test(TM_Score)
tm_normality

rmsd_normality <- tm_data %>%
  group_by(Method) %>%
  shapiro_test(RMSD)
rmsd_normality


# Decide the statistical tests based on normality
if (all(tm_normality$p > 0.05)) {
  # If TM_Score is normally distributed, use ANOVA and Tukey's HSD
  tm_anova <- anova_test(TM_Score ~ Method, data = tm_data)
  tm_posthoc <- tm_data %>%
    tukey_hsd(TM_Score ~ Method) %>%
    mutate(Significance = sapply(p.adj, get_significance))
} else {
  # If TM_Score is not normally distributed, use Kruskal-Wallis and Dunn's test
  tm_test <- kruskal_test(TM_Score ~ Method, data = tm_data)
  tm_posthoc <- tm_data %>%
    dunn_test(TM_Score ~ Method, p.adjust.method = "BH") %>%
    mutate(Significance = sapply(p.adj, get_significance))
}

if (all(rmsd_normality$p > 0.05)) {
  # If RMSD is normally distributed, use ANOVA and Tukey's HSD
  rmsd_anova <- anova_test(RMSD ~ Method, data = tm_data)
  rmsd_posthoc <- tm_data %>%
    tukey_hsd(RMSD ~ Method) %>%
    mutate(Significance = sapply(p.adj, get_significance))
} else {
  # If RMSD is not normally distributed, use Kruskal-Wallis and Dunn's test
  rmsd_test <- kruskal_test(RMSD ~ Method, data = tm_data)
  rmsd_posthoc <- tm_data %>%
    dunn_test(RMSD ~ Method, p.adjust.method = "BH") %>%
    mutate(Significance = sapply(p.adj, get_significance))
}


#########################################
# Generate letter panel from pairwise comparisons

generate_letters <- function(p_values) {
  p_matrix <- matrix(p_values, nrow = 3, ncol = 3, byrow = TRUE)
  rownames(p_matrix) <- colnames(p_matrix) <- c("AlphaFold2 vs ESMFold", "PDB vs ESMFold", "PDB vs AlphaFold2")
  letters <- multcompLetters(p_matrix, compare = "<=")$Letters
  return(letters)
}

tm_p_values <- c(NA, tm_posthoc$p.adj[1], tm_posthoc$p.adj[2],
                 tm_posthoc$p.adj[1], NA, tm_posthoc$p.adj[3],
                 tm_posthoc$p.adj[2], tm_posthoc$p.adj[3], NA)
tm_p_values

tm_letters <- generate_letters(tm_p_values)
tm_letters 


rmsd_p_values <- c(NA, rmsd_posthoc$p.adj[1], rmsd_posthoc$p.adj[2],
                   rmsd_posthoc$p.adj[1], NA, rmsd_posthoc$p.adj[3],
                   rmsd_posthoc$p.adj[2], rmsd_posthoc$p.adj[3], NA)

rmsd_letters <- generate_letters(rmsd_p_values)
rmsd_letters


# Add letters to TM-score plot
plot_mean_TMscore <- plot_mean_TMscore +
  geom_text(aes(x = 1, y = max(tm_data$TM_Score, na.rm = TRUE) + 0.12, label = tm_letters["AlphaFold2 vs ESMFold"]),
            size = 6, color = "black") +
  geom_text(aes(x = 2, y = max(tm_data$TM_Score, na.rm = TRUE) + 0.12, label = tm_letters["PDB vs ESMFold"]),
            size = 6, color = "black") +
  geom_text(aes(x = 3, y = max(tm_data$TM_Score, na.rm = TRUE) + 0.12, label = tm_letters["PDB vs AlphaFold2"]),
            size = 6, color = "black")

# Add letters to RMSD plot
plot_mean_RMSD <- plot_mean_RMSD +
  geom_text(aes(x = 1, y = max(tm_data$RMSD, na.rm = TRUE) + 1.3, label = rmsd_letters["AlphaFold2 vs ESMFold"]),
            size = 6, color = "black") +
  geom_text(aes(x = 2, y = max(tm_data$RMSD, na.rm = TRUE) + 1.3, label = rmsd_letters["PDB vs ESMFold"]),
            size = 6, color = "black") +
  geom_text(aes(x = 3, y = max(tm_data$RMSD, na.rm = TRUE) + 1.3, label = rmsd_letters["PDB vs AlphaFold2"]),
            size = 6, color = "black")


plot_mean_TMscore
plot_mean_RMSD


ggsave("Comparison_TMscore_AF2_ESMFold_PDB_Annotated.pdf", plot = plot_mean_TMscore, device = cairo_pdf, width = 5, height = 8)
ggsave("Comparison_RMSD_AF2_ESMFold_PDB_Annotated.pdf", plot = plot_mean_RMSD, device = cairo_pdf, width = 5, height = 8)


#########################################
#########################################
#Create plots for  TM-score and mean RMSD  by protein length 

library(seqinr)

# Load the FASTA file for each protein from PDB and calculate sequence lengths
fasta_file <- "Protein_structures/Protein_Data_Bank/seqs_from_PDB.fasta"
sequences <- read.fasta(file = fasta_file, seqtype = "AA", as.string = TRUE)

sequence_lengths <- data.frame(
  Protein.ID = names(sequences),
  Length = sapply(sequences, function(x) nchar(x[[1]]))
)

# Extract the first four letters from Protein.ID to match with the other data
sequence_lengths$Protein.ID <- substr(sequence_lengths$Protein.ID, 1, 4)

af2_vs_esmfold$Protein.ID <- substr(af2_vs_esmfold$Protein_Name, 1, 4)
pdb_vs_af2$Protein.ID <- substr(pdb_vs_af2$Protein_Name, 1, 4)
pdb_vs_esmfold$Protein.ID <- substr(pdb_vs_esmfold$Protein_Name, 1, 4)


tm_data <- bind_rows(
  af2_vs_esmfold %>% select(Protein.ID, TM_Score, RMSD) %>% mutate(Method = "AF2_vs_ESMFold"),
  pdb_vs_af2 %>% select(Protein.ID, TM_Score, RMSD) %>% mutate(Method = "PDB_vs_AF2"),
  pdb_vs_esmfold %>% select(Protein.ID, TM_Score, RMSD) %>% mutate(Method = "PDB_vs_ESMFold")
)

tm_data <- merge(tm_data, sequence_lengths, by = "Protein.ID")

tm_data$TM_Score <- as.numeric(as.character(tm_data$TM_Score))
tm_data$RMSD <- as.numeric(as.character(tm_data$RMSD))

tm_data$Method <- factor(tm_data$Method,
                         levels = c("AF2_vs_ESMFold", "PDB_vs_ESMFold", "PDB_vs_AF2"),
                         labels = c("AlphaFold2 vs ESMFold", "PDB vs ESMFold", "PDB vs AlphaFold2"))


colors <- c("AlphaFold2 vs ESMFold" = "grey50", "PDB vs ESMFold" = "#2F79B5", "PDB vs AlphaFold2" = "#C13639")

tm_data$Length_Category <- cut(tm_data$Length,
                               breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, Inf),
                               labels = c("1-100", "101-200", "201-300", "301-400", "401-500", "501-600", "601-700", "701-800", "801-900", "901-1000", ">1000"))


#########################################
# Create the initial plots for TM_Score and RMSD
plot_TMscore_by_length <- ggplot(tm_data, aes(x = Length_Category, y = TM_Score, fill = Method)) +
  geom_boxplot(aes(color = Method), position = position_dodge(0.8), outlier.colour = NA, alpha = 0.4) +
  geom_point(aes(color = Method), position = position_dodge(0.8), alpha = 0.4, shape = 16) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(title = "", x = "Mature protein length", y = "TM-score") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
plot_TMscore_by_length

plot_RMSD_by_length <- ggplot(tm_data, aes(x = Length_Category, y = RMSD, fill = Method)) +
  geom_boxplot(aes(color = Method), position = position_dodge(0.8), outlier.colour = NA, alpha = 0.4) +
  geom_point(aes(color = Method), position = position_dodge(0.8), alpha = 0.4, shape = 16) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(title = "", x = "Mature protein length", y = "RMSD") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

plot_TMscore_by_length
plot_RMSD_by_length


#########################################
# Filter data to include only groups with 3 or more observations
filtered_tm_data <- tm_data %>%
  group_by(Length_Category, Method) %>%
  filter(n() >= 3) %>%
  ungroup()

#########################################
# Test for normality 
TM_score_normality <- filtered_tm_data %>%
  group_by(Length_Category, Method) %>%
  shapiro_test(TM_Score)

RMSD_normality <- filtered_tm_data %>%
  group_by(Length_Category, Method) %>%
  shapiro_test(RMSD)

TM_score_normality
RMSD_normality


#########################################
# Perform test based on normality

TM_test_results <- list()
RMSD_test_results <- list()

for(length_category in unique(tm_data$Length_Category)) {
  cat("Processing Length Category:", length_category, "\n")  # Debugging info
  
  data_filtered_tm <- filter(tm_data, Length_Category == length_category)
  
  # Check if there are at least 3 observations for each method in this length category
  valid_methods <- data_filtered_tm %>%
    group_by(Method) %>%
    summarise(n = n()) %>%
    filter(n >= 3) %>%
    pull(Method)
  
  if (length(valid_methods) < 3) {
    cat("Skipping Length Category:", length_category, "- not enough data for all methods\n")
    next
  }
  
  data_filtered_tm <- filter(data_filtered_tm, Method %in% valid_methods)
  

  # Perform the appropriate test for TM_Score
  if(all(TM_score_normality$p[TM_score_normality$Length_Category == length_category] > 0.05)) {
    cat("Performing ANOVA for TM_Scores in Length Category:", length_category, "\n")
    TM_test_result <- anova_test(TM_Score ~ Method, data = data_filtered_tm)
  } else {
    cat("Performing Kruskal-Wallis for TM_Scores in Length Category:", length_category, "\n")
    TM_test_result <- kruskal_test(TM_Score ~ Method, data = data_filtered_tm)
  }
  TM_test_results[[length_category]] <- TM_test_result
  

  # Perform the appropriate test for RMSD
  if(all(RMSD_normality$p[RMSD_normality$Length_Category == length_category] > 0.05)) {
    cat("Performing ANOVA for RMSD in Length Category:", length_category, "\n")
    RMSD_test_result <- anova_test(RMSD ~ Method, data = data_filtered_tm)
  } else {
    cat("Performing Kruskal-Wallis for RMSD in Length Category:", length_category, "\n")
    RMSD_test_result <- kruskal_test(RMSD ~ Method, data = data_filtered_tm)
  }
  RMSD_test_results[[length_category]] <- RMSD_test_result
}

# Display results
cat("TM Score Test Results:\n")
print(TM_test_results)

cat("\nRMSD Test Results:\n")
print(RMSD_test_results)


TM_pairwise_comparisons <- list()
RMSD_pairwise_comparisons <- list()

# Perform pairwise comparisons for significant length categories
for(length_category in names(TM_test_results)) {
  if(TM_test_results[[length_category]]$p < 0.05) {
    TM_pairwise_comparisons[[length_category]] <- dunn_test(TM_Score ~ Method, data = filter(tm_data, Length_Category == length_category), p.adjust.method = "BH")
  } else {
    TM_pairwise_comparisons[[length_category]] <- data.frame(Length_Category = length_category, Significance = "ns")
  }
}

for(length_category in names(RMSD_test_results)) {
  if(RMSD_test_results[[length_category]]$p < 0.05) {
    RMSD_pairwise_comparisons[[length_category]] <- dunn_test(RMSD ~ Method, data = filter(tm_data, Length_Category == length_category), p.adjust.method = "BH")
  } else {
    RMSD_pairwise_comparisons[[length_category]] <- data.frame(Length_Category = length_category, Significance = "ns")
  }
}

cat("TM Score Pairwise Comparisons:\n")
print(TM_pairwise_comparisons)

cat("\nRMSD Pairwise Comparisons:\n")
print(RMSD_pairwise_comparisons)

TM_letters_by_length <- list()
RMSD_letters_by_length <- list()

#########################################

generate_letters <- function(p_values) {
  p_matrix <- matrix(p_values, nrow = 3, ncol = 3, byrow = TRUE)
  rownames(p_matrix) <- colnames(p_matrix) <- c("AlphaFold2 vs ESMFold", "PDB vs ESMFold", "PDB vs AlphaFold2")
  letters <- multcompLetters(p_matrix, compare = "<=")$Letters
  return(letters)
}

# Process TM score pairwise comparisons to generate letters
for (length_category in names(TM_pairwise_comparisons)) {
  comparison_data <- TM_pairwise_comparisons[[length_category]]
  
  if ("p.adj" %in% colnames(comparison_data)) {
    tm_p_values <- c(NA, comparison_data$p.adj[1], comparison_data$p.adj[2],
                     comparison_data$p.adj[1], NA, comparison_data$p.adj[3],
                     comparison_data$p.adj[2], comparison_data$p.adj[3], NA)
  
    TM_letters_by_length[[length_category]] <- generate_letters(tm_p_values)
  }
}

# Process RMSD pairwise comparisons to generate letters
for (length_category in names(RMSD_pairwise_comparisons)) {
  comparison_data <- RMSD_pairwise_comparisons[[length_category]]
  
  if ("p.adj" %in% colnames(comparison_data)) {
    rmsd_p_values <- c(NA, comparison_data$p.adj[1], comparison_data$p.adj[2],
                       comparison_data$p.adj[1], NA, comparison_data$p.adj[3],
                       comparison_data$p.adj[2], comparison_data$p.adj[3], NA)
    
    RMSD_letters_by_length[[length_category]] <- generate_letters(rmsd_p_values)
  }
}


#########################################
# Define the y-position for the significance annotations
y_position_tm <- max(tm_data$TM_Score, na.rm = TRUE) + 0.1
y_position_rmsd <- max(tm_data$RMSD, na.rm = TRUE) + 0.1

bar_positions_tm <- data.frame(x_start = numeric(), x_end = numeric(), y = numeric())
bar_positions_rmsd <- data.frame(x_start = numeric(), x_end = numeric(), y = numeric())

# Add letters and "ns" to TM-score plot
for (length_category in levels(tm_data$Length_Category)) {
  x_position <- as.numeric(which(levels(tm_data$Length_Category) == length_category))
  
  if (length_category %in% names(TM_letters_by_length)) {
    letters <- TM_letters_by_length[[length_category]]
    
    x_position_af2_vs_esmfold <- x_position - 0.25
    x_position_pdb_vs_esmfold <- x_position
    x_position_pdb_vs_af2 <- x_position + 0.25
    
    plot_TMscore_by_length <- plot_TMscore_by_length +
      annotate("text", x = x_position_af2_vs_esmfold, y = y_position_tm, label = letters["AlphaFold2 vs ESMFold"], size = 6, color = "black") +
      annotate("text", x = x_position_pdb_vs_esmfold, y = y_position_tm, label = letters["PDB vs ESMFold"], size = 6, color = "black") +
      annotate("text", x = x_position_pdb_vs_af2, y = y_position_tm, label = letters["PDB vs AlphaFold2"], size = 6, color = "black")
    
    bar_positions_tm <- rbind(bar_positions_tm, data.frame(x_start = x_position_af2_vs_esmfold - 0.1, x_end = x_position_pdb_vs_af2 + 0.1, y = y_position_tm - 0.02))
    
  } else if (length_category %in% names(TM_pairwise_comparisons)) {
    if (all(TM_pairwise_comparisons[[length_category]]$Significance == "ns")) {
      plot_TMscore_by_length <- plot_TMscore_by_length +
        annotate("text", x = x_position, y = y_position_tm, label = "ns", size = 6, color = "black")
      
      bar_positions_tm <- rbind(bar_positions_tm, data.frame(x_start = x_position - 0.35, x_end = x_position + 0.35, y = y_position_tm - 0.02))
    }
  }
}

plot_TMscore_by_length <- plot_TMscore_by_length +
  geom_segment(data = bar_positions_tm, aes(x = x_start, xend = x_end, y = y, yend = y), color = "black", inherit.aes = FALSE)

# Add letters and "ns" to RMSD plot 
for (length_category in levels(tm_data$Length_Category)) {
  x_position <- as.numeric(which(levels(tm_data$Length_Category) == length_category))
  
  if (length_category %in% names(RMSD_letters_by_length)) {
    letters <- RMSD_letters_by_length[[length_category]]
    
    x_position_af2_vs_esmfold <- x_position - 0.25
    x_position_pdb_vs_esmfold <- x_position
    x_position_pdb_vs_af2 <- x_position + 0.25
    
    plot_RMSD_by_length <- plot_RMSD_by_length +
      annotate("text", x = x_position_af2_vs_esmfold, y = y_position_rmsd, label = letters["AlphaFold2 vs ESMFold"], size = 6, color = "black") +
      annotate("text", x = x_position_pdb_vs_esmfold, y = y_position_rmsd, label = letters["PDB vs ESMFold"], size = 6, color = "black") +
      annotate("text", x = x_position_pdb_vs_af2, y = y_position_rmsd, label = letters["PDB vs AlphaFold2"], size = 6, color = "black")
    
    bar_positions_rmsd <- rbind(bar_positions_rmsd, data.frame(x_start = x_position_af2_vs_esmfold - 0.1, x_end = x_position_pdb_vs_af2 + 0.1 , y = y_position_rmsd - 0.20))
    
  } else if (length_category %in% names(RMSD_pairwise_comparisons)) {
    if (all(RMSD_pairwise_comparisons[[length_category]]$Significance == "ns")) {
      plot_RMSD_by_length <- plot_RMSD_by_length +
        annotate("text", x = x_position, y = y_position_rmsd, label = "ns", size = 6, color = "black")
      
      bar_positions_rmsd <- rbind(bar_positions_rmsd, data.frame(x_start = x_position - 0.35, x_end = x_position + 0.35, y = y_position_rmsd - 0.20))
    }
  }
}

plot_RMSD_by_length <- plot_RMSD_by_length +
  geom_segment(data = bar_positions_rmsd, aes(x = x_start, xend = x_end, y = y, yend = y), color = "black", inherit.aes = FALSE)

plot_TMscore_by_length
plot_RMSD_by_length


ggsave("TMscore_by_length_category.pdf", plot = plot_TMscore_by_length, device = cairo_pdf, width = 12, height = 8.5)
ggsave("RMSD_by_length_category.pdf", plot = plot_RMSD_by_length, device = cairo_pdf, width = 12, height = 8.5)






