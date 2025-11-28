library(dplyr)
library(ggplot2)
library(tidyr)
library(extrafont)
library(rstatix)
library(gridExtra)
library(grid)
library(ggnewscale)  
library(multcompView)

loadfonts(device = "pdf")

###################################################################
# Functions

get_significance <- function(p_value) {
  if (is.na(p_value)) return("ns")
  if (p_value < 0.001) return("***")
  else if (p_value < 0.01) return("**")
  else if (p_value < 0.05) return("*")
  else return("ns")
}


letters_from_pairwise <- function(pair_df, data_df,
                                  group_col = "Group",
                                  value_col = "value",
                                  alpha = 0.05,
                                  augment_non_sig = FALSE) {
  # group order
  grp <- data_df[[group_col]]
  if (is.factor(grp)) levs <- levels(droplevels(grp)) else levs <- unique(as.character(grp))
  if (length(levs) < 2L) return(setNames("a", levs))
  
  pcol <- if ("p.adj" %in% names(pair_df)) "p.adj"
  else if ("p_adj" %in% names(pair_df)) "p_adj"
  else stop("Adjusted p-value column (p.adj or p_adj) not found.")
  
  # compatibility matrix: TRUE = not different (can share letters)
  comp <- matrix(TRUE, nrow = length(levs), ncol = length(levs),
                 dimnames = list(levs, levs))
  for (i in seq_len(nrow(pair_df))) {
    g1 <- as.character(pair_df$group1[i]); g2 <- as.character(pair_df$group2[i])
    if (!nzchar(g1) || !nzchar(g2) || !(g1 %in% levs) || !(g2 %in% levs)) next
    p  <- suppressWarnings(as.numeric(pair_df[[pcol]][i]))
    is_ndiff <- is.na(p) || p >= alpha
    comp[g1, g2] <- is_ndiff
    comp[g2, g1] <- is_ndiff
  }
  
  # minimal CLD (greedy)
  alph <- c(letters, as.vector(outer(letters, letters, paste0)))
  out <- setNames(rep("", length(levs)), levs)
  out[levs[1]] <- alph[1]
  next_i <- 2
  
  holders_of <- function(L) names(out)[vapply(out, function(s) grepl(L, s, fixed = TRUE), logical(1))]
  safe_to_add <- function(g, L) all(comp[g, holders_of(L)])
  
  for (g in levs[-1]) {
    used_letters <- unique(unlist(strsplit(paste0(out, collapse = ""), "")))
    assigned <- character(0)
    if (length(used_letters)) {
      for (L in used_letters) {
        if (safe_to_add(g, L)) assigned <- c(assigned, L)
      }
    }
    if (length(assigned) == 0) {
      out[g] <- alph[next_i]; next_i <- next_i + 1
    } else {
      out[g] <- paste0(assigned, collapse = "")
    }
  }
  
  if (!augment_non_sig) return(out)
  
  ns_pairs <- which(comp & upper.tri(comp), arr.ind = TRUE)
  letters_of <- function(g) unlist(strsplit(out[g], ""))
  
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    if (nrow(ns_pairs) == 0) break
    for (k in seq_len(nrow(ns_pairs))) {
      g1 <- rownames(comp)[ns_pairs[k, 1]]
      g2 <- colnames(comp)[ns_pairs[k, 2]]
      if (length(intersect(letters_of(g1), letters_of(g2))) > 0) next  # already share
      
      # try reusing an existing letter safe for both
      used_letters <- unique(unlist(strsplit(paste0(out, collapse = ""), "")))
      reused <- FALSE
      if (length(used_letters)) {
        for (L in used_letters) {
          if (safe_to_add(g1, L) && safe_to_add(g2, L)) {
            out[g1] <- paste0(out[g1], L)
            out[g2] <- paste0(out[g2], L)
            changed <- TRUE; reused <- TRUE; break
          }
        }
      }
      if (reused) next
      
      # otherwise, introduce a new letter and assign to both
      newL <- alph[next_i]; next_i <- next_i + 1
      out[g1] <- paste0(out[g1], newL)
      out[g2] <- paste0(out[g2], newL)
      changed <- TRUE
    }
  }
  out
}


get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) {
    legend <- tmp$grobs[[leg]]
    return(legend)
  } else {
    return(NULL)
  }
}

###################################################################
# Read and prepare data

# Input files
annot <- read.delim("~/Desktop/Posdoc/CAU/annotations/Zpa796_secretome_metadata.tsv")
props <- read.csv("~/Desktop/Posdoc/CAU/annotations/amapec/intermediate_files/Zpa796_AF2bestmodels_amapec_seqProperties.csv")
struct_props <- read.csv("~/Desktop/Posdoc/CAU/annotations/amapec/intermediate_files/Zpa796_AF2bestmodels_amapec_structProperties.csv")

annot$ID <- as.character(annot$Protein.ID)
props$ID <- as.character(props$Name)
colnames(struct_props)
struct_props$ID <- as.character(struct_props$X)
struct_props <- struct_props %>%
  select(ID, Surface_Hydrophobicity = Surface.hydrophobicity)
colnames(struct_props)

# Define protein groups
secreted_ids <- annot %>% filter(Secretome == "Yes") %>% pull(ID)
effector_ids <- annot %>% filter(EffectorP3 != "Non-effector" & !is.na(EffectorP3)) %>% pull(ID)
amp_ids <- annot %>% filter(AM.prediction..AF2.AMAPEC. == "Antimicrobial") %>% pull(ID)

# Filter for proteins in structural clusters
cluster_data <- annot %>%
  filter(!is.na(Structural.subgraph) & grepl("G\\.", Structural.subgraph)) %>%
  select(ID, Structural.subgraph, Secretome, EffectorP3, AM.prediction..AF2.AMAPEC.) %>%
  left_join(props, by = "ID") %>%
  left_join(struct_props, by = "ID")

cat("Total proteins in structural clusters:", nrow(cluster_data), "\n\n")

num_cols <- c("Net.Charge", "Hydrophobicity", "Hydrophobic.Moment", "Length",
              "Surface_Hydrophobicity")
for (cn in intersect(num_cols, colnames(cluster_data))) {
  cluster_data[[cn]] <- suppressWarnings(as.numeric(as.character(cluster_data[[cn]])))
}

# Calculate charge density
if ("Length" %in% colnames(cluster_data) & "Net.Charge" %in% colnames(cluster_data)) {
  cluster_data$Charge.density <- ifelse(!is.na(cluster_data$Length) & cluster_data$Length > 0,
                                        cluster_data$Net.Charge / cluster_data$Length, NA_real_)
}

cluster_data <- cluster_data %>%
  mutate(
    Group = case_when(
      ID %in% effector_ids & ID %in% amp_ids ~ "Antimicrobial effectors",
      !(ID %in% effector_ids) & ID %in% amp_ids ~ "AMPs only",
      ID %in% effector_ids & !(ID %in% amp_ids) ~ "Effectors only",
      TRUE ~ "Other secreted proteins"
    ),
    Group = factor(Group, levels = c("Effectors only", "Antimicrobial effectors",
                                     "AMPs only", "Other secreted proteins"))
  )

# Cluster statistics
cluster_stats <- cluster_data %>%
  group_by(Structural.subgraph) %>%
  summarise(
    n_proteins = n(),
    mean_length = mean(Length, na.rm = TRUE),
    mean_charge = mean(Net.Charge, na.rm = TRUE),
    mean_charge_density = mean(Charge.density, na.rm = TRUE),
    mean_hydro = mean(Hydrophobicity, na.rm = TRUE),
    mean_moment = mean(Hydrophobic.Moment, na.rm = TRUE),
    mean_surf_hydro = mean(Surface_Hydrophobicity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_proteins))

cat("=== Cluster Statistics ===\n")
cat("Total clusters:", nrow(cluster_stats), "\n")
cat("Cluster size range:", min(cluster_stats$n_proteins), "-",
    max(cluster_stats$n_proteins), "\n\n")
print(head(cluster_stats, 10))


outdir <- "physicochemical_analysis"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(cluster_stats, file.path(outdir, "cluster_statistics.csv"), row.names = FALSE)

# Categorize clusters by size

# For pairwise scatter plots
large_clusters <- cluster_stats %>% filter(n_proteins >= 6) %>% pull(Structural.subgraph)
small_clusters <- cluster_stats %>% filter(n_proteins < 6) %>% pull(Structural.subgraph)

cat("\nClusters with >=6 proteins:", length(large_clusters), "\n")
cat("Clusters with <6 proteins:", length(small_clusters), "\n\n")

# For statistical tests
min_cluster_size <- 3  
clusters_for_stats <- cluster_stats %>% 
  filter(n_proteins >= min_cluster_size) %>% 
  pull(Structural.subgraph)

cat("Clusters with >=", min_cluster_size, "proteins (for statistical analysis):", 
    length(clusters_for_stats), "\n")

n_proteins_in_filtered <- cluster_data %>%
  filter(Structural.subgraph %in% clusters_for_stats) %>%
  nrow()
pct_proteins_in_filtered <- round(100 * n_proteins_in_filtered / nrow(cluster_data), 1)

cat("Proteins in these clusters:", n_proteins_in_filtered, 
    "(", pct_proteins_in_filtered, "% of proteins in network)\n\n")

# Filtered dataset
cluster_data_for_stats <- cluster_data %>%
  filter(Structural.subgraph %in% clusters_for_stats)

cluster_data <- cluster_data %>%
  mutate(
    Cluster_size_category = ifelse(Structural.subgraph %in% large_clusters,
                                   "≥6 proteins", "<6 proteins"),
    Cluster_display = ifelse(Structural.subgraph %in% large_clusters,
                             as.character(Structural.subgraph),
                             "Other clusters"),
    Cluster_display = factor(Cluster_display,
                             levels = c(large_clusters, "Other clusters"))
  )


# Plot config
group_colors <- c(
  "Effectors only"     = "#C55A11",
  "Antimicrobial effectors"         = "#9E3A8C",
  "AMPs only"= "#2171B5",
  "Other secreted proteins"         = "grey55"
)

n_large_clusters <- length(large_clusters)
shape_codes <- c(19, 17, 15, 18, 20, 8, 3, 4, 7, 9, 10, 11, 12, 13, 14, 6, 5, 2, 0)
cluster_shapes <- c(
  setNames(shape_codes[1:n_large_clusters], large_clusters),
  "Other clusters" = 1
)

base_theme <- theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    plot.title = element_text(size = 22)
  )

base_theme_violin <- theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

###################################################################
# Compare properties by protein groups

# Function for group statistics
run_group_stats <- function(df, feature_col, feature_label, file_prefix, outdir) {
  stopifnot(feature_col %in% names(df))
  
  df_feat <- df %>% 
    dplyr::select(Group, value = dplyr::all_of(feature_col)) %>% 
    tidyr::drop_na(value)
  
  df_feat <- as.data.frame(df_feat)
  df_feat$Group <- droplevels(df_feat$Group)
  
  if (nrow(df_feat) == 0 || n_distinct(df_feat$Group) < 2L) {
    message("Insufficient data for ", feature_label)
    return(invisible(NULL))
  }
  
  cat("\n=== Testing", feature_label, "===\n")
  
  # Kruskal-Wallis test 
  kw_result <- kruskal.test(value ~ Group, data = df_feat)
  
  n <- nrow(df_feat)
  k <- n_distinct(df_feat$Group)
  H <- kw_result$statistic
  epsilon_squared <- (H - k + 1) / (n - k) # effective size
  
  main_tbl <- data.frame(
    Effect = "Group",
    statistic = kw_result$statistic,
    df = kw_result$parameter,
    p = kw_result$p.value,
    epsilon_squared = epsilon_squared,
    test = "Kruskal-Wallis",
    stringsAsFactors = FALSE
  )
  
  # Dunn's test for post-hoc (BH adjusted)
  posthoc <- rstatix::dunn_test(value ~ Group, data = df_feat, p.adjust.method = "BH")
  posthoc <- as.data.frame(posthoc)
  posthoc$p.adj <- as.numeric(posthoc$p.adj)
  
  # Cohen's d for each pairwise comparison
  group_levels <- levels(df_feat$Group)
  posthoc$cohens_d <- NA
  
  for (i in 1:nrow(posthoc)) {
    g1 <- as.character(posthoc$group1[i])
    g2 <- as.character(posthoc$group2[i])
    
    if (g1 %in% group_levels && g2 %in% group_levels) {
      vals1 <- df_feat$value[df_feat$Group == g1]
      vals2 <- df_feat$value[df_feat$Group == g2]
      
      # Pooled standard deviation
      n1 <- length(vals1)
      n2 <- length(vals2)
      sd1 <- sd(vals1, na.rm = TRUE)
      sd2 <- sd(vals2, na.rm = TRUE)
      pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
      
      mean_diff <- mean(vals1, na.rm = TRUE) - mean(vals2, na.rm = TRUE)
      cohens_d <- mean_diff / pooled_sd
      
      posthoc$cohens_d[i] <- cohens_d
    }
  }
  
  descriptive_stats <- df_feat %>%
    group_by(Group) %>%
    summarise(
      N = n(),
      Mean = mean(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      SE = SD / sqrt(N),
      .groups = "drop"
    )
  descriptive_stats$feature <- feature_label
  

  letters_vec <- letters_from_pairwise(
    pair_df = posthoc, 
    data_df = df_feat,
    group_col = "Group",
    value_col = "value",
    alpha = 0.05,
    augment_non_sig = (feature_label == "Hydrophobicity")  # TRUE only for this feature
  )
  
  letters_tbl <- data.frame(
    group = names(letters_vec),
    letter = unname(letters_vec),
    stringsAsFactors = FALSE
  )
  
  descriptive_stats <- descriptive_stats %>%
    left_join(letters_tbl, by = c("Group" = "group"))
  
  main_tbl$feature <- feature_label
  posthoc$feature  <- feature_label
  letters_tbl$feature <- feature_label
  
  write.csv(main_tbl, file.path(outdir, paste0(file_prefix, "_main_test.csv")), row.names = FALSE)
  write.csv(posthoc, file.path(outdir, paste0(file_prefix, "_posthoc.csv")), row.names = FALSE)
  write.csv(descriptive_stats, file.path(outdir, paste0(file_prefix, "_descriptive_stats.csv")), row.names = FALSE)
  
  cat("p-value:", kw_result$p.value, "\n")
  cat("Overall effect size (ε²):", round(epsilon_squared, 4), "\n")
  
  invisible(list(main = main_tbl, posthoc = posthoc, letters = letters_tbl, descriptive = descriptive_stats))
}


make_violin_plot <- function(df, feature_col, feature_label, letters_df,
                             y_label, file_prefix, outdir) {
  
  df_plot <- df %>% 
    dplyr::select(Group, value = dplyr::all_of(feature_col)) %>% 
    tidyr::drop_na(value)
  
  df_plot$Group <- droplevels(df_plot$Group)
  
  # Letter position
  letter_pos <- df_plot %>%
    group_by(Group) %>%
    summarise(max_val = max(value, na.rm = TRUE), .groups = "drop") %>%
    left_join(letters_df, by = c("Group" = "group"))
  
  global_max <- max(df_plot$value, na.rm = TRUE)
  y_position_letters <- global_max * 1.5
  
  bar_positions <- data.frame(
    Group = letter_pos$Group,
    x_start = as.numeric(letter_pos$Group) - 0.4,
    x_end = as.numeric(letter_pos$Group) + 0.4,
    y = y_position_letters - (global_max * 0.01)
  )
  
  p <- ggplot(df_plot, aes(x = Group, y = value, fill = Group, color = Group)) +
    geom_violin(trim = FALSE, alpha = 0.4, aes(color = Group)) +
    geom_boxplot(width = 0.1, aes(fill = Group), color = "black", outlier.shape = NA) +
    geom_segment(data = bar_positions,
                 aes(x = x_start, xend = x_end, y = y, yend = y),
                 color = "black", inherit.aes = FALSE) +
    geom_text(data = letter_pos,
              aes(x = Group, y = y_position_letters, label = letter),
              size = 6, color = "black", inherit.aes = FALSE) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(x = NULL, y = y_label) +
    base_theme_violin +
    expand_limits(y = y_position_letters)  
  ggsave(file.path(outdir, paste0(file_prefix, "_violin.pdf")),
         plot = p, width = 10, height = 6, device = cairo_pdf)
  
  return(p)
}

# Run statistics and create violin plots
stats_length <- run_group_stats(cluster_data, "Length", "Protein length", "Length", outdir)
stats_charge <- run_group_stats(cluster_data, "Net.Charge", "Net charge", "Net_Charge", outdir)
stats_charge_density <- run_group_stats(cluster_data, "Charge.density", "Charge density", "Charge_density", outdir)
stats_muH <- run_group_stats(cluster_data, "Hydrophobic.Moment", "Hydrophobic moment", "Hydrophobic_Moment", outdir)
stats_hydro <- run_group_stats(cluster_data, "Hydrophobicity", "Hydrophobicity", "Hydrophobicity", outdir)
stats_surf_hydro <- run_group_stats(cluster_data, "Surface_Hydrophobicity", "Surface hydrophobicity", "Surface_Hydrophobicity", outdir)

plot_length <- if (!is.null(stats_length)) make_violin_plot(cluster_data, "Length", "Protein length", stats_length$letters, "Mature protein length", "Length", outdir)
plot_charge <- if (!is.null(stats_charge)) make_violin_plot(cluster_data, "Net.Charge", "Net charge", stats_charge$letters, "Net charge", "Net_Charge", outdir)
plot_charge_density <- if (!is.null(stats_charge_density)) make_violin_plot(cluster_data, "Charge.density", "Charge density", stats_charge_density$letters, "Charge density", "Charge_density", outdir)
plot_muH <- if (!is.null(stats_muH)) make_violin_plot(cluster_data, "Hydrophobic.Moment", "Hydrophobic moment", stats_muH$letters, "Hydrophobic moment", "Hydrophobic_Moment", outdir)
plot_hydro <- if (!is.null(stats_hydro)) make_violin_plot(cluster_data, "Hydrophobicity", "Hydrophobicity", stats_hydro$letters, "Mean hydrophobicity", "Hydrophobicity", outdir)
plot_surf_hydro <- if (!is.null(stats_surf_hydro)) make_violin_plot(cluster_data, "Surface_Hydrophobicity", "Surface hydrophobicity", stats_surf_hydro$letters, "Surface hydrophobicity", "Surface_Hydrophobicity", outdir)

plot_length
plot_charge 
plot_charge_density
plot_muH
plot_hydro
plot_surf_hydro

                      
###################################################################
# Scatter plots 

# Mature protein lentgh vs property
plot_data_1 <- cluster_data %>% filter(!is.na(Length) & !is.na(Net.Charge))
p1 <- ggplot(plot_data_1, aes(x = Length, y = Net.Charge, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Length, y = Net.Charge)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Mature protein length", y = "Net charge") +
  base_theme
ggsave(file.path(outdir, "length_vs_charge_by_cluster.pdf"),
       plot = p1, width = 12, height = 8, device = cairo_pdf)
p1

plot_data_1b <- cluster_data %>% filter(!is.na(Length) & !is.na(Charge.density))
p1b <- ggplot(plot_data_1b, aes(x = Length, y = Charge.density, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Length, y = Charge.density)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Mature protein length", y = "Charge density") +
  base_theme
ggsave(file.path(outdir, "length_vs_charge_density_by_cluster.pdf"),
       plot = p1b, width = 12, height = 8, device = cairo_pdf)
p1b

plot_data_2 <- cluster_data %>% filter(!is.na(Length) & !is.na(Hydrophobicity))
p2 <- ggplot(plot_data_2, aes(x = Length, y = Hydrophobicity, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Length, y = Hydrophobicity)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Mature protein length", y = "Mean hydrophobicity") +
  base_theme
ggsave(file.path(outdir, "length_vs_hydrophobicity_by_cluster.pdf"),
       plot = p2, width = 12, height = 8, device = cairo_pdf)

p2

plot_data_2c <- cluster_data %>% filter(!is.na(Length) & !is.na(Surface_Hydrophobicity))
p2c <- ggplot(plot_data_2c, aes(x = Length, y = Surface_Hydrophobicity, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Length, y = Surface_Hydrophobicity)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Mature protein length", y = "Surface hydrophobicity") +
  base_theme
ggsave(file.path(outdir, "length_vs_surface_hydrophobicity_by_cluster.pdf"),
       plot = p2c, width = 12, height = 8, device = cairo_pdf)
p2c

plot_data_3 <- cluster_data %>% filter(!is.na(Length) & !is.na(Hydrophobic.Moment))
p3 <- ggplot(plot_data_3, aes(x = Length, y = Hydrophobic.Moment, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Length, y = Hydrophobic.Moment)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Mature protein length", y = "Hydrophobic moment") +
  base_theme
ggsave(file.path(outdir, "length_vs_moment_by_cluster.pdf"),
       plot = p3, width = 12, height = 8, device = cairo_pdf)
p3


# Property vs property (color by group, shape by cluster)

# Net charge vs mean hydrophobicity
plot_data_4 <- cluster_data %>% filter(!is.na(Net.Charge) & !is.na(Hydrophobicity))
p4 <- ggplot(plot_data_4, aes(x = Net.Charge, y = Hydrophobicity, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Net.Charge, y = Hydrophobicity)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Net charge", y = "Mean hydrophobicity") +
  base_theme
ggsave(file.path(outdir, "charge_vs_hydrophobicity_by_cluster.pdf"),
       plot = p4, width = 12, height = 8, device = cairo_pdf)
p4

# Charge density vs mean hydrophobicity
plot_data_4_cd <- cluster_data %>% filter(!is.na(Charge.density) & !is.na(Hydrophobicity))
p4_cd <- ggplot(plot_data_4_cd, aes(x = Charge.density, y = Hydrophobicity, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Charge.density, y = Hydrophobicity)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Charge density", y = "Mean hydrophobicity") +
  base_theme
ggsave(file.path(outdir, "charge_density_vs_hydrophobicity_by_cluster.pdf"),
       plot = p4_cd, width = 12, height = 8, device = cairo_pdf)
p4_cd 

# Net charge vs surface hydrophobicity
plot_data_4c <- cluster_data %>% filter(!is.na(Net.Charge) & !is.na(Surface_Hydrophobicity))
p4c <- ggplot(plot_data_4c, aes(x = Net.Charge, y = Surface_Hydrophobicity, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Net.Charge, y = Surface_Hydrophobicity)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Net charge", y = "Surface hydrophobicity") +
  base_theme
ggsave(file.path(outdir, "charge_vs_surface_hydrophobicity_by_cluster.pdf"),
       plot = p4c, width = 12, height = 8, device = cairo_pdf)
p4c

# Charge density vs surface hydrophobicity
plot_data_4c_cd <- cluster_data %>% filter(!is.na(Charge.density) & !is.na(Surface_Hydrophobicity))
p4c_cd <- ggplot(plot_data_4c_cd, aes(x = Charge.density, y = Surface_Hydrophobicity, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Charge.density, y = Surface_Hydrophobicity)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Charge density", y = "Surface hydrophobicity") +
  base_theme
ggsave(file.path(outdir, "charge_density_vs_surface_hydrophobicity_by_cluster.pdf"),
       plot = p4c_cd, width = 12, height = 8, device = cairo_pdf)
p4c_cd 

# Net charge vs hydrophobic moment
plot_data_5 <- cluster_data %>% filter(!is.na(Net.Charge) & !is.na(Hydrophobic.Moment))
p5 <- ggplot(plot_data_5, aes(x = Net.Charge, y = Hydrophobic.Moment, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Net.Charge, y = Hydrophobic.Moment)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Net charge", y = "Hydrophobic moment") +
  base_theme
ggsave(file.path(outdir, "charge_vs_moment_by_cluster.pdf"),
       plot = p5, width = 12, height = 8, device = cairo_pdf)
p5

# Charge density vs hydrophobic moment
plot_data_5_cd <- cluster_data %>% filter(!is.na(Charge.density) & !is.na(Hydrophobic.Moment))
p5_cd <- ggplot(plot_data_5_cd, aes(x = Charge.density, y = Hydrophobic.Moment, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Charge.density, y = Hydrophobic.Moment)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Charge density", y = "Hydrophobic moment") +
  base_theme
ggsave(file.path(outdir, "charge_density_vs_moment_by_cluster.pdf"),
       plot = p5_cd, width = 12, height = 8, device = cairo_pdf)
p5_cd 

# Net charge vs charge density
plot_data_nc_cd <- cluster_data %>% filter(!is.na(Net.Charge) & !is.na(Charge.density))
p_nc_cd <- ggplot(plot_data_nc_cd, aes(x = Net.Charge, y = Charge.density, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Net.Charge, y = Charge.density)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Net charge", y = "Charge density") +
  base_theme
ggsave(file.path(outdir, "netcharge_vs_charge_density_by_cluster.pdf"),
       plot = p_nc_cd, width = 12, height = 8, device = cairo_pdf)
p_nc_cd 

# Mean hydrophobicity vs hydrophobic moment
plot_data_6 <- cluster_data %>% filter(!is.na(Hydrophobicity) & !is.na(Hydrophobic.Moment))
p6 <- ggplot(plot_data_6, aes(x = Hydrophobicity, y = Hydrophobic.Moment, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Hydrophobicity, y = Hydrophobic.Moment)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Mean hydrophobicity", y = "Hydrophobic moment") +
  base_theme
ggsave(file.path(outdir, "hydrophobicity_vs_moment_by_cluster.pdf"),
       plot = p6, width = 12, height = 8, device = cairo_pdf)
p6

# Surface hydrophobicity vs hydrophobic moment
plot_data_6c <- cluster_data %>% filter(!is.na(Surface_Hydrophobicity) & !is.na(Hydrophobic.Moment))
p6c <- ggplot(plot_data_6c, aes(x = Surface_Hydrophobicity, y = Hydrophobic.Moment, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Surface_Hydrophobicity, y = Hydrophobic.Moment)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Surface hydrophobicity", y = "Hydrophobic moment") +
  base_theme
ggsave(file.path(outdir, "surface_hydrophobicity_vs_moment_by_cluster.pdf"),
       plot = p6c, width = 12, height = 8, device = cairo_pdf)
p6c

# Mean hydrophobicity vs surface hydrophobicity
plot_data_8 <- cluster_data %>% filter(!is.na(Hydrophobicity) & !is.na(Surface_Hydrophobicity))
p8 <- ggplot(plot_data_8, aes(x = Hydrophobicity, y = Surface_Hydrophobicity, color = Group, shape = Cluster_display)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, mapping = aes(x = Hydrophobicity, y = Surface_Hydrophobicity)) +
  scale_color_manual(values = group_colors, name = "Secreted protein") +
  scale_shape_manual(values = cluster_shapes, name = "Structural cluster") +
  labs(x = "Mean hydrophobicity", y = "Surface hydrophobicity") +
  base_theme
ggsave(file.path(outdir, "eisenberg_vs_surface_hydrophobicity_by_cluster.pdf"),
       plot = p8, width = 12, height = 8, device = cairo_pdf)
p8


###################################################################
# Boxplots for properties by structural clusters

run_cluster_stats <- function(df, feature_col, feature_label) {
  df_feat <- df %>% 
    dplyr::select(Structural.subgraph, value = dplyr::all_of(feature_col)) %>% 
    tidyr::drop_na(value)
  
  df_feat$Structural.subgraph <- factor(df_feat$Structural.subgraph)
  
  if (nrow(df_feat) == 0 || n_distinct(df_feat$Structural.subgraph) < 2) {
    return(NULL)
  }
  
  # Kruskal-Wallis test
  kw_result <- kruskal.test(value ~ Structural.subgraph, data = df_feat)
  
  # Effect size (Epsilon squared, ε²) for overall test
  n <- nrow(df_feat)
  k <- n_distinct(df_feat$Structural.subgraph)
  H <- kw_result$statistic
  epsilon_squared <- (H - k + 1) / (n - k)
  
  cat("Overall effect size (ε²):", round(epsilon_squared, 4), "\n")
  
  if (kw_result$p.value < 0.05) {
    # Dunn's test with p-value correction
    posthoc <- rstatix::dunn_test(value ~ Structural.subgraph, data = df_feat,
                                  p.adjust.method = "BH")
    posthoc <- as.data.frame(posthoc)
    
    # Cohen's d for each pairwise comparison
    cluster_levels <- levels(df_feat$Structural.subgraph)
    posthoc$cohens_d <- NA
    
    for (i in 1:nrow(posthoc)) {
      g1 <- as.character(posthoc$group1[i])
      g2 <- as.character(posthoc$group2[i])
      
      if (g1 %in% cluster_levels && g2 %in% cluster_levels) {
        vals1 <- df_feat$value[df_feat$Structural.subgraph == g1]
        vals2 <- df_feat$value[df_feat$Structural.subgraph == g2]
        
        # Calculate pooled standard deviation
        n1 <- length(vals1)
        n2 <- length(vals2)
        sd1 <- sd(vals1, na.rm = TRUE)
        sd2 <- sd(vals2, na.rm = TRUE)
        pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
        
        # Calculate Cohen's d
        mean_diff <- mean(vals1, na.rm = TRUE) - mean(vals2, na.rm = TRUE)
        cohens_d <- mean_diff / pooled_sd
        
        posthoc$cohens_d[i] <- cohens_d
      }
    }
    
    # Count significant comparisons
    n_sig <- sum(posthoc$p.adj < 0.05, na.rm = TRUE)
    n_total <- nrow(posthoc)
    
    cat("  Significant pairwise comparisons:", n_sig, "out of", n_total, "\n")
    
    # Only calculate letters if there are significant comparisons and not too many groups
    letters_tbl <- NULL
    n_groups <- n_distinct(df_feat$Structural.subgraph)
    
    if (n_sig > 0 && n_groups <= 30) {
      tryCatch({
        letters_vec <- letters_from_pairwise(
          pair_df = posthoc,
          data_df = df_feat,
          group_col = "Structural.subgraph",
          value_col = "value",
          alpha = 0.05
        )
        
        letters_tbl <- data.frame(
          cluster = names(letters_vec),
          letter = unname(letters_vec),
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        cat("  Warning: Could not generate compact letter display\n")
      })
    } else {
      cat("  Too many groups (", n_groups, ") - skipping compact letter display\n")
    }
    
    return(list(kw = kw_result,
                epsilon_squared = epsilon_squared,   # <- renamed
                posthoc = posthoc,
                letters = letters_tbl,
                n_sig = n_sig,
                n_total = n_total,
                feature = feature_label))
  }
  
  return(list(kw = kw_result,
              epsilon_squared = epsilon_squared,   # <- renamed
              posthoc = NULL, 
              letters = NULL,
              n_sig = 0, 
              n_total = 0, 
              feature = feature_label))
}

make_cluster_boxplot <- function(df, feature_col, y_label, file_prefix, outdir,
                                 stats_result = NULL, clusters_in_stats = NULL) {
  df_plot <- df %>% filter(!is.na(!!sym(feature_col)))
  
  p <- ggplot(df_plot, aes(x = Structural.subgraph, y = !!sym(feature_col))) +
    geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
    # Points colored by protein group
    geom_point(aes(color = Group), size = 2.5, alpha = 0.7,
               position = position_jitter(width = 0, height = 0)) +
    scale_color_manual(values = group_colors, name = "Secreted protein") +
    labs(x = "", y = y_label) +
    base_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          legend.position = "bottom",
          legend.text = element_text(size = 14))
  
    ggsave(file.path(outdir, paste0(file_prefix, "_by_cluster_ALL.pdf")),
         plot = p, width = 24, height = 8, device = cairo_pdf)
  
  return(p)
}

# Store all KW results
all_kw_results <- list()

# Run statistics and create boxplots for each property
cat("\nTesting Net Charge by cluster...\n")
stats_cluster_charge <- run_cluster_stats(cluster_data_for_stats, "Net.Charge", "Net charge")
if (!is.null(stats_cluster_charge)) {
  cat("Kruskal-Wallis p-value:", stats_cluster_charge$kw$p.value, "\n")
  all_kw_results$Net_Charge <- data.frame(
    Property = "Net charge",
    Statistic = stats_cluster_charge$kw$statistic,
    df = stats_cluster_charge$kw$parameter,
    p_value = stats_cluster_charge$kw$p.value,
    epsilon_squared = stats_cluster_charge$epsilon_squared, 
    Significant = ifelse(stats_cluster_charge$kw$p.value < 0.05, "Yes", "No")
  )
  if (!is.null(stats_cluster_charge$posthoc)) {
    write.csv(stats_cluster_charge$posthoc,
              file.path(outdir, "netcharge_by_cluster_posthoc.csv"), row.names = FALSE)
    cat("Significant differences found. Posthoc results saved.\n")
  }
  ###if (!is.null(stats_cluster_charge$letters)) {
  ###  write.csv(stats_cluster_charge$letters,
  ###            file.path(outdir, "netcharge_by_cluster_letters.csv"), row.names = FALSE)
  ###}
}

# Plot ALL clusters, not just the ones used in stats
p14 <- make_cluster_boxplot(cluster_data, "Net.Charge", "Net charge", "netcharge",
                            outdir, stats_cluster_charge, clusters_for_stats)
p14

cat("\nTesting Charge Density by cluster...\n")
stats_cluster_charge_density <- run_cluster_stats(cluster_data_for_stats, "Charge.density", "Charge density")
if (!is.null(stats_cluster_charge_density)) {
  cat("Kruskal-Wallis p-value:", stats_cluster_charge_density$kw$p.value, "\n")
  all_kw_results$Charge_Density <- data.frame(
    Property = "Charge density",
    Statistic = stats_cluster_charge_density$kw$statistic,
    df = stats_cluster_charge_density$kw$parameter,
    p_value = stats_cluster_charge_density$kw$p.value,
    epsilon_squared = stats_cluster_charge_density$epsilon_squared,  
    Significant = ifelse(stats_cluster_charge_density$kw$p.value < 0.05, "Yes", "No")
  )
  if (!is.null(stats_cluster_charge_density$posthoc)) {
    write.csv(stats_cluster_charge_density$posthoc,
              file.path(outdir, "charge_density_by_cluster_posthoc.csv"), row.names = FALSE)
    cat("Significant differences found. Posthoc results saved.\n")
  }
  ###if (!is.null(stats_cluster_charge_density$letters)) {
  ###  write.csv(stats_cluster_charge_density$letters,
  ###            file.path(outdir, "charge_density_by_cluster_letters.csv"), row.names = FALSE)
  ###}
}

p14b <- make_cluster_boxplot(cluster_data, "Charge.density", "Charge density",
                             "charge_density", outdir, stats_cluster_charge_density, clusters_for_stats)
p14b

cat("\nTesting Mean Hydrophobicity by cluster...\n")
stats_cluster_hydro <- run_cluster_stats(cluster_data_for_stats, "Hydrophobicity", "Hydrophobicity")
if (!is.null(stats_cluster_hydro)) {
  cat("Kruskal-Wallis p-value:", stats_cluster_hydro$kw$p.value, "\n")
  all_kw_results$Hydrophobicity <- data.frame(
    Property = "Mean hydrophobicity",
    Statistic = stats_cluster_hydro$kw$statistic,
    df = stats_cluster_hydro$kw$parameter,
    p_value = stats_cluster_hydro$kw$p.value,
    epsilon_squared = stats_cluster_hydro$epsilon_squared,  
    Significant = ifelse(stats_cluster_hydro$kw$p.value < 0.05, "Yes", "No")
  )
  if (!is.null(stats_cluster_hydro$posthoc)) {
    write.csv(stats_cluster_hydro$posthoc,
              file.path(outdir, "hydrophobicity_by_cluster_posthoc.csv"), row.names = FALSE)
    cat("Significant differences found. Posthoc results saved.\n")
  }
  ###if (!is.null(stats_cluster_hydro$letters)) {
  ###  write.csv(stats_cluster_hydro$letters,
  ###            file.path(outdir, "hydrophobicity_by_cluster_letters.csv"), row.names = FALSE)
  ###}
}

p15 <- make_cluster_boxplot(cluster_data, "Hydrophobicity", "Mean hydrophobicity",
                            "hydrophobicity", outdir, stats_cluster_hydro, clusters_for_stats)
p15 

cat("\nTesting Surface Hydrophobicity by cluster...\n")
stats_cluster_surf <- run_cluster_stats(cluster_data_for_stats, "Surface_Hydrophobicity", "Surface hydrophobicity")
if (!is.null(stats_cluster_surf)) {
  cat("Kruskal-Wallis p-value:", stats_cluster_surf$kw$p.value, "\n")
  all_kw_results$Surface_Hydrophobicity <- data.frame(
    Property = "Surface hydrophobicity",
    Statistic = stats_cluster_surf$kw$statistic,
    df = stats_cluster_surf$kw$parameter,
    p_value = stats_cluster_surf$kw$p.value,
    epsilon_squared = stats_cluster_surf$epsilon_squared, 
    Significant = ifelse(stats_cluster_surf$kw$p.value < 0.05, "Yes", "No")
  )
  if (!is.null(stats_cluster_surf$posthoc)) {
    write.csv(stats_cluster_surf$posthoc,
              file.path(outdir, "surface_hydrophobicity_by_cluster_posthoc.csv"), row.names = FALSE)
    cat("Significant differences found. Posthoc results saved.\n")
  }
  ###if (!is.null(stats_cluster_surf$letters)) {
  ###  write.csv(stats_cluster_surf$letters,
  ###            file.path(outdir, "surface_hydrophobicity_by_cluster_letters.csv"), row.names = FALSE)
  ###}
}

p15c <- make_cluster_boxplot(cluster_data, "Surface_Hydrophobicity", "Surface hydrophobicity",
                             "surface_hydrophobicity", outdir, stats_cluster_surf, clusters_for_stats)
p15c

cat("\nTesting Hydrophobic Moment by cluster...\n")
stats_cluster_moment <- run_cluster_stats(cluster_data_for_stats, "Hydrophobic.Moment", "Hydrophobic moment")
if (!is.null(stats_cluster_moment)) {
  cat("Kruskal-Wallis p-value:", stats_cluster_moment$kw$p.value, "\n")
  all_kw_results$Hydrophobic_Moment <- data.frame(
    Property = "Hydrophobic moment",
    Statistic = stats_cluster_moment$kw$statistic,
    df = stats_cluster_moment$kw$parameter,
    p_value = stats_cluster_moment$kw$p.value,
    epsilon_squared = stats_cluster_moment$epsilon_squared, 
    Significant = ifelse(stats_cluster_moment$kw$p.value < 0.05, "Yes", "No")
  )
  if (!is.null(stats_cluster_moment$posthoc)) {
    write.csv(stats_cluster_moment$posthoc,
              file.path(outdir, "hydrophobic_moment_by_cluster_posthoc.csv"), row.names = FALSE)
    cat("Significant differences found. Posthoc results saved.\n")
  }
  ###if (!is.null(stats_cluster_moment$letters)) {
  ###  write.csv(stats_cluster_moment$letters,
  ###            file.path(outdir, "hydrophobic_moment_by_cluster_letters.csv"), row.names = FALSE)
  ###}
}

p16 <- make_cluster_boxplot(cluster_data, "Hydrophobic.Moment", "Hydrophobic moment",
                            "hydrophobic_moment", outdir, stats_cluster_moment, clusters_for_stats)
p16

cat("\nTesting Length by cluster...\n")
stats_cluster_length <- run_cluster_stats(cluster_data_for_stats, "Length", "Length")
if (!is.null(stats_cluster_length)) {
  cat("Kruskal-Wallis p-value:", stats_cluster_length$kw$p.value, "\n")
  all_kw_results$Length <- data.frame(
    Property = "Mature protein length",
    Statistic = stats_cluster_length$kw$statistic,
    df = stats_cluster_length$kw$parameter,
    p_value = stats_cluster_length$kw$p.value,
    epsilon_squared = stats_cluster_length$epsilon_squared,  # renamed
    Significant = ifelse(stats_cluster_length$kw$p.value < 0.05, "Yes", "No")
  )
  if (!is.null(stats_cluster_length$posthoc)) {
    write.csv(stats_cluster_length$posthoc,
              file.path(outdir, "length_by_cluster_posthoc.csv"), row.names = FALSE)
    cat("Significant differences found. Posthoc results saved.\n")
  }
  ###if (!is.null(stats_cluster_length$letters)) {
  ###  write.csv(stats_cluster_length$letters,
  ###            file.path(outdir, "length_by_cluster_letters.csv"), row.names = FALSE)
  ###}
}

p17 <- make_cluster_boxplot(cluster_data, "Length", "Mature protein length",
                            "length", outdir, stats_cluster_length, clusters_for_stats)
p17

# Combine all KW results into a summary table
if (length(all_kw_results) > 0) {
  kw_summary <- do.call(rbind, all_kw_results)
  rownames(kw_summary) <- NULL
  write.csv(kw_summary, file.path(outdir, "kruskal_wallis_summary_by_cluster.csv"), row.names = FALSE)
  cat("\n=== Kruskal-Wallis Summary (by Cluster) ===\n")
  print(kw_summary)
  cat("\nSummary saved to: kruskal_wallis_summary_by_cluster.csv\n")
}

# Plot cluster size distribution
p_size <- ggplot(cluster_stats, aes(x = reorder(Structural.subgraph, -n_proteins), y = n_proteins)) +
  geom_bar(stat = "identity", fill = "#2F79B5", alpha = 0.7) +
  geom_text(aes(label = n_proteins), vjust = -0.5, size = 3) +
  labs(x = "Structural cluster", y = "Number of proteins") +
  base_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))
ggsave(file.path(outdir, "cluster_size_distribution.pdf"),
       plot = p_size, width = 14, height = 8, device = cairo_pdf)
p_size


cat("\n\n\n")
cat("Done!!\n")


