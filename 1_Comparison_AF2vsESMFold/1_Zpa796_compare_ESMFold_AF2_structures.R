library(dplyr)
library(ggplot2)
library(tidyr)
library(colorspace)
library(rstatix)
library(effsize)
library(multcompView)
library(extrafont)


font_import() 
loadfonts(device = "pdf") 



df <- read.delim("Zpa796_secretome_metadata.tsv")

##################################################################################################
##################################################################################################
# pLDDT values

df_melted <- df %>%
  select(Protein.ID, `pLDDT..ESMFold.`, `pLDDT..AF2.`) %>%
  pivot_longer(cols = -Protein.ID, names_to = "Method", values_to = "pLDDT")

df_melted$Method <- factor(df_melted$Method,
                           levels = c("pLDDT..ESMFold.", "pLDDT..AF2."),
                           labels = c("ESMFold", "AlphaFold2"))

colors <- c("ESMFold" = "#2F79B5", "AlphaFold2" = "#C13639")


##################################################################################################
# Plot pLDDTs values as boxplot with density distribution (violin)

plot_mean_pLDDT <- ggplot(df_melted, aes(x = Method, y = pLDDT, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.4, aes(color = Method)) +
  geom_boxplot(width = 0.1, fill = colors) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_classic() +
  labs(title = "",
       x = "", y = "pLDDT") +
  theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_blank(),       
    axis.text.x = element_blank(),        
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22)
  )
plot_mean_pLDDT


# Statistical test for pLDDT
pLDDT_long <- df %>%
  select(Protein.ID, pLDDT..ESMFold., pLDDT..AF2.) %>%
  pivot_longer(cols = starts_with("pLDDT"), names_to = "Method", values_to = "pLDDT")

pLDDT_long$Method <- recode(pLDDT_long$Method,
                            "pLDDT..ESMFold." = "ESMFold",
                            "pLDDT..AF2." = "AlphaFold2")

pLDDT_normality <- pLDDT_long %>%
  group_by(Method) %>%
  shapiro_test(pLDDT)
pLDDT_normality 


# Choose test based on normality
if (any(pLDDT_normality$p < 0.05)) {
  # Non-parametric test
  pLDDT_test <- wilcox_test(pLDDT ~ Method, data = pLDDT_long)
  pLDDT_effect_size <- effsize::cohen.d(pLDDT ~ Method, data = pLDDT_long, paired = FALSE)
} else {
  # Parametric test
  pLDDT_test <- t_test(pLDDT ~ Method, data = pLDDT_long)
  pLDDT_effect_size <- effsize::cohen.d(pLDDT ~ Method, data = pLDDT_long, paired = FALSE)
}

pLDDT_test
pLDDT_effect_size


p_value <- pLDDT_test$p
significance <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns")))

max_y <- max(df_melted$pLDDT, na.rm = TRUE) 
offset <- 0.2 * max_y  

plot_mean_pLDDT <- plot_mean_pLDDT +
  geom_segment(aes(x = 1, xend = 2, y = max_y + offset, yend = max_y + offset), color = "black") +
  annotate("text", x = 1.5, y = max_y + 1.1 * offset, label = significance, size = 6, color = "black")
plot_mean_pLDDT 


##################################################################################################
# Plot pLDDTs according to protein length range


df$Length_Category <- cut(df$`Mature.protein...AA.`,
                          breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, Inf),
                          labels = c("1-100", "101-200", "201-300", "301-400", "401-500", "501-600", "601-700", "701-800", "801-900", "901-1000", ">1000"))

df_melted_length <- df %>%
  select(Protein.ID, `pLDDT..ESMFold.`, `pLDDT..AF2.`, Length_Category) %>%
  pivot_longer(cols = c(`pLDDT..ESMFold.`, `pLDDT..AF2.`), names_to = "Method", values_to = "pLDDT") %>%
  mutate(Method = factor(Method, levels = c("pLDDT..ESMFold.", "pLDDT..AF2."),
                         labels = c("ESMFold", "AlphaFold2")))

df_melted_length$Method <- factor(df_melted_length$Method,
                                  levels = c("ESMFold", "AlphaFold2"),
                                  labels = c("ESMFold", "AlphaFold2"))

# Plotting boxplots side by side for each length category 
plot_pLDDT_by_length <- ggplot(df_melted_length, aes(x = Length_Category, y = pLDDT)) +
  geom_boxplot(
    aes(fill = Method, color = Method), # Set fill and color based on Method
    position = position_dodge(0.8),
    outlier.colour = NA,  # Hide the default outliers
    alpha = 0.4  # Set transparency for the boxplot fill
  ) +
  scale_fill_manual(values = colors) + # Apply custom fill colors
  scale_color_manual(values = colors) + # Apply custom outline colors
  theme_classic() +
  labs(title = "",
       x = "Mature protein length", y = "pLDDT") +
  geom_point(
    aes(color = Method), # Add the outliers manually with matching color
    position = position_dodge(0.8),
    alpha = 0.4,
    shape = 16 #solid circle
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
plot_pLDDT_by_length 


pLDDT_long <- df %>%
  select(Protein.ID, `pLDDT..ESMFold.`, `pLDDT..AF2.`, Length_Category) %>% 
  pivot_longer(cols = -c(Protein.ID, Length_Category), names_to = "Method", values_to = "pLDDT") %>% 
  mutate(Method = recode(Method, `pLDDT..ESMFold.` = "ESMFold", `pLDDT..AF2.` = "AlphaFold2"))

# Normality test
pLDDT_normality_by_length <- pLDDT_long %>%
  group_by(Length_Category, Method) %>%
  shapiro_test(pLDDT)
pLDDT_normality_by_length

# Loop through each length category to perform tests
# Assuming non-normality and proceeding with Wilcoxon tests
pLDDT_comparisons <- list()

unique_length_categories <- unique(pLDDT_long$Length_Category)

for(length_category in unique_length_categories) {
  data_filtered <- filter(pLDDT_long, Length_Category == length_category)
  test_result <- wilcox_test(pLDDT ~ Method, data = data_filtered)
  effect_size <- effsize::cohen.d(filter(pLDDT_long, Length_Category == length_category)$pLDDT ~ filter(pLDDT_long, Length_Category == length_category)$Method, paired = FALSE)
  pLDDT_comparisons[[length_category]] <- list(test_result = test_result, effect_size = effect_size)
}

for(length_category in unique_length_categories) {
  cat("Length Category:", length_category, "\n")
  print(pLDDT_comparisons[[length_category]]$test_result)
  print(pLDDT_comparisons[[length_category]]$effect_size)
  cat("\n")
}


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

str(pLDDT_comparisons[[1]]$test_result)
pLDDT_comparisons


# Generate significance_data dataframe based on pLDDT_comparisons
significance_data <- data.frame(Length_Category = character(), Significance = character(), stringsAsFactors = FALSE)
significance_data 

for (length_category in names(pLDDT_comparisons)) {
  p_value <- pLDDT_comparisons[[length_category]]$test_result$p
  significance <- get_significance(p_value)
  significance_data <- rbind(significance_data, data.frame(Length_Category = length_category, Significance = significance))
}

y_position <- 100 

# Annotate plot_pLDDT_by_length with significance data
plot_pLDDT_by_length <- plot_pLDDT_by_length +
  geom_text(data = significance_data, aes(x = Length_Category, y = y_position, label = Significance), position = position_dodge(width = 0.75), size = 6, vjust = 0)

dodge_width <- 0.75 
segment_y_position <- 99

plot_pLDDT_by_length <- plot_pLDDT_by_length +
  geom_segment(data = significance_data, aes(x = as.numeric(as.factor(Length_Category)) - dodge_width/2, y = segment_y_position, xend = as.numeric(as.factor(Length_Category)) + dodge_width/2, yend = segment_y_position), color = "black")

plot_pLDDT_by_length


##################################################################################################
##################################################################################################
#### TM-scores

df_melted <- df %>%
  select(Protein.ID, `TM.score..AF2.vs.ESMFold.`) %>%
  pivot_longer(cols = -Protein.ID, names_to = "Method", values_to = "TM_score")

df_melted$TM_score <- as.numeric(as.character(df_melted$TM_score))


##################################################################################################
# Create boxplot with density distribution (violin) for TM-score

plot_mean_TMscore <- ggplot(df_melted, aes(x = Method, y = TM_score, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.4, fill = "grey") + 
  geom_boxplot(width = 0.1, fill = "darkgrey", color = "black") +
  theme_classic() +
  labs(title = "",
       x = "", y = "TM-score") +
  theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_blank(),       
    axis.text.x = element_blank(),        
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    legend.position = "none"
  )
plot_mean_TMscore 


##################################################################################################
# Plot TM-scores according to protein length range

df$Length_Category <- cut(df$`Protein.length...AA.`,
                          breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, Inf),
                          labels = c("1-100", "101-200", "201-300", "301-400", "401-500", "501-600", "601-700", "701-800", "801-900", "901-1000", ">1000"))


df_melted_length <- df %>%
  select(Protein.ID, `TM.score..AF2.vs.ESMFold.`, Length_Category) %>%
  pivot_longer(cols = `TM.score..AF2.vs.ESMFold.`, names_to = "Method", values_to = "TM_score") %>%
  mutate(Method = "TM-score") # Since we only have one score, we name the method directly

df_melted_length$TM_score <- as.numeric(as.character(df_melted_length$TM_score))

# Plotting boxplots for TM-scores by protein length category
plot_TMscore_by_length <- ggplot(df_melted_length, aes(x = Length_Category, y = TM_score)) +
  geom_boxplot(
    fill = "grey", color = "grey50",
    alpha = 0.4,  
    outlier.shape = NA  
  ) +
  theme_classic() +
  labs(title = "",
       x = "Mature protein length", y = "TM-score") +
  geom_point(
    color = "grey10", 
    position = position_dodge(0.8),
    alpha = 0.4,
    shape = 16  # a solid circle
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_TMscore_by_length


# Test for normality within each length category for TM-score
TMscore_normality <- df_melted_length %>%
  group_by(Length_Category) %>%
  shapiro_test(TM_score)
TMscore_normality

# If non-normality is found in any group, use Kruskal-Wallis, otherwise use ANOVA
if(any(TMscore_normality$p < 0.05, na.rm = TRUE)) {
  TMscores_comparison <- 
    kruskal_test(TM_score ~ Length_Category, data = df_melted_length)
} else {
  TMscores_comparison <- anova_test(TM_score ~ Length_Category, data = df_melted_length)
}
TMscores_comparison

# Perform the Dunn test for multiple comparisons
TM_scores_posthoc <- df_melted_length %>%
  dunn_test(TM_score ~ Length_Category, p.adjust.method = "BH")
TM_scores_posthoc 

tm_p_values <- TM_scores_posthoc$p.adj
tm_p_values 

generate_letters <- function(p_values) {
  p_matrix <- matrix(p_values, nrow = length(unique(df_melted_length$Length_Category)), ncol = length(unique(df_melted_length$Length_Category)), byrow = TRUE)
  rownames(p_matrix) <- colnames(p_matrix) <- levels(df_melted_length$Length_Category)
  letters <- multcompLetters(p_matrix, compare = "<=")$Letters
  return(letters)
}

TM_letters_by_length <- generate_letters(tm_p_values)
TM_letters_by_length


y_position_tm <- max(df_melted_length$TM_score, na.rm = TRUE) + 0.1
for (length_category in levels(df_melted_length$Length_Category)) {
  plot_TMscore_by_length <- plot_TMscore_by_length +
    annotate("text", x = which(levels(df_melted_length$Length_Category) == length_category), y = y_position_tm, label = TM_letters_by_length[length_category], size = 5, color = "black")
}


plot_TMscore_by_length


#########################################################################################
# Check plots and save them
plot_mean_pLDDT
plot_pLDDT_by_length
plot_mean_TMscore
plot_TMscore_by_length


ggsave("Zpa796_AF2-ESMFdold_plot_mean_pLDDT.pdf", plot = plot_mean_pLDDT, device = cairo_pdf, width = 8, height = 6)
ggsave("Zpa796_AF2-ESMFdold_plot_pLDDT_by_length.pdf", plot = plot_pLDDT_by_length, device = cairo_pdf, width = 8, height = 8.5)
ggsave("Zpa796_AF2-ESMFdold_plot_mean_TMscore.pdf", plot = plot_mean_TMscore, device = cairo_pdf, width = 5, height = 6)
ggsave("Zpa796_AF2-ESMFdold_plot_TMscore_by_length.pdf", plot = plot_TMscore_by_length, device = cairo_pdf, width = 8, height = 6)


########################################################################################
########################################################################################
# Compare pLDDT values from AlphaFold2 and ESMFold with a scatter plot and test the correlation of both methods

library(ggplot2)
library(plotly)

# Load the data and define columns
data <- df
print(colnames(data))

alphaFold_col <- "pLDDT..AF2."
esmFold_col <- "pLDDT..ESMFold."
protein_id_col <- "Protein.ID"  


data[[alphaFold_col]][is.na(data[[alphaFold_col]])] <- 0
data[[esmFold_col]][is.na(data[[esmFold_col]])] <- 0


data$Method <- ifelse(data[[alphaFold_col]] > 0 & data[[esmFold_col]] > 0, "Both",
                      ifelse(data[[alphaFold_col]] > 0, "AlphaFold2",
                             ifelse(data[[esmFold_col]] > 0, "ESMFold", "None")))


print(unique(data$Method))

# Filter out rows where both pLDDT values are zero (i.e., Source is "None")
data <- data[data$Method != "None", ]

# Filter data for "Both" methods
data_both <- data[data$Method == "Both", ]

# Test for normality
shapiro_esmFold <- shapiro.test(data_both[[esmFold_col]])
shapiro_alphaFold <- shapiro.test(data_both[[alphaFold_col]])
print(shapiro_esmFold)
print(shapiro_alphaFold)

# Perform correlation test and determine the method
if (shapiro_esmFold$p.value > 0.05 & shapiro_alphaFold$p.value > 0.05) {
  # Both distributions are normal
  correlation_test_both <- cor.test(data_both[[esmFold_col]], data_both[[alphaFold_col]], method = "pearson")
  method_used <- "Pearson's test"
  lm_fit <- lm(data_both[[alphaFold_col]] ~ data_both[[esmFold_col]])  # Fit linear model
  lm_equation <- paste0("y = ", round(coef(lm_fit)[2], 2), "x + ", round(coef(lm_fit)[1], 2))
  regression_line <- geom_smooth(method = "lm", se = FALSE, color = "orange")  # Linear line
} else {
  # At least one distribution is not normal
  correlation_test_both <- cor.test(data_both[[esmFold_col]], data_both[[alphaFold_col]], method = "spearman")
  method_used <- "Spearman's test"
  
  # Fit LOESS model and approximate with polynomial regression
  # LOESS: Locally Estimated Scatterplot Smoothing is a non-parametric regression method 
  loess_fit <- loess(data_both[[alphaFold_col]] ~ data_both[[esmFold_col]], span = 0.5)
  poly_fit <- lm(loess_fit$fitted ~ poly(data_both[[esmFold_col]], 2, raw = TRUE))  # Polynomial fit (degree 2)
  
  # Extract coefficients of the polynomial
  poly_coeff <- coef(poly_fit)
  lm_equation <- paste0("y = ", round(poly_coeff[1], 2), 
                        " + ", round(poly_coeff[2], 2), "x",
                        " + ", round(poly_coeff[3], 4), "x²")
  regression_line <- stat_function(fun = function(x) poly_coeff[1] + poly_coeff[2] * x + poly_coeff[3] * x^2,
                                   color = "orange")  # Add polynomial curve
}


correlation_value_both <- round(correlation_test_both$estimate, 2)
p_value_both <- format.pval(correlation_test_both$p.value, digits = 2)
print(correlation_value_both)
print(p_value_both)

# Define colors
colors <- c("ESMFold" = "#2F79B5", "AlphaFold2" = "#C13639", "Both" = "black")

# Create scatter plot
scatter_plot_pLDDT <- ggplot(data, aes_string(x = esmFold_col, y = alphaFold_col, color = "Method", text = protein_id_col)) +
  geom_point(size = 3.5, alpha = 0.4, stroke = 0.5) +
  scale_color_manual(values = colors) +
  labs(title = "",
       x = "ESMFold",
       y = "AlphaFold2") +
  theme_classic() +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = 50, linetype = "dotted", color = "grey40") +
  regression_line + 
  annotate("text", x = 0, y = 100, label = paste0(
    lm_equation,
    "\nCorr. coeff.: ", correlation_value_both,
    "\nP-value: ", p_value_both), 
    hjust = 0, vjust = 1, size = 5, color = "black")+
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.position = "bottom"
  )

scatter_plot_pLDDT

# Interactive plot
#ggplotly(scatter_plot_pLDDT, tooltip = c("x", "y", "text"))

# Save the plot
ggsave("Zpa796_AF2-ESMFold_scatterplot_pLDDT.pdf", plot = scatter_plot_pLDDT, device = cairo_pdf, width = 8, height = 8)


##################################################################################################
##################################################################################################
#### RMSD analysis

df_melted <- df %>%
  select(Protein.ID, `RMSD.AF2.vs.ESMFold.`) %>%
  pivot_longer(cols = -Protein.ID, names_to = "Method", values_to = "RMSD")

df_melted$RMSD <- as.numeric(as.character(df_melted$RMSD))


##################################################################################################
# Create boxplot with density distribution (violin) for RMSD and tRMSD

plot_mean_RMSD <- ggplot(df_melted, aes(x = Method, y = RMSD, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.4, fill = "grey") + 
  geom_boxplot(width = 0.1, fill = "darkgrey", color = "black") +
  theme_classic() +
  labs(title = "",
       x = "", y = "RMSD") +
  theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_blank(),       
    axis.text.x = element_blank(),        
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    legend.position = "none"
  )
plot_mean_RMSD


##################################################################################################
# Plot RMSDs according to protein length range

df$Length_Category <- cut(df$`Protein.length...AA.`,
                          breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, Inf),
                          labels = c("1-100", "101-200", "201-300", "301-400", "401-500", "501-600", "601-700", "701-800", "801-900", "901-1000", ">1000"))

# Melting the dataframe for plotting by protein length
df_melted_length <- df %>%
  select(Protein.ID, `RMSD.AF2.vs.ESMFold.`, Length_Category) %>%
  pivot_longer(cols = `RMSD.AF2.vs.ESMFold.`, names_to = "Method", values_to = "RMSD") %>%
  mutate(Method = "RMSD") 

df_melted_length$RMSD <- as.numeric(as.character(df_melted_length$RMSD))


# Plotting boxplots for RMSDs by protein length category
plot_RMSD_by_length <- ggplot(df_melted_length, aes(x = Length_Category, y = RMSD)) +
  geom_boxplot(
    fill = "grey", color = "grey50", 
    alpha = 0.4,  
    outlier.shape = NA  
  ) +
  theme_classic() +
  labs(title = "",
       x = "Mature protein length", y = "RMSD") +
  geom_point(
    color = "grey10", 
    position = position_dodge(0.8),
    alpha = 0.4,
    shape = 16  
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_RMSD_by_length 


#########################################################################################
#Statistical tests

# Test for normality within each length category for RMSD
RMSD_normality <- df_melted_length %>%
  group_by(Length_Category) %>%
  shapiro_test(RMSD)
RMSD_normality

# If non-normality is found in any group, use Kruskal-Wallis, otherwise use ANOVA
if(any(RMSD_normality$p < 0.05, na.rm = TRUE)) {
  RMSDs_comparison <- 
    kruskal_test(RMSD ~ Length_Category, data = df_melted_length)
} else {
  RMSDs_comparison <- anova_test(RMSD ~ Length_Category, data = df_melted_length)
}
RMSDs_comparison

# Perform the Dunn test for multiple comparisons
RMSDs_posthoc <- df_melted_length %>%
  dunn_test(RMSD ~ Length_Category, p.adjust.method = "BH")
RMSDs_posthoc


generate_letters <- function(p_values) {
  p_matrix <- matrix(p_values, nrow = length(unique(df_melted_length$Length_Category)), ncol = length(unique(df_melted_length$Length_Category)), byrow = TRUE)
  rownames(p_matrix) <- colnames(p_matrix) <- levels(df_melted_length$Length_Category)
  letters <- multcompLetters(p_matrix, compare = "<=")$Letters
  return(letters)
}

rmsd_p_values <- RMSDs_posthoc$p.adj
RMSD_letters_by_length <- generate_letters(rmsd_p_values)
RMSD_letters_by_length

y_position_rmsd <- max(df_melted_length$RMSD, na.rm = TRUE) + 0.15

# Add letters to RMSD by length plot
for (length_category in levels(df_melted_length$Length_Category)) {
  plot_RMSD_by_length <- plot_RMSD_by_length +
    annotate("text", x = which(levels(df_melted_length$Length_Category) == length_category), y = y_position_rmsd, label = RMSD_letters_by_length[length_category], size = 5, color = "black")
}

plot_mean_RMSD
plot_RMSD_by_length


ggsave("Zpa796_AF2-ESMFold_plot_mean_RMSD.pdf", plot = plot_mean_RMSD, device = cairo_pdf, width = 5, height = 6)
ggsave("Zpa796_AF2-ESMFold_plot_RMSD_by_length.pdf", plot = plot_RMSD_by_length, device = cairo_pdf, width = 8, height = 6)

