library(ggplot2)
library(extrafont)

font_import() 
loadfonts(device = "pdf") 

secretome_data <- read.csv("Zpa796_secretome_metadata.tsv", sep="\t")
#str(secretome_data)
colnames(secretome_data)


categorize_length <- function(length) {
  if (length <= 100) {
    return('1-100')
  } else if (length <= 200) {
    return('101-200')
  } else if (length <= 300) {
    return('201-300')
  } else if (length <= 400) {
    return('301-400')
  } else if (length <= 500) {
    return('401-500')
  } else if (length <= 600) {
    return('501-600')
  } else if (length <= 700) {
    return('601-700')
  } else if (length <= 800) {
    return('701-800')
  } else if (length <= 900) {
    return('801-900')
  } else if (length <= 1000) {
    return('901-1000')
  } else {
    return('>1000')
  }
}

#############################################
######## AlphaFold2-based prediction ########

secretome_data$Protein_Size_Group <- sapply(secretome_data$"Mature.protein...AA.", categorize_length)

secretome_data$pLDDT <- cut(secretome_data$pLDDT..AF2., 
                            breaks = c(-Inf, 50, 70, 90, Inf),
                            labels = c("<50", "50-70", "70-90", ">90"))


secretome_data$pLDDT <- as.character(secretome_data$pLDDT)
secretome_data$pLDDT[is.na(secretome_data$pLDDT..AF2.)] <- "NA"


secretome_data$pLDDT <- factor(secretome_data$pLDDT, levels = c("<50", "50-70", "70-90", ">90", "NA"))

secretome_data$Protein_Size_Group <- factor(secretome_data$Protein_Size_Group, 
                                            levels = c('1-100', '101-200', '201-300', '301-400', '401-500', 
                                                       '501-600', '601-700', '701-800', '801-900', '901-1000', '>1000'))

# Generate the distribution plot for pLDDT (AlphaFold2-based) 
protein_size_barplot <- ggplot(secretome_data, aes(x = Protein_Size_Group, fill = pLDDT)) +
  geom_bar(width = 0.7) +
  labs(title = "",
    x = "Mature protein legnth", y = "Number of proteins") +
  scale_fill_manual(values = c("<50" = "darkorange", "50-70" = "gold", "70-90" = "dodgerblue", ">90" = "darkblue", "NA" = "grey")) +
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
protein_size_barplot

ggsave("Zpa796_pLDDT_AF2_mature_secretome_protein_length_distribution.pdf", plot = protein_size_barplot, device = cairo_pdf, width = 8, height = 8)

###########################################################
# Generate the distribution plot for antimicrobial activity

secretome_data$AM.prediction..AF2.AMAPEC. <- factor(secretome_data$AM.prediction..AF2.AMAPEC., levels = c("Non-antimicrobial", "Antimicrobial", "NA"))

secretome_data$AM.prediction..AF2.AMAPEC.[is.na(secretome_data$AM.prediction..AF2.AMAPEC.)] <- "NA"

secretome_data$AM.prediction..AF2.AMAPEC. <- factor(secretome_data$AM.prediction..AF2.AMAPEC., 
                                                    levels = c("Non-antimicrobial", "Antimicrobial", "NA"))

# Generate the antimicrobial barplot (AlphaFold2-based) 
antimicrobial_barplot <- ggplot(secretome_data, aes(x = Protein_Size_Group, fill = AM.prediction..AF2.AMAPEC.)) +
  geom_bar(width = 0.7) +
  labs(title = "",
    x = "Mature protein legnth", y = "Number of proteins") +
  scale_fill_manual(values = c("Antimicrobial" = "#2171B5", "Non-antimicrobial" = "#41AB5D", "NA" = "grey")) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )+
  labs(fill = "Antimicrobial activity")

antimicrobial_barplot

ggsave("Zpa796_AM-activity_AF2_mature_secretome_protein_length_distribution.pdf", plot = antimicrobial_barplot,  device = cairo_pdf, width = 8, height = 8)


#############################################
######## ESMFold-based prediction ########

secretome_data$Protein_Size_Group <- sapply(secretome_data$"Mature.protein...AA.", categorize_length)

secretome_data$pLDDT <- cut(secretome_data$pLDDT..ESMFold., 
                            breaks = c(-Inf, 50, 70, 90, Inf),
                            labels = c("<50", "50-70", "70-90", ">90"))


secretome_data$pLDDT <- as.character(secretome_data$pLDDT)
secretome_data$pLDDT[is.na(secretome_data$pLDDT..ESMFold.)] <- "NA"


secretome_data$pLDDT <- factor(secretome_data$pLDDT, levels = c("<50", "50-70", "70-90", ">90", "NA"))

secretome_data$Protein_Size_Group <- factor(secretome_data$Protein_Size_Group, 
                                            levels = c('1-100', '101-200', '201-300', '301-400', '401-500', 
                                                       '501-600', '601-700', '701-800', '801-900', '901-1000', '>1000'))

# Generate the distribution plot for pLDDT (ESMFold-based) 
protein_size_barplot <- ggplot(secretome_data, aes(x = Protein_Size_Group, fill = pLDDT)) +
  geom_bar(width = 0.7) +
  labs(title = "",
    x = "Mature protein legnth", y = "Number of proteins") +
  scale_fill_manual(values = c("<50" = "darkorange", "50-70" = "gold", "70-90" = "dodgerblue", ">90" = "darkblue", "NA" = "grey")) +
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

protein_size_barplot

ggsave("Zpa796_pLDDT_ESMFold_mature_secretome_protein_length_distribution.pdf", plot = protein_size_barplot, device = cairo_pdf, width = 8, height = 8)


###########################################################
###########################################################
### Correlation between pLDDT values and mature protein length 

library(dplyr)
library(ggplot2)
library(tidyr)
library(mgcv)


#############################################
######## AlphaFold2-based prediction ########


df_filtered_AF2 <- secretome_data %>%
  select(Protein.ID, pLDDT_AF2 = `pLDDT..AF2.`, Mature_length_AF2 = `Mature.protein...AA.`) %>%
  filter(!is.na(pLDDT_AF2) & !is.na(Mature_length_AF2))

df_filtered_AF2 <- df_filtered_AF2 %>%
  mutate(pLDDT_range_AF2 = case_when(
    pLDDT_AF2 < 50 ~ "<50",
    pLDDT_AF2 >= 50 & pLDDT_AF2 < 70 ~ "50-70",
    pLDDT_AF2 >= 70 & pLDDT_AF2 < 90 ~ "70-90",
    pLDDT_AF2 >= 90 ~ ">90"
  ))

df_filtered_AF2$pLDDT_range_AF2 <- factor(df_filtered_AF2$pLDDT_range_AF2, levels = c("<50", "50-70", "70-90", ">90"))


linear_model_AF2 <- lm(pLDDT_AF2 ~ Mature_length_AF2, data = df_filtered_AF2)
poly_model_AF2 <- lm(pLDDT_AF2 ~ poly(Mature_length_AF2, 2), data = df_filtered_AF2)
gam_model_AF2 <- gam(pLDDT_AF2 ~ s(Mature_length_AF2), data = df_filtered_AF2)

# Calculate BIC for each model
bic_values_AF2 <- c(
  linear = BIC(linear_model_AF2),
  polynomial = BIC(poly_model_AF2),
  gam = BIC(gam_model_AF2)
)

# Choose the best model
best_model_name_AF2 <- names(which.min(bic_values_AF2))
print(bic_values_AF2)
print(paste("Best model (AF2):", best_model_name_AF2))

# Create scatter plot
AF2_pLDDT_scatterplot <- ggplot(df_filtered_AF2, aes(x = Mature_length_AF2, y = pLDDT_AF2, color = pLDDT_range_AF2)) +
  geom_point(size = 3.5, alpha = 0.4, stroke = 0.5) +
  scale_color_manual(values = c(
    "<50" = "darkorange",
    "50-70" = "gold",
    "70-90" = "dodgerblue",
    ">90" = "darkblue"
  )) +
  labs(title = "",
       x = "Mature protein length",
       y = "pLDDT",
       color = "pLDDT") +
  ylim(20, 100) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.position = "bottom"
  )
AF2_pLDDT_scatterplot


if (best_model_name_AF2 == "linear") {
  AF2_pLDDT_scatterplot <- AF2_pLDDT_scatterplot +
    geom_smooth(method = "lm", color = "grey40", se = FALSE)
} else if (best_model_name_AF2 == "polynomial") {
  AF2_pLDDT_scatterplot <- AF2_pLDDT_scatterplot +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "grey40", se = FALSE)
} else if (best_model_name_AF2 == "gam") {
  AF2_pLDDT_scatterplot <- AF2_pLDDT_scatterplot +
    geom_smooth(method = "gam", formula = y ~ s(x), color = "grey40", se = FALSE)
}
AF2_pLDDT_scatterplot

# Check the distribution of pLDDT_AF2
ggplot(df_filtered_AF2, aes(x = pLDDT_AF2)) +
  geom_histogram(binwidth = 5, fill = "blue", alpha = 0.7) +
  labs(title = "Distribution of pLDDT values", x = "pLDDT", y = "Frequency")

# Check the distribution of Mature_length
ggplot(df_filtered_AF2, aes(x = Mature_length_AF2)) +
  geom_histogram(binwidth = 10, fill = "green", alpha = 0.7) +
  labs(title = "Distribution of Mature Protein Lengths", x = "Mature Length (#AA)", y = "Frequency")

# Shapiro-Wilk normality test for pLDDT_AF2
shapiro_pLDDT_AF2 <- shapiro.test(df_filtered_AF2$pLDDT_AF2)
shapiro_pLDDT_AF2

# Shapiro-Wilk normality test for Mature_length
shapiro_Mature_length_AF2 <- shapiro.test(df_filtered_AF2$Mature_length_AF2)
shapiro_Mature_length_AF2

# Perform correlation test based on normality results
if (shapiro_pLDDT_AF2$p.value > 0.05 & shapiro_Mature_length_AF2$p.value > 0.05) {
  # Both distributions are normal
  correlation_test_AF2 <- cor.test(df_filtered_AF2$pLDDT_AF2, df_filtered_AF2$Mature_length_AF2, method = "pearson")
  method_used_AF2<- "Pearson's test"
} else {
  # At least one distribution is not normal
  correlation_test_AF2 <- cor.test(df_filtered_AF2$pLDDT_AF2, df_filtered_AF2$Mature_length_AF2, method = "spearman")
  method_used_AF2 <- "Spearman's test"
}

correlation_value_AF2 <- round(correlation_test_AF2$estimate, 2)
p_value_AF2 <- format.pval(correlation_test_AF2$p.value, digits = 2)

# Create the final scatter plot with annotation and best fit curve
AF2_pLDDT_scatterplot <- AF2_pLDDT_scatterplot + 
  annotate("text", x = Inf, y = Inf, label = paste("Corr. coeff.:", correlation_value_AF2, "\nP-value:", p_value_AF2),
           hjust = 1.5, vjust = 1.5, size = 6, color = "black")

AF2_pLDDT_scatterplot

ggsave("Zpa796_pLDDT_AF2_mature_secretome_protein_length_scatterplot.pdf", plot =AF2_pLDDT_scatterplot, device = cairo_pdf, width = 8, height = 8)


#############################################
######## ESMFold-based prediction ########


df_filtered_ESMFold <- secretome_data %>%
  select(Protein.ID, pLDDT_ESMFold = `pLDDT..ESMFold.`, Mature_length_ESMFold = `Mature.protein...AA.`) %>%
  filter(!is.na(pLDDT_ESMFold) & !is.na(Mature_length_ESMFold))

df_filtered_ESMFold <- df_filtered_ESMFold %>%
  mutate(pLDDT_range_ESMFold = case_when(
    pLDDT_ESMFold < 50 ~ "<50",
    pLDDT_ESMFold >= 50 & pLDDT_ESMFold < 70 ~ "50-70",
    pLDDT_ESMFold >= 70 & pLDDT_ESMFold < 90 ~ "70-90",
    pLDDT_ESMFold >= 90 ~ ">90"
  ))

df_filtered_ESMFold$pLDDT_range_ESMFold <- factor(df_filtered_ESMFold$pLDDT_range_ESMFold, levels = c("<50", "50-70", "70-90", ">90"))


linear_model_ESMFold <- lm(pLDDT_ESMFold ~ Mature_length_ESMFold, data = df_filtered_ESMFold)
poly_model_ESMFold <- lm(pLDDT_ESMFold ~ poly(Mature_length_ESMFold, 2), data = df_filtered_ESMFold)
gam_model_ESMFold <- gam(pLDDT_ESMFold ~ s(Mature_length_ESMFold), data = df_filtered_ESMFold)

# Calculate BIC for each model
bic_values_ESMFold <- c(
  linear = BIC(linear_model_ESMFold),
  polynomial = BIC(poly_model_ESMFold),
  gam = BIC(gam_model_ESMFold)
)

# Choose the best model
best_model_name_ESMFold <- names(which.min(bic_values_ESMFold))
print(bic_values_ESMFold)
print(paste("Best model (ESMFold):", best_model_name_ESMFold))

# Create scatter plot
ESMFold_pLDDT_scatterplot <- ggplot(df_filtered_ESMFold, aes(x = Mature_length_ESMFold, y = pLDDT_ESMFold, color = pLDDT_range_ESMFold)) +
  geom_point(size = 3.5, alpha = 0.4, stroke = 0.5) +
  scale_color_manual(values = c(
    "<50" = "darkorange",
    "50-70" = "gold",
    "70-90" = "dodgerblue",
    ">90" = "darkblue"
  )) +
  labs(title = "",
       x = "Mature protein length",
       y = "pLDDT",
       color = "pLDDT") +
  ylim(20, 100) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.position = "bottom"
  )
ESMFold_pLDDT_scatterplot


if (best_model_name_ESMFold == "linear") {
  ESMFold_pLDDT_scatterplot <- ESMFold_pLDDT_scatterplot +
    geom_smooth(method = "lm", color = "grey40", se = FALSE)
} else if (best_model_name_ESMFold == "polynomial") {
  ESMFold_pLDDT_scatterplot <- ESMFold_pLDDT_scatterplot +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "grey40", se = FALSE)
} else if (best_model_name_ESMFold == "gam") {
  ESMFold_pLDDT_scatterplot <- ESMFold_pLDDT_scatterplot +
    geom_smooth(method = "gam", formula = y ~ s(x), color = "grey40", se = FALSE)
}
ESMFold_pLDDT_scatterplot

# Check the distribution of pLDDT_ESMFold
ggplot(df_filtered_ESMFold, aes(x = pLDDT_ESMFold)) +
  geom_histogram(binwidth = 5, fill = "blue", alpha = 0.7) +
  labs(title = "Distribution of pLDDT values", x = "pLDDT", y = "Frequency")

# Check the distribution of Mature_length
ggplot(df_filtered_ESMFold, aes(x = Mature_length_ESMFold)) +
  geom_histogram(binwidth = 10, fill = "green", alpha = 0.7) +
  labs(title = "Distribution of Mature Protein Lengths", x = "Mature Length (#AA)", y = "Frequency")

# Shapiro-Wilk normality test for pLDDT_ESMFold
shapiro_pLDDT_ESMFold <- shapiro.test(df_filtered_ESMFold$pLDDT_ESMFold)
shapiro_pLDDT_ESMFold

# Shapiro-Wilk normality test for Mature_length
shapiro_Mature_length_ESMFold <- shapiro.test(df_filtered_ESMFold$Mature_length_ESMFold)
shapiro_Mature_length_ESMFold

# Perform correlation test based on normality results
if (shapiro_pLDDT_ESMFold$p.value > 0.05 & shapiro_Mature_length_ESMFold$p.value > 0.05) {
  # Both distributions are normal
  correlation_test_ESMFold <- cor.test(df_filtered_ESMFold$pLDDT_ESMFold, df_filtered_ESMFold$Mature_length_ESMFold, method = "pearson")
  method_used_ESMFold<- "Pearson's test"
} else {
  # At least one distribution is not normal
  correlation_test_ESMFold <- cor.test(df_filtered_ESMFold$pLDDT_ESMFold, df_filtered_ESMFold$Mature_length_ESMFold, method = "spearman")
  method_used_ESMFold <- "Spearman's test"
}

correlation_value_ESMFold <- round(correlation_test_ESMFold$estimate, 2)
p_value_ESMFold <- format.pval(correlation_test_ESMFold$p.value, digits = 2)

# Create the final scatter plot with annotation and best fit curve
ESMFold_pLDDT_scatterplot <- ESMFold_pLDDT_scatterplot + 
  annotate("text", x = Inf, y = Inf, label = paste("Corr. coeff.:", correlation_value_ESMFold, "\nP-value:", p_value_ESMFold),
           hjust = 1.5, vjust = 1.5, size = 6, color = "black")

ESMFold_pLDDT_scatterplot

ggsave("Zpa796_pLDDT_ESMFold_mature_secretome_protein_length_scatterplot.pdf", plot =ESMFold_pLDDT_scatterplot,  device = cairo_pdf, width = 8, height = 8)


