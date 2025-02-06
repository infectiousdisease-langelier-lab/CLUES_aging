
library(dplyr)
library(tibble)
library(limma)
library(DESeq2)
library(tidyr)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(GEOquery)
library(data.table)
library(readxl)
library(readr)


## Read in data

ab_list <- c("ENSG00000152778, ENSG00000111331, ENSG00000135114, ENSG00000111335, ENSG00000119922, ENSG00000185745, ENSG00000185507, ENSG00000160710, ENSG00000119917, ENSG00000187608, ENSG00000134321, ENSG00000183486, ENSG00000126709, ENSG00000213928, ENSG00000157601, ENSG00000142089, ENSG00000185885, ENSG00000184979, ENSG00000089127, ENSG00000132530, ENSG00000165949, ENSG00000185338, ENSG00000170581, ENSG00000115415, ENSG00000130303, ENSG00000068079, ENSG00000172183, ENSG00000206503")

g_list <- c("ENSG00000111331, ENSG00000135114, ENSG00000111335, ENSG00000185507, ENSG00000213928, ENSG00000112343, ENSG00000132274, ENSG00000089127, ENSG00000132256, ENSG00000106785, ENSG00000067066, ENSG00000185338, ENSG00000140464, ENSG00000121236, ENSG00000115415, ENSG00000162654, ENSG00000117228, ENSG00000121060, ENSG00000206503, ENSG00000132109, ENSG00000162645, ENSG00000125148, ENSG00000137265, ENSG00000166710, ENSG00000234127, ENSG00000204525, ENSG00000162931, ENSG00000154451, ENSG00000168310, ENSG00000204592, ENSG00000204632, ENSG00000116030, ENSG00000184557, ENSG00000116525, ENSG00000198019, ENSG00000117226, ENSG00000234745")


#DE results from CLUES bulk
all_genes_CLUES <- read_csv()


# IFN ab and g analysis ===================

ab_list <- trimws(unlist(strsplit(ab_list, ",")))

g_list <- trimws(unlist(strsplit(g_list, ",")))

gene_IFNpathway_list <- c(ab_list, g_list)

IFN_SLE <- all_genes_CLUES %>%
  dplyr::filter(ID %in% gene_IFNpathway_list) %>%
  filter(adj.P.Val < 0.05)

IFN_ab_SLE_list <- IFN_SLE %>%
  dplyr::filter(ID %in% ab_list) %>%
  pull(ID)

IFN_g_SLE_list <- IFN_SLE %>%
  dplyr::filter(ID %in% g_list) %>%
  pull(ID)

IFN_ab <- all_genes_CLUES %>%
  dplyr::filter(ID %in% ab_list) %>%
  pull(gene_name)

IFN_g <- all_genes_CLUES %>%
  dplyr::filter(ID %in% g_list) %>%
  pull(gene_name)


# Read cleaned genecounts data -----
SLE_counts <- read.csv(
  "cleaned_CLUES_bulkcounts.csv",
  row.names = 1,
  check.names = FALSE
)


SLE_norm <- edgeR::cpm(SLE_counts, log=TRUE)

# Read cleaned metadata -----

SLE_cleaned_metadata <- read.csv(
  "dxage_meta.csv", row.names = 1, check.names = FALSE
)

# SLE analysis ----

SLE_score <- data.frame(
  sampleid = colnames(SLE_norm),
  score_ab = colMeans(SLE_norm[IFN_ab_SLE_list, , drop = FALSE]),
  score_g = colMeans(SLE_norm[IFN_g_SLE_list, , drop = FALSE])
)

SLE_cleaned_metadata$sampleid <- SLE_cleaned_metadata$SampleID

SLE_score <- SLE_score %>%
  inner_join(., SLE_cleaned_metadata, by = "sampleid")


# Regression
# Fit linear models for score_ab and score_g

SLE_score$raceeth <- as.factor(SLE_score$RaceEth)
SLE_score$female <- as.factor(SLE_score$female)
SLE_score$AGEDX_factor <- as.factor(SLE_score$AGEDX_factor)


lm_ab <- lm(score_ab ~ Age + female + RaceEth + AGEDX_factor, data = SLE_score)
lm_g <- lm(score_g ~ Age + female + RaceEth + AGEDX_factor, data = SLE_score)

# Extract p-values
age_coef_ab <- summary(lm_ab)$coefficients["Age", "Estimate"]
p_value_ab <- summary(lm_ab)$coefficients["Age", "Pr(>|t|)"]

age_coef_g <- summary(lm_g)$coefficients["Age", "Estimate"]
p_value_g <- summary(lm_g)$coefficients["Age", "Pr(>|t|)"]

# Plot

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )


g_SLE_plot <- ggplot(SLE_score, aes(x = Age)) +
  geom_point(aes(y = score_g)) +    # Points for score_g
  labs(
    title = "SLE, adjusting for diagnosis age",
    x = "Age",
    y = "IFN-γ Score"
  ) +
  geom_smooth(aes(y = score_g, color = "#5e3c99"), method = "lm", se = TRUE, fill = "#5e3c99", color = "#5e3c99") +
  annotate("text", x = max(SLE_score$Age), y = max(SLE_score$score_g),
           label = paste("Coef:", signif(age_coef_g, 3), 
                         "\nP:", signif(p_value_g, 3)), 
           hjust = 1) +
  my.theme

ab_SLE_plot <- ggplot(SLE_score, aes(x = Age)) +
  geom_point(aes(y = score_ab)) +    # Points for score_g
  labs(
    title = "SLE, adjusting for diagnosis age",
    x = "Age",
    y = "IFN-α/β Score"
  ) +
  geom_smooth(aes(y = score_g, color = "#5e3c99"), method = "lm", se = TRUE, fill = "#5e3c99", color = "#5e3c99") +
  annotate("text", x = max(SLE_score$Age), y = max(SLE_score$score_ab),
           label = paste("Coef:", signif(age_coef_ab, 3), 
                         "\nP:", signif(p_value_ab, 3)), 
           hjust = 1) +
  my.theme



# Boxplot

df_long <- data.frame(
  sampleid = rep(SLE_score$sampleid, 2),  # Repeat sample IDs
  AGEDX_factor = rep(SLE_score$AGEDX_factor, 2),  # Repeat AGEDX factors
  Score_Type = rep(c("IFN-α/β", "IFN-γ"), each = nrow(SLE_score)),  # Score type
  Score = c(SLE_score$score_ab, SLE_score$score_g))  # Combine scores


ggplot(df_long, aes(x = Score_Type, y = Score, fill = AGEDX_factor)) +
  geom_boxplot() +
  labs(
    title = "",
    x = "Pathway",
    y = "Score",
    fill = "Diagnosis Age"
  ) +
  my.theme
