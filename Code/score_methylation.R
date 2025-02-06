
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


# DE results of bulk analysis

all_genes_RS <- read_csv()
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

IFN_RS_ab <- all_genes_RS %>%
  dplyr::select(-c(gene_name...9)) %>%
  dplyr::filter(gene_name...1 %in% IFN_ab) %>%
  filter(adj.P.Val < 0.05) %>% pull(gene_name...1)

IFN_RS_g <- all_genes_RS %>%
  dplyr::select(-c(gene_name...9)) %>%
  dplyr::filter(gene_name...1 %in% IFN_g) %>%
  filter(adj.P.Val < 0.05) %>% pull(gene_name...1)


# Read cleaned genecounts data -----
SLE_counts <- read.csv(
  "./cleaned_CLUES_bulkcounts.csv",
  row.names = 1,
  check.names = FALSE
)


SLE_norm <- edgeR::cpm(SLE_counts, log=TRUE)

# Read cleaned metadata -----

SLE_cleaned_metadata <- read.csv(
  "./cleaned_CLUES_metadata.csv", row.names = 1, check.names = FALSE
)

# SLE analysis ----

SLE_score <- data.frame(
  sampleid = colnames(SLE_norm),
  score_ab = colMeans(SLE_norm[IFN_ab_SLE_list, , drop = FALSE]),
  score_g = colMeans(SLE_norm[IFN_g_SLE_list, , drop = FALSE])
)

SLE_score <- SLE_score %>%
  inner_join(., SLE_cleaned_metadata, by = "sampleid")


# Regression
# Fit linear models for score_ab and score_g

SLE_score$raceeth <- as.factor(SLE_score$raceeth)
SLE_score$female <- as.factor(SLE_score$female)


lm_ab <- lm(score_ab ~ age + female + raceeth, data = SLE_score)
lm_g <- lm(score_g ~ age + female + raceeth, data = SLE_score)

# Extract p-values
age_coef_ab <- summary(lm_ab)$coefficients["age", "Estimate"]
p_value_ab <- summary(lm_ab)$coefficients["age", "Pr(>|t|)"]

age_coef_g <- summary(lm_g)$coefficients["age", "Estimate"]
p_value_g <- summary(lm_g)$coefficients["age", "Pr(>|t|)"]

# Plot

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )


ggplot(SLE_score, aes(x = age)) +
  geom_point(aes(y = score_ab, color = "IFNab")) +  # Points for score_ab
  geom_point(aes(y = score_g, color = "IFNg")) +    # Points for score_g
  # Trend lines
  geom_smooth(aes(y = score_ab, color = "IFNab"), method = "lm", se = FALSE) +
  geom_smooth(aes(y = score_g, color = "IFNg"), method = "lm", se = FALSE) +
  annotate("text", x = max(SLE_score$age), y = max(SLE_score$score_ab),
           label = paste("Coef (AB):", signif(age_coef_ab, 3), 
                         "\nP-value (AB):", signif(p_value_ab, 3)), 
           hjust = 1) +
  annotate("text", x = max(SLE_score$age), y = max(SLE_score$score_g),
           label = paste("Coef (G):", signif(age_coef_g, 3), 
                         "\nP-value (G):", signif(p_value_g, 3)), 
           hjust = 1) +
  labs(
    title = "Age vs. Interferon Scores (SLE)",
    x = "Age",
    y = "Score (average expression)",
    color = "Pathway"
  ) +
  my.theme


g_SLE_plot <- ggplot(SLE_score, aes(x = age)) +
  geom_point(aes(y = score_g)) +    # Points for score_g
  labs(
    title = "SLE",
    x = "Age",
    y = "IFNg Score"
  ) +
  geom_smooth(aes(y = score_g, color = "#5e3c99"), method = "lm", se = TRUE, fill = "#5e3c99", color = "#5e3c99") +
  annotate("text", x = max(SLE_score$age), y = max(SLE_score$score_g),
           label = paste("Coef:", signif(age_coef_g, 3), 
                         "\nP:", signif(p_value_g, 3)), 
           hjust = 1) +
  my.theme

ab_SLE_plot <- ggplot(SLE_score, aes(x = age)) +
  geom_point(aes(y = score_ab)) +    # Points for score_g
  labs(
    title = "SLE",
    x = "Age",
    y = "IFNab Score"
  ) +
  geom_smooth(aes(y = score_g, color = "#5e3c99"), method = "lm", se = TRUE, fill = "#5e3c99", color = "#5e3c99") +
  annotate("text", x = max(SLE_score$age), y = max(SLE_score$score_ab),
           label = paste("Coef:", signif(age_coef_ab, 3), 
                         "\nP:", signif(p_value_ab, 3)), 
           hjust = 1) +
  my.theme

ggsave("/Users/rithwikn/Documents/CZBiohub_2024/revisions/Reviewer2/ifn_score/SLE_ab.svg", 
       ab_SLE_plot,
       width = 4,
       height = 4)

ggsave("/Users/rithwikn/Documents/CZBiohub_2024/revisions/Reviewer2/ifn_score/SLE_g.svg", 
       g_SLE_plot,
       width = 4,
       height = 4)


# RS analysis ====

# Read metadata ====

RS <- read.csv(
  "./cleaned_RS_metadata.csv"
)

RS <- as.data.frame(RS)

rownames(RS) <- RS$RS3.ID

# Read cleaned microarray data -----

expression = fread("./cleaned_RS_microarray.csv", header = TRUE)


## Map Probe IDs to gene symbol

gse <- getGEO(GEO = "GSE33828", GSEMatrix = TRUE)

## Fetch feature data to complete mapping

feature.data <- gse$GSE33828_series_matrix.txt.gz@featureData@data
feature.data <- feature.data[,c(1,6)]


genex <- as.data.frame(expression)

nams = genex$V1
rownames(genex) <- make.names(nams, unique = TRUE)
genex <- genex[,-c(1)]

genex$V1 <- rownames(genex)

genex <- genex %>% 
  inner_join(., feature.data, by = c('V1' = 'ID'))


genex <- dplyr::select(genex, -c("V1"))

rownames(genex) <- make.names(genex$ILMN_Gene, unique = TRUE)


genex <- dplyr::select(genex, -c("ILMN_Gene"))
qsave(genex, "./RS_cleaned_counts.qs")

# RS score

RS_score <- data.frame(
  RS3.ID = as.integer(colnames(genex)),
  score_ab = colMeans(genex[IFN_RS_ab, , drop = FALSE]),
  score_g = colMeans(genex[IFN_RS_g, , drop = FALSE])
)


RS_score <- RS_score %>%
  inner_join(., RS, by = "RS3.ID")



# Regression =====

RS$sex <- as.factor(RS$sex)


lm_ab_RS <- lm(score_ab ~ age + sex, data = RS_score)
lm_g_RS <- lm(score_g ~ age + sex, data = RS_score)

# Extract p-values
age_coef_ab_RS <- summary(lm_ab_RS)$coefficients["age", "Estimate"]
p_value_ab_RS <- summary(lm_ab_RS)$coefficients["age", "Pr(>|t|)"]

age_coef_g_RS <- summary(lm_g_RS)$coefficients["age", "Estimate"]
p_value_g_RS <- summary(lm_g_RS)$coefficients["age", "Pr(>|t|)"]

# plot ====

ggplot(RS_score, aes(x = age)) +
  geom_point(aes(y = score_ab, color = "IFNab")) +  # Points for score_ab
  geom_point(aes(y = score_g, color = "IFNg")) +    # Points for score_g
  labs(
    title = "Age vs. Interferon Scores (controls)",
    x = "Age",
    y = "Score (average expression)",
    color = ""
  ) +
  geom_smooth(aes(y = score_ab, color = "IFNab"), method = "lm", se = FALSE) +
  geom_smooth(aes(y = score_g, color = "IFNg"), method = "lm", se = FALSE) +
  annotate("text", x = max(RS_score$age), y = max(RS_score$score_ab),
           label = paste("Coef (AB):", signif(age_coef_ab_RS, 3), 
                         "\nP-value (AB):", signif(p_value_ab_RS, 3)), 
           hjust = 1) +
  annotate("text", x = max(RS_score$age), y = max(RS_score$score_g),
           label = paste("Coef (G):", signif(age_coef_g_RS, 3), 
                         "\nP-value (G):", signif(p_value_g_RS, 3)), 
           hjust = 1) +
  my.theme



g_rs_plot <- ggplot(RS_score, aes(x = age)) +
  geom_point(aes(y = score_g)) +    # Points for score_g
  labs(
    title = "Controls",
    x = "Age",
    y = "IFNg Score"
  ) +
  geom_smooth(aes(y = score_g), method = "lm", se = TRUE, fill = "#ca0020", color = "#ca0020") +
  annotate("text", x = max(RS_score$age), y = max(RS_score$score_g),
           label = paste("Coef:", signif(age_coef_g_RS, 3), 
                         "\nP:", signif(p_value_g_RS, 3)), 
           hjust = 1) +
  my.theme

ab_rs_plot <- ggplot(RS_score, aes(x = age)) +
  geom_point(aes(y = score_ab)) +    # Points for score_g
  labs(
    title = "Controls",
    x = "Age",
    y = "IFNab Score"
  ) +
  geom_smooth(aes(y = score_ab), method = "lm", se = TRUE, fill = "#ca0020", color = "#ca0020") +
  annotate("text", x = max(RS_score$age), y = max(RS_score$score_ab),
           label = paste("Coef:", signif(age_coef_ab_RS, 3), 
                         "\nP:", signif(p_value_ab_RS, 3)), 
           hjust = 1) +
  my.theme


ggsave("/Users/rithwikn/Documents/CZBiohub_2024/revisions/Reviewer2/ifn_score/SLE_g.svg", 
       g_rs_plot,
       width = 4,
       height = 4)

ggsave("/Users/rithwikn/Documents/CZBiohub_2024/revisions/Reviewer2/ifn_score/SLE_ab.svg", 
       ab_rs_plot,
       width = 4,
       height = 4)


# Methlyation analysis SLE patients ========

data_dir <- "~/Documents/CZBiohub_2024/revisions/methylation_analysis"
interferome <- read.csv(paste0(data_dir, "/interferome.csv"))
Cpgs_overlap2 <- read.csv(paste0(data_dir, "/significant_Cpgs_overlap.csv"))
M_fltered <- readRDS(paste0(data_dir, "/Mvalues_filtered.RDS"))
meta_clean <- read.csv(paste0(data_dir, "/metadata_clean.csv"))
down_cpgs <- read.csv(paste0(data_dir, "/HyperCpgs_DownRNA.csv"))



SLE_norm_interferome <- SLE_norm %>%
  as.data.frame(.) %>%
  rownames_to_column("ensg") %>%
  dplyr::filter(ensg %in% interferome$Ensembl.Id) %>%
  column_to_rownames("ensg")

cpg_mapping <- Cpgs_overlap %>%
  dplyr::select(gene) %>%
  inner_join(., all_genes_CLUES, by = c("gene" = "gene_name")) %>%
  dplyr::select(gene, ID)

Cpgs_overlap <- Cpgs_overlap %>%
  dplyr::filter(gene %in% cpg_mapping$gene)

Cpgs_overlap <- cbind(Cpgs_overlap, cpg_mapping)

genex_genes <- rownames(SLE_norm_interferome)

# Fix non-unique column names
colnames(Cpgs_overlap) <- make.names(colnames(Cpgs_overlap), unique = TRUE)

Cpgs_overlap <- Cpgs_overlap %>%
  dplyr::filter(ID %in% genex_genes)

## read in and edit methylation data
selected_v <- as.data.frame(t(M_fltered[Cpgs_overlap$X,]))

sampleids <- as.data.frame(colnames(SLE_counts))
sampleids <- sampleids %>%
  dplyr::mutate(mapped = sub("([^D]+)D.*", "\\1", sampleids$`colnames(SLE_counts)`)) 

colnames(sampleids)[1] <- "count_names"

M_metadata$SubID <- as.character(M_metadata$SubID)

M_metadata <- M_metadata %>%
  inner_join(., sampleids, by = c("SubID" = "mapped"))


# Calculate M values
Cpgs_overlap_filtered <- Cpgs_overlap[Cpgs_overlap$X %in% colnames(selected_v), ]
cpg_to_gene <- setNames(Cpgs_overlap_filtered$ID, Cpgs_overlap_filtered$X)


gene_averages <- split(colnames(selected_v), cpg_to_gene)
gene_matrix <- sapply(names(gene_averages), function(gene) {
  cpgs_for_gene <- gene_averages[[gene]] 
  rowMeans(selected_v[, cpgs_for_gene, drop = FALSE], na.rm = TRUE)
})

final_gene_matrix <- as.data.frame(gene_matrix)
rownames(final_gene_matrix) <- rownames(selected_v)

# Rename columns

m_ids_mapped <- M_metadata %>%
  dplyr::select(c(sampleID_r, count_names))

final_gene_matrix <- final_gene_matrix %>%
  rownames_to_column("sampleID_r") %>%
  inner_join(., m_ids_mapped, by = "sampleID_r") %>%
  column_to_rownames("count_names") %>%
  dplyr::select(-c(sampleID_r))

final_gene_matrix_tranposed <- as.data.frame(t(final_gene_matrix))
final_SLE_norm_data <- SLE_norm_interferome[rownames(final_gene_matrix_tranposed), ]


common_cols <- intersect(colnames(final_SLE_norm_data), colnames(final_gene_matrix_tranposed)) 


final_SLE_norm_data <- final_SLE_norm_data[, common_cols]
final_gene_matrix_tranposed <- final_gene_matrix_tranposed[, common_cols]


correlation_results <- data.frame(Gene = rownames(final_SLE_norm_data), 
                                  Correlation = NA, 
                                  P_value = NA)

# Loop through each gene
for (gene in rownames(final_SLE_norm_data)) {
  
  # Extract expression and methylation values for the gene
  expr_values <- as.numeric(final_SLE_norm_data[gene, ])
  meth_values <- as.numeric(final_gene_matrix_tranposed[gene, ])
  
  # Perform correlation test (returns both correlation coefficient and p-value)
  cor_test <- cor.test(meth_values, expr_values, method = "pearson")
  
  # Store correlation coefficient and p-value in the dataframe
  correlation_results[correlation_results$Gene == gene, "Correlation"] <- cor_test$estimate
  correlation_results[correlation_results$Gene == gene, "P_value"] <- cor_test$p.value
}

correlation_results$Adj_P_Value <- p.adjust(correlation_results$P_value, method = "BH")

correlation_results$NegLogPadj <- -log10(correlation_results$Adj_P_Value)

ggplot(correlation_results, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.1, fill = "#5e3c99", color = "black", alpha = 0.7) +
  labs(
    title = "Distribution of correlation coefficients",
    x = "Correlation coefficient",
    y = "Count"
  ) + 
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     limits = c(-1, 1)) + 
  my.theme


ensg_mapping <- all_genes_CLUES %>%
  dplyr::select(c(ID, gene_name))

correlation_results <- correlation_results %>%
  inner_join(., ensg_mapping, by = c("Gene" = "ID")) 

correlation_results <- correlation_results %>%
  mutate(label = ifelse(NegLogPadj > 20, gene_name, NA))

volcano <- ggplot(correlation_results, aes(x = Correlation, y = NegLogPadj)) +
  geom_point(aes(color = NegLogPadj > 1.3), size = 2) +  
  scale_color_manual(values = c("gray", "#5e3c99")) +   
  geom_text_repel(aes(label = label), 
                  max.overlaps = Inf,  
                  box.padding = 0.8,  
                  point.padding = 0.5, 
                  segment.color = "grey50") +
  labs(x = "Pearson correlation", y = expression(-log[10](P[adj]))) +
  my.theme + 
  theme(legend.position = "none")  +
  scale_x_continuous(limits = c(-1, 1)) 

ggsave(paste0(data_dir, "/volcano_M.svg"), volcano, width = 5, height = 5)


## Same analysis for downregulated hypermethylated cpgs

cpg_mapping <- cpg_mapping %>% distinct(gene, .keep_all = TRUE)

down_cpgs <- down_cpgs %>% merge(., cpg_mapping, by = "gene")

down_cpgs_genes <- down_cpgs %>%
  dplyr::distinct(ID) %>%
  dplyr::pull(ID)

SLE_norm_down_cpgs <- SLE_norm %>%
  as.data.frame(.) %>%
  rownames_to_column("ensg") %>%
  dplyr::filter(ensg %in% down_cpgs_genes) %>%
  column_to_rownames("ensg")


selected_cps <- as.data.frame(t(M_fltered[down_cpgs$X,]))
# Calculate M values
down_hyper_Cpgs <- down_cpgs[down_cpgs$X %in% colnames(selected_cps), ]
dh_cpg_to_gene <- setNames(down_hyper_Cpgs$ID, down_hyper_Cpgs$X)


dh_gene_averages <- split(colnames(selected_cps), dh_cpg_to_gene)
dh_gene_matrix <- sapply(names(dh_gene_averages), function(gene) {
  cpgs_for_gene <- dh_gene_averages[[gene]] 
  rowMeans(selected_cps[, cpgs_for_gene, drop = FALSE], na.rm = TRUE)
})


dh_final_gene_matrix <- as.data.frame(dh_gene_matrix)
#rownames(dh_final_gene_matrix) <- rownames(selected_cps)



dh_final_gene_matrix <- dh_final_gene_matrix %>%
  rownames_to_column("sampleID_r") %>%
  inner_join(., m_ids_mapped, by = "sampleID_r") %>%
  column_to_rownames("count_names") %>%
  dplyr::select(-c(sampleID_r))



dh_final_gene_matrix_tranposed <- as.data.frame(t(dh_final_gene_matrix))
dh_SLE_norm <- SLE_norm_down_cpgs[rownames(dh_final_gene_matrix_tranposed), ]

dh_common_cols <- intersect(colnames(dh_SLE_norm), colnames(dh_final_gene_matrix_tranposed)) 


dh_SLE_norm <- dh_SLE_norm[, dh_common_cols]
dh_final_gene_matrix_tranposed <- dh_final_gene_matrix_tranposed[, dh_common_cols]





dh_correlation_results <- data.frame(Gene = rownames(dh_SLE_norm), 
                                  Correlation = NA, 
                                  P_value = NA)

# Loop through each gene
for (gene in rownames(dh_SLE_norm)) {
  
  # Extract expression and methylation values for the gene
  expr_values <- as.numeric(dh_SLE_norm[gene, ])
  meth_values <- as.numeric(dh_final_gene_matrix_tranposed[gene, ])
  
  # Perform correlation test (returns both correlation coefficient and p-value)
  cor_test <- cor.test(meth_values, expr_values, method = "pearson")
  
  # Store correlation coefficient and p-value in the dataframe
  dh_correlation_results[dh_correlation_results$Gene == gene, "Correlation"] <- cor_test$estimate
  dh_correlation_results[dh_correlation_results$Gene == gene, "P_value"] <- cor_test$p.value
}

dh_correlation_results$Adj_P_Value <- p.adjust(dh_correlation_results$P_value, method = "BH")

dh_correlation_results$NegLogPadj <- -log10(dh_correlation_results$Adj_P_Value)

dh_correlation_results <- dh_correlation_results %>%
  inner_join(., ensg_mapping, by = c("Gene" = "ID")) 

dh_correlation_results <- dh_correlation_results %>%
  mutate(label = ifelse(NegLogPadj > 20, gene_name, NA))

dh_volcano <- ggplot(dh_correlation_results, aes(x = Correlation, y = NegLogPadj)) +
  geom_point(aes(color = NegLogPadj > 1.3), size = 2) +  
  scale_color_manual(values = c("gray", "#5e3c99")) +   
  geom_text_repel(aes(label = label), 
                  max.overlaps = Inf,  
                  box.padding = 0.8,  
                  point.padding = 0.5, 
                  segment.color = "grey50") +
  labs(x = "Pearson correlation", y = expression(-log[10](P[adj]))) +
  my.theme + 
  theme(legend.position = "none")  +
  scale_x_continuous(limits = c(-1, 1))  

ggsave(paste0(data_dir, "/dh_volcano_M.svg"), dh_volcano, width = 5, height = 5)

correlation_results <- correlation_results %>%
  mutate(dh = ifelse(Gene %in% down_cpgs_genes, "dh", NA))

dh_correlation_results <- dh_correlation_results %>%
  mutate(dh = ifelse(Gene %in% correlation_results$Gene, "dh", NA))

dh_correlation_results <- dh_correlation_results %>%
  mutate(
    dh_status = case_when(
      NegLogPadj > 1.3 & dh == "dh" ~ "Both",
      NegLogPadj > 1.3 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

volcano_interferome_dh <- ggplot(dh_correlation_results, aes(x = Correlation, y = NegLogPadj)) +
  geom_point(aes(color = NegLogPadj > 1.3), size = 2) +
  geom_point(
    data = dh_correlation_results %>% filter(dh_status == "Both"),
    aes(x = Correlation, y = NegLogPadj),
    shape = 21,       
    size = 2.5,       
    color = "black",  
    fill = NA,
    stroke = 1.25
  ) +
    scale_color_manual(values = c("gray", "#5e3c99")) + 
  # Add text labels
  geom_text_repel(
    aes(label = label),
    max.overlaps = Inf,
    box.padding = 0.8,
    point.padding = 0.5,
    segment.color = "grey50"
  ) +
  # Add axis labels and theme
  labs(x = "Pearson correlation", y = expression(-log[10](P[adj]))) +
  my.theme +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(-1, 1))


ggsave(paste0(data_dir, "/dh_volcano_M_interferome_labeled.svg"), volcano_interferome_dh, width = 5, height = 5)

