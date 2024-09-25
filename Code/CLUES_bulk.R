# Here I analyze the CLUES bulk data

# Are main-idea commments

## Are details/steps

# Load packages ----

library(dplyr)
library(tibble)
library(limma)
library(DESeq2)
library(tidyr)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(ggeffects)

setwd("/Users/rithwikn/Documents/CZBiohub_2023/Lupus/")

# Read cleaned genecounts data -----
counts <- read.csv(
  "../../cleaned_CLUES_bulkcounts.csv",
  row.names = 1,
  check.names = FALSE
)

# Read cleaned metadata -----

cleaned_metadata <- read.csv(
  "../../cleaned_CLUES_metadata.csv", row.names = 1, check.names = FALSE
)

## Convert categorial variables into factors

for (i in c("female","raceeth")) {
  cleaned_metadata[,i] <- as.factor(cleaned_metadata[,i])
}


# Load GSEA pathways -----
## Reactome gene sets
reactome <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = "C2", subcategory = "CP:REACTOME")
reactome.list <- split(x = reactome$ensembl_gene, f = reactome$gs_name)

# SLEDAI regression plots -----

## Set up ggplot2 theme

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## Plot SLEDAI regression

SLEDAI <- lm(sledaiscore ~ (age + female + raceeth), data = cleaned_metadata)
SLEDAI_slope <- summary(SLEDAI)$coefficients["age", "Estimate"]
SLEDAI_pval <- summary(SLEDAI)$coefficients["age", "Pr(>|t|)"]
SLEDAI_CI <- ggpredict(SLEDAI, terms = "age")
plot <- ggplot(cleaned_metadata, aes(x=age, y=sledaiscore)) + geom_point(color = "black", size = 2) +
  labs(x = "Age", y = "SLEDAI Score")
plot <- plot +
  geom_line(aes(x = x, y = predicted), linewidth = 2, data = SLEDAI_CI, color = "#5e3c99") +
  geom_ribbon(inherit.aes = FALSE,
              aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
              data = SLEDAI_CI,
              alpha = 0.1
  ) +
  ggtitle("SLEDAI Score vs. Age") + my.theme

# DE analysis -----

## Make design matrix, controlling for ethnicity, sex
design <- model.matrix(~age + female + raceeth,
                       data = cleaned_metadata)

## limma-voom normalization
vwts <- voom(counts, 
             design = design,
             normalize.method = "quantile",
             plot = T)

# DE calculation -----
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)


## Extract genes
top_DE <- topTable(vfit, coef = "age", sort.by = "none",
                   number = Inf, p.value = 1)

top_DE$Significance <- ifelse(top_DE$adj.P.Val < 0.05, "Significant (FDR < 0.05)", "Not Significant")

top_DE$diffexpressed <- "Not significant"
top_DE$diffexpressed[top_DE$logFC > 0 & top_DE$adj.P.Val < 0.05] <- "Significantly upregulated"
top_DE$diffexpressed[top_DE$logFC < 0 & top_DE$adj.P.Val < 0.05] <- "Significantly downregulated"

# Volcano plot -----

volcano <- ggplot(top_DE, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = diffexpressed), size = 0.8) +
  scale_color_manual(name = "", values = c("grey", "#5e3c99", "#e66101")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("CLUES") + xlab("Slope") +
  ylab("-log10(P  )") + my.theme + theme(axis.text.y = element_text(size = 11), legend.position = "none", legend.box = "vertical",
                                         legend.spacing.x = unit(1, 'mm'), legend.spacing.y = unit(1, 'mm'), legend.margin=margin(0,0,0,0),
                                         legend.box.margin=margin(-10,-10,-10,-10), legend.text = element_text(size=12)) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) + xlim(c(-0.075, 0.075))


# Regression plots for expression of genes of interest with age -----

norm <- edgeR::cpm(counts, log=TRUE)

OAS3 <- data.frame(cleaned_metadata$age, norm["ENSG00000111331",])
ISG15 <- data.frame(cleaned_metadata$age, norm["ENSG00000187608",])

goi <- list(ISG15, OAS3)
names(goi) <- c("ISG15", "OAS3")

cols = c('age', 'expression')

goi <- lapply(goi, setNames, cols)

ISG15_CLUES_ex <- ggplot(data = goi[[1]],aes(x = age, y = expression)) + geom_jitter(size = 0.8) + geom_smooth(method = 'lm', se = TRUE, color = "#5e3c99", fill = "#5e3c99") + 
  ggtitle(label = "ISG15 Expression vs. Age") + labs(y = "Normalized expression", x = "Age") + my.theme

OAS3_CLUES_ex <- ggplot(data = goi[[2]],aes(x = age, y = expression)) + geom_jitter(size = 0.8) + geom_smooth(method = 'lm', se = TRUE, color = "#5e3c99", fill = "#5e3c99") + 
  ggtitle(label = "OAS3 Expression vs. Age") + labs(y = "Normalized expression", x = "Age") + my.theme






# fGSEA ------

gene.ranks <- -log10(top_DE$P.Value) * sign(top_DE$logFC)
names(gene.ranks) <- rownames(top_DE)
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

## Focusing on reactome pathways

set.seed(19) ## Reproducibility
reactome.gsea <- fgseaMultilevel(
  pathways = reactome.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

unique_pathways2 <- c(
  "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_INTERFERON_SIGNALING",
  "REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING",
  "REACTOME_DDX58_IFIH1_MEDIATED_INDUCTION_OF_INTERFERON_ALPHA_BETA",
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
  "REACTOME_SIGNALING_BY_INTERLEUKINS",
  "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS",
  "REACTOME_INTERLEUKIN_6_FAMILY_SIGNALING", "REACTOME_INTERLEUKIN_17_SIGNALING"
)

## Replace the pathway column in reactome.gsea with the modified list
reactome.gsea2 <- reactome.gsea[reactome.gsea$pathway %in% unique_pathways2, ]

pathway_plot <- ggplot(data = reactome.gsea2, aes(x=reorder(pathway,NES), y=NES)) +
  geom_segment(aes(xend=reorder(pathway,NES)), yend=0) +
  geom_point(aes(fill = ifelse(padj<0.05, NES > 0, '6')),
             pch=21, size=4, stroke=1) +
  scale_fill_manual(name = "As Age Increases,
  the Pathway is", values = c('TRUE' = "#e66101", 'FALSE' = "#5e3c99", '6' = 'white'), drop = FALSE,
                    breaks = c('TRUE', 'FALSE'), 
                    labels = c("TRUE" = "Upregulated","FALSE" = "Downregulated")) +
  labs(x = element_blank(), y="NES", 
       title="CLUES")  +
  coord_flip() + my.theme + theme(panel.border = element_rect(color = "black", fill = NA), axis.text.y = element_text(size = 11), legend.position = "none", legend.box = "vertical",
                                  legend.spacing.x = unit(1, 'mm'), legend.spacing.y = unit(1, 'mm'), legend.margin=margin(0,0,0,0),
                                  legend.box.margin=margin(-10,-10,-10,-10), legend.text = element_text(size=10)) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) + ylim(c(-3.0, 3.0)) + geom_hline(yintercept = 0)





