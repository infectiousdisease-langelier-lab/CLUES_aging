# Here I analyze the CLUES bulk data for female patients

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

setwd("")

# Read cleaned genecounts data -----
counts <- read.csv(
  "",
  row.names = 1,
  check.names = FALSE
)

# Read cleaned metadata -----

cleaned_metadata <- read.csv(
  "", row.names = 1, check.names = FALSE
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





