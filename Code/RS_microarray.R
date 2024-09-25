# Here I look for gene expression with age in healthy individuals using microarray data

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
library(GEOquery)
library(data.table)

# Read cleaned microarray data -----

expression = fread("./cleaned_RS_microarray.csv", header = TRUE)


## Map Probe IDs to gene symbol

gse <- getGEO(GEO = "GSE33828", GSEMatrix = TRUE)

## Fetch feature data to complete mapping

feature.data <- gse$GSE33828_series_matrix.txt.gz@featureData@data
feature.data <- feature.data[,c(1,6)]

genex <- expression %>% 
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data, by = 'ID')

genex <- as.data.frame(genex)

nams = genex$ILMN_Gene
rownames(genex) <- make.names(nams, unique = TRUE)
genex <- genex[,-c(1)]
genex <- dplyr::select(genex, -c("ILMN_Gene"))


# Read metadata -----


RS <- read.csv(
  "./cleaned_RS_metadata.csv"
)

RS <- as.data.frame(RS)

rownames(RS) <- RS$RS3.ID

# Metadata manipulation ----

## Convert categorial variables into factors

for (i in c("sex")) {
  RS[,i] <- as.factor(RS[,i])
}


# DE ------

##Controlling for sex
set.seed(19)
design <- model.matrix(~age + sex,
                       data = RS)


## DE analysis ------
fit <- lmFit(genex, design)
fit <- eBayes(fit)



## Extract top DE genes
top_DE <- topTable(fit, coef = "age", sort.by = "none",
                   number = Inf, p.value = 1)

# Volcano Plot -----

## Set up ggplot2 theme

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## Plotting and wrangling

top_DE$Significance <- ifelse(top_DE$adj.P.Val < 0.05, "Significant (FDR < 0.05)", "Not Significant")

top_DE$gene_name = rownames(top_DE)

top_DE$diffexpressed <- "Not significant"
top_DE$diffexpressed[top_DE$logFC > 0 & top_DE$adj.P.Val < 0.05] <- "Significantly upregulated"
top_DE$diffexpressed[top_DE$logFC < 0 & top_DE$adj.P.Val < 0.05] <- "Significantly downregulated"

volcano <- ggplot(top_DE, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = diffexpressed), size = 0.8) +
  scale_color_manual(name = "", values = c("grey", "#0571b0", "#ca0020")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("Healthy individuals") + xlab("Slope") +
  ylab("-log10(P   )") + my.theme + theme(axis.text.y = element_text(size = 12), legend.position = "none", legend.box = "vertical",
                                          legend.spacing.x = unit(1, 'mm'), legend.spacing.y = unit(1, 'mm'), legend.margin=margin(0,0,0,0),
                                          legend.box.margin=margin(-10,-10,-10,-10), legend.text = element_text(size = 12)) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) + scale_x_continuous(breaks = c(-0.15, 0.0, 0.15), limits = c(-0.15, 0.15)) + my.theme

# Regression plots for expression of genes of interest with age -----

set.seed(19)

expression_values <- as.matrix(genex)

OAS3 <- data.frame(RS$age, expression_values["OAS3",])
ISG15 <- data.frame(RS$age, expression_values["ISG15",])


goi <- list(ISG15, OAS3)
names(goi) <- c("ISG15", "OAS3")


cols = c('age', 'expression')

goi <- lapply(goi, setNames, cols)


ISG15_RS_ex <- ggplot(data = goi[[1]],aes(x = age, y = expression)) + geom_jitter(size = 0.8) + geom_smooth(method = 'lm', se = TRUE, color = "#ca0020", fill = "#ca0020") + 
  ggtitle(label = "ISG15 Expression vs. Age") + labs(y = "Normalized expression", x = "Age") + my.theme

OAS3_RS_ex <- ggplot(data = goi[[2]],aes(x = age, y = expression)) + geom_jitter(size = 0.8) + geom_smooth(method = 'lm', se = TRUE, color = "#ca0020", fill = "#ca0020") + 
  ggtitle(label = "OAS3 Expression vs. Age") + labs(y = "Normalized expression", x = "Age") + my.theme


# Load GSEA pathways -----
## Reactome
reactome <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = "C2", subcategory = "CP:REACTOME")
reactome.list <- split(x = reactome$ensembl_gene, f = reactome$gs_name)



# fGSEA ------

## Convert gene names to ensembl ids
library("AnnotationDbi")
library("org.Hs.eg.db")
top_DE$ensid = mapIds(org.Hs.eg.db,
                      keys=rownames(top_DE), 
                      column="ENSEMBL",
                      keytype="SYMBOL",
                      multiVals="first")

keep_DE <- na.omit(top_DE)

gene.ranks <- -log10(keep_DE$P.Value) * sign(keep_DE$logFC)
names(gene.ranks) <- keep_DE$ensid
gene.ranks <- sort(gene.ranks, decreasing = TRUE)


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
reactome.gsea2 <- reactome.gsea[reactome.gsea$pathway %in% unique_pathways2, ]

pathway_plots <- ggplot(data = reactome.gsea2, aes(x=reorder(pathway,NES), y=NES)) +
  geom_segment(aes(xend=reorder(pathway,NES), yend=0)) +
  geom_point(aes(fill = ifelse(padj<0.05, NES > 0, '6')),
             pch=21, size=4, stroke=1) +
  scale_fill_manual(drop = "false", name = "As Age Increases, 
  the Pathway is", values = c('TRUE' = "#ca0020", 'FALSE' = "#0571b0", '6' = 'white'),
                    labels = c("TRUE" = "Upregulated","FALSE" = "Downregulated"), 
                    breaks = c('FALSE', 'TRUE')) + 
  labs(x = element_blank(), y="NES", 
       title="Healthy Controls")  +
  coord_flip() + my.theme + 
  theme(panel.border = element_rect(color = "black", fill = NA), axis.text.y = element_text(size = 11), legend.position = "none", legend.box = "vertical",
        legend.spacing.x = unit(1, 'mm'), legend.spacing.y = unit(1, 'mm'), legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10), legend.text = element_text(size = 12)) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) + ylim(c(-3.0, 3.0)) + geom_hline(yintercept = 0)



