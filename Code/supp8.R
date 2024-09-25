# Here I look for genes DE with age

# Are main-idea comments

## Are details/steps

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

# Read genecounts data -----
bulk <- read.csv(
  "",
  row.names = 1,
  check.names = FALSE
)

# Bulk data cleaning ----- 

gene <- rownames(bulk)

#QC -----

## Remove samples with fewer than 10,000 genes
bulk <- bulk[, colSums(bulk>0) >= 10000]

## Keep genes with >=10 counts in >=20% of samples
keep <- rowSums(bulk >= 10) >= 0.2*ncol(bulk)
bulk <- bulk[keep, ]


# Read metadata -----

metadata <- read.csv(
  ""
)

# Metadata manipulation and cleaning ----

## We are interested in these variables
metadata <- metadata[,c("subject_id","sex","age","raceeth", "Group")]


## Age > 50
metadata <- subset(metadata, metadata$age > 65) # change based on age filter
metadata <- subset(metadata, metadata$sex == "female")


## Convert categorial variables into factors

for (i in c("raceeth", "Group")) {
  metadata[,i] <- as.factor(metadata[,i])
}

## Set up sample id 

cleaned_metadata <- metadata

rownames(cleaned_metadata) <- cleaned_metadata$subject_id


# Load GSEA pathways -----
## Reactome
reactome <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = "C2", subcategory = "CP:REACTOME")
reactome.list <- split(x = reactome$gene_symbol, f = reactome$gs_name)


# Prepare for DE -----

## Select the samples from bulk for which metadata exists in the cleaned metadata dataframe
cleaned_metadata <- cleaned_metadata %>%
  subset((subject_id %in% colnames(bulk)))

counts <- bulk[,cleaned_metadata$subject_id]
stopifnot(colnames(counts)==cleaned_metadata$subject_id)


##Controlling for ethnicity, sex
design <- model.matrix(~Group + raceeth,
                       data = cleaned_metadata)

## limma-voom
vwts <- voom(counts, 
             design = design,
             normalize.method = "quantile",
             plot = T)

# DE analysis ------
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)

## Set up ggplot2 theme

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )



## Extract genes
top_DE <- topTable(vfit, coef = "GroupLupus", sort.by = "none",
                   number = Inf, p.value = 1)

top_DE$Significance <- ifelse(top_DE$adj.P.Val < 0.05, "Significant (FDR < 0.05)", "Not Significant")

top_DE$diffexpressed <- "Not significant"
top_DE$diffexpressed[top_DE$logFC > 0 & top_DE$adj.P.Val < 0.05] <- "Significantly upregulated"
top_DE$diffexpressed[top_DE$logFC < 0 & top_DE$adj.P.Val < 0.05] <- "Significantly downregulated"
top_DE$gene_name <- rownames(top_DE)
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

modified_list <- c(
  "Antigen activates BCR leading to second messengers",
  "Class I MHC antigen processing presentation",
  "DDX58 IFIH1 induction of IFN alpha beta",
  "IFN alpha beta signaling",
  "IFN gamma signaling",
  "IFN signaling",
  "IL 17 signaling",
  "IL 1 family signaling",
  "MHC class II antigen presentation",
  "Signaling by interleukins"
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
# Replace the pathway column in reactome.gsea with the modified list
reactome.gsea2 <- reactome.gsea[reactome.gsea$pathway %in% unique_pathways2, ]
reactome.gsea2$pathway2 <- modified_list

# change based on age

young <- ggplot(data = reactome.gsea2, aes(x=reorder(pathway2,NES), y=NES)) +
  geom_segment(aes(xend=reorder(pathway2,NES)), yend=0) +
  geom_point(aes(fill = ifelse(padj<0.05, NES > 0, '6')),
             pch=21, size=4, stroke=1) +
  scale_fill_manual(name = "", values = c('TRUE' = "#5e3c99", 'FALSE' = "#ca0020", '6' = 'white'), drop = FALSE,
                    breaks = c('TRUE', 'FALSE'), 
                    labels = c("TRUE" = "Upregulated with SLE","FALSE" = "Downregulated with SLE")) +
  labs(x = element_blank(), element_blank(), 
       title=element_blank())  +
  coord_flip() + my.theme + ylim(c(-3.0, 3.0)) + geom_hline(yintercept = 0)

not_as_young <- ggplot(data = reactome.gsea2, aes(x=reorder(pathway2,NES), y=NES)) +
  geom_segment(aes(xend=reorder(pathway2,NES)), yend=0) +
  geom_point(aes(fill = ifelse(padj<0.05, NES > 0, '6')),
             pch=21, size=4, stroke=1) +
  scale_fill_manual(name = "", values = c('TRUE' = "#5e3c99", 'FALSE' = "#ca0020", '6' = 'white'), drop = FALSE,
                    breaks = c('TRUE', 'FALSE'), 
                    labels = c("TRUE" = "Upregulated with SLE","FALSE" = "Downregulated with SLE")) +
  labs(x = element_blank(), element_blank(), 
       title=element_blank())  +
  coord_flip() + my.theme + ylim(c(-3.0, 3.0)) + geom_hline(yintercept = 0)

slightly_older <- ggplot(data = reactome.gsea2, aes(x=reorder(pathway2,NES), y=NES)) +
  geom_segment(aes(xend=reorder(pathway2,NES)), yend=0) +
  geom_point(aes(fill = ifelse(padj<0.05, NES > 0, '6')),
             pch=21, size=4, stroke=1) +
  scale_fill_manual(name = "", values = c('TRUE' = "#5e3c99", 'FALSE' = "#ca0020", '6' = 'white'), drop = FALSE,
                    breaks = c('TRUE', 'FALSE'), 
                    labels = c("TRUE" = "Upregulated with SLE","FALSE" = "Downregulated with SLE")) +
  labs(x = element_blank(), element_blank(), 
       title=element_blank())  +
  coord_flip() + my.theme + ylim(c(-3.0, 3.0)) + geom_hline(yintercept = 0)

young_volc <- ggplot(top_DE, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = diffexpressed), size = 0.8) +
  scale_color_manual(name = "", values = c("grey", "#ca0020", "#5e3c99")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle(element_blank()) + xlab("Slope") +
  ylab(expression(-log[10](P[adj]))) + my.theme + theme(legend.position = "none") 

mid_volc <- ggplot(top_DE, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = diffexpressed), size = 0.8) +
  scale_color_manual(name = "", values = c("grey", "#ca0020", "#5e3c99")) +
  theme_bw(base_size = 12) +
  ggtitle(element_blank()) + xlab("Slope") +
  ylab(expression(-log[10](P[adj]))) + my.theme + theme(legend.position = "none") 

old_volc <- ggplot(top_DE, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = diffexpressed), size = 0.8) +
  scale_color_manual(name = "", values = c("grey", "#ca0020", "#5e3c99")) +
  theme_bw(base_size = 12) +
  ggtitle(element_blank()) + xlab("Slope") +
  ylab(expression(-log[10](P[adj]))) + my.theme + theme(legend.position = "none") 







######################### ========

norm <- edgeR::cpm(counts, log=TRUE)


expression_df <- data.frame(
  sample_id = names(norm["ISG15", ]),
  ISG15 = as.numeric(norm["ISG15", ]),
  OAS3 = as.numeric(norm["OAS3", ])
)
expression_data <- merge(cleaned_metadata, expression_df, by.x = "subject_id", by.y = "sample_id")

# Filter out "Hispanic" raceeth if needed
expression_data <- expression_data %>%
  filter(raceeth != "Hispanic or Latin American")


