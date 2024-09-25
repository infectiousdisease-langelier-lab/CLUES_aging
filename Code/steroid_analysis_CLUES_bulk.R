# Here I look for genes DE with age in bulk RNA-seq data for CLUES patients, this time controllling for steroid dosage

# Are main-idea commments

## Are details/steps

library(dplyr)
library(tibble)
library(limma)
library(DESeq2)
library(tidyr)
library(fgsea)
library(ggplot2)
library(ggrepel)

# Read genecounts data -----
bulk <- read.csv(
  "./ImmGeneCounts.csv",
  row.names = 1,
  check.names = FALSE
)

# Bulk data cleaning ----- 

## Save gene_name column as row data
rData <- bulk[,"gene_name",drop=FALSE]
bulk <- bulk[,-c(1)]

## Only keep samples with "DNA" in their names
bulk <- bulk[,grepl("DNA",colnames(bulk))]
## Remove samples with "-1" in their names
bulk <- bulk[,!grepl("-1",colnames(bulk), fixed=TRUE)]

# Read metadata -----

all_immune <- read.csv(
  "./bulk_cohort_metadata_022724.csv"
)

## Metadata splitting by time

t1_meta <- subset(all_immune, substr(sampleid, nchar(sampleid)-1, nchar(sampleid)) == "01", select = -steroral3mohigh_t3)
t3_meta <- subset(all_immune, substr(sampleid, nchar(sampleid)-1, nchar(sampleid)) == "11", select = -steroral3mohigh_t1)

# Metadata manipulation and cleaning ----

## We are interested in SterOral3moHigh; drop rows without SterOral3moHigh

t1_meta <- t1_meta %>% drop_na(steroral3mohigh_t1)
t3_meta <- t3_meta %>% drop_na(steroral3mohigh_t3)

### leads to four samples being dropped from t1

## Add 2 to the age of t3 samples as age data was collected at t1

t3_meta <- t3_meta %>% mutate(age = age + 2)

# Prepare for merge -----
dataframes <- list(t1_meta, t3_meta)

## Rename the steroid column in each dataframe to allow for merging
RenameIDs = function(x){
  names(x)[7] <- "steroid"
  return(x)} 

dataframes <- lapply(dataframes, RenameIDs)

t1_meta_cleaned <- as.data.frame(dataframes[1])
t3_meta_cleaned <- as.data.frame(dataframes[2])

meta_final <- rbind(t1_meta_cleaned, t3_meta_cleaned)

# Final cleaning and preparation for analysis

## Convert categorial variables into factors

for (i in c("female","raceeth", "steroid")) {
  meta_final[,i] <- as.factor(meta_final[,i])
}

cleaned_metadata <- meta_final

# Prepare for DE -----

## Select the samples from bulk for which metadata exists in the cleaned metadata dataframe
cleaned_metadata <- cleaned_metadata %>%
  subset((sampleid %in% colnames(bulk)))

counts <- bulk[,cleaned_metadata$sampleid]
stopifnot(colnames(counts)==cleaned_metadata$sampleid)

#QC -----

## Remove samples with fewer than 10,000 genes
counts <- counts[, colSums(counts>0) >= 10000]

## Keep genes with >=10 counts in >=20% of samples
keep <- rowSums(counts >= 10) >= 0.2*ncol(counts)
counts <- counts[keep, ]
stopifnot(cleaned_metadata$sampleid==colnames(counts))


## Make design matrix -----

##Controlling for ethnicity, sex, oral steroid dosage
design <- model.matrix(~age + female + raceeth + steroid,
                       data = cleaned_metadata)

## limma-voom
vwts <- voom(counts, 
             design = design,
             normalize.method = "quantile",
             plot = T)

# DE analysis ------
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)

## Extract top DE genes
top_DE <- topTable(vfit, coef = "age", sort.by = "none",
                   number = Inf, p.value = 1)

# Set up ggplot2 theme -----

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

# Load GSEA pathways -----
## Reactome
reactome <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = "C2", subcategory = "CP:REACTOME")
reactome.list <- split(x = reactome$ensembl_gene, f = reactome$gs_name)

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




