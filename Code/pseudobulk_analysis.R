# Here I look for genes DE with age based on single cell pseudobulk data

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
library(stringr)
library(RColorBrewer)

setwd("./new_files/")

# Read all data -----

my_files <- list.files()
all_csv <- lapply(my_files,read.csv,header = 1)
## Set the name of each list element to its respective file name. Note full.names = FALSE to get
## only the file names, not the full path.
names(all_csv) <- gsub(".csv","",
                       list.files("./",full.names = FALSE),
                       fixed = TRUE)


pseudo <- all_csv[-c(13)]

# Read metadata -----

metadata <- as.data.frame(all_csv["metadata"])

# Metadata manipulation and cleaning ----

## Drop index column and rename columns

metadata = subset(metadata, select = -c(1))

column_headers <- c("age", "sex", "raceeth", "subject_id", "cc", "group", "study")

colnames(metadata) = c(column_headers)

## We are interested in sex, age, and ethnicity/race
metadata <- metadata[,c("sex","age","raceeth", "subject_id", "group")]

## Drop raceeth = hispanic due to low # of hispanic patients

metadata <- subset(metadata, metadata$raceeth != "Hispanic or Latin American")

## Separate metadata into control and lupus

control_metadata <- subset(metadata,group == 'Control')

control_metadata <- droplevels(control_metadata)
lupus_metadata <- subset(metadata,group == 'Lupus')

## Convert categorial variables into factors

for (i in c("sex","raceeth", "group")) {
  control_metadata[,i] <- as.factor(control_metadata[,i])
  lupus_metadata[,i] <- as.factor(lupus_metadata[,i])
}

# Load GSEA pathways -----

## Reactome gene sets
reactome <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = "C2", subcategory = "CP:REACTOME")
reactome.list <- split(x = reactome$gene_symbol, f = reactome$gs_name)

#QC for pseudobulk data -----

## Input: bulk refers to the pseudobulk data for a given cell type
## Output: the pseudobulk data after quality control has been performed

quality_control <- function(bulk) {
  
  ## Remove samples with fewer than 5,000 genes
  
  bulk <- bulk[, colSums(bulk>0) >= 5000]
  
  ## Keep genes with >=5 counts in >=15% of samples
  
  keep <- rowSums(bulk >= 5) >= 0.15*ncol(bulk)
  bulk <- bulk[keep, ]
}

## Apply quality control to relevant data

pseudo_QC <- lapply(pseudo, quality_control)


# Analysis ---------

## Set up ggplot2 theme

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## Organize pseudobulk data
T8_temp <- as.data.frame(pseudo_QC["annotated_CLUES_T8_naive_pseudobulk_count"])
B_mem <- as.data.frame(pseudo_QC["annotated_CLUES_B_mem_pseudobulk_count"])
B_naive <- as.data.frame(pseudo_QC["annotated_CLUES_B_naive_pseudobulk_count"])
NK <- as.data.frame(pseudo_QC["annotated_CLUES_NK_pseudobulk_count"])
T4_em <- as.data.frame(pseudo_QC["annotated_CLUES_T4_em_pseudobulk_count"])
T4_naive <- as.data.frame(pseudo_QC["annotated_CLUES_T4_naive_pseudobulk_count"])
T4_reg <- as.data.frame(pseudo_QC["annotated_CLUES_T4_reg_pseudobulk_count"])
ncM <- as.data.frame(pseudo_QC["annotated_CLUES_ncM_pseudobulk_count"])
cM <- as.data.frame(pseudo_QC["annotated_CLUES_cM_pseudobulk_count"])
pDC <- as.data.frame(pseudo_QC["annotated_CLUES_pDC_pseudobulk_count"])
Prolif <- as.data.frame(pseudo_QC["annotated_CLUES_Prolif_pseudobulk_count"])
cDC <- as.data.frame(pseudo_QC["annotated_CLUES_cDC_pseudobulk_count"])


# Clean pseudobulk data -----

## Inputs: current_df refers to the pseudobulk data for a given cell type, name refers to the cell type name
## Output: returns the pseudobulk data after it has been cleaned; the row names are
## set to the value of the first column, the first column is removed, HC. is coverted to HC-
clean <- function(current_df, name) {
  
  rownames(current_df) <- current_df[,1]
  current_df = subset(current_df, select = -c(1))
  names(current_df) <- sub(paste('annotated_CLUES', name, 'pseudobulk_count.X', sep = '_'), '', names(current_df))
  names(current_df) <- sub(paste('annotated_CLUES', name, 'pseudobulk_count.', sep = '_'), '', names(current_df))
  names(current_df) <- str_replace(names(current_df), "HC.", "HC-")
  
  return(current_df)
}

# Select the samples from pseudobulk of interest for which metadata exists in the cleaned metadata dataframe -----

## Input: current_df refers to the pseudobulk data for a given cell type
## Output: A list containing the metadata for SLE patients and controls, pseudobulk counts
## for SLE patients and controls, and the number of samples for both categories for that cell type

select <- function(cleaned_df) {
  
  ## Work with control samples first
  
  metadata_control <- control_metadata %>%
    subset((subject_id %in% colnames(cleaned_df)))
  counts_control <- cleaned_df[,metadata_control$subject_id]
  stopifnot(colnames(counts_control)==metadata_control$subject_id)
  
  ## Now, work with lupus samples
  
  metadata_lupus <- lupus_metadata %>%
    subset((subject_id %in% colnames(cleaned_df)))
  counts_lupus <- cleaned_df[,metadata_lupus$subject_id]
  stopifnot(colnames(counts_lupus)==metadata_lupus$subject_id)
  
  
  relevant_data <- list(metadata_control, metadata_lupus, counts_control, counts_lupus, nrow(metadata_control), nrow(metadata_lupus))
  names(relevant_data) <- c("metadata_control", "metadata_lupus", "counts_control", "counts_lupus", "ncont", "nlup")
  return(relevant_data)
}

# Design matrix and DE analysis -----

## Input: list_of_data is a list containing the metadata for SLE patients and controls, pseudobulk counts
## for SLE patients and controls, and the number of samples for both categories for that cell type

## Output: A list containing the output of the DE analysis for the cell-specific pseudobulk
## data for SLE patients and controls

analysis <- function(list_of_data) {
  
  ## Set up design matrices 
  
  design_lupus <- model.matrix(~age + sex + raceeth,
                               data = list_of_data[["metadata_lupus"]])
  
  design_control <- model.matrix(~age + sex + raceeth,
                                 data = list_of_data[["metadata_control"]])
  
  ## limma-voom
  
  vwts_control <- voom(list_of_data[["counts_control"]], 
                       design = design_control,
                       normalize.method = "quantile",
                       plot = T)
  
  vwts_lupus <- voom(list_of_data[["counts_lupus"]], 
                     design = design_lupus,
                     normalize.method = "quantile",
                     plot = T)
  
  
  ## DE analysis
  
  vfit_control <- lmFit(vwts_control)
  vfit_control <- eBayes(vfit_control)
  
  vfit_lupus <- lmFit(vwts_lupus)
  vfit_lupus <- eBayes(vfit_lupus)
  
  DE_results <- list(vfit_control, vfit_lupus)
  names(DE_results) <- c("vfit_control", "vfit_lupus")
  return(DE_results)
  
  
}

# Extract top tables from linear regression analysis -----

## Input: DE_results is a list containing the output of the DE analysis for the cell-specific pseudobulk
## data for SLE patients and controls

## Output: A list "tables." top_control is the result of toptable on the DE analysis for the cell-specific
## pseudobulk data for controls. sig_control is the same data but only for significantly DE genes. sig2_control
## is similar to sig_control, but with a less stringent threshold for DE genes nsigcont
## is the number of rows in Sig_control, which can be interepreted as the number of DE genes. nsig2cont
## is the number of DE genes with a less stringent p.value threshold

## The same naming scheme applies to _lupus

extract <- function(DE_results) {
  
  ## Toptable for significant genes
  
  top_control <- topTable(DE_results[["vfit_control"]], coef = "age", sort.by = "none",
                          number = Inf, p.value = 1)
  
  sig_control <- topTable(DE_results[["vfit_control"]], coef = "age", sort.by = "none",
                          number = Inf, p.value = 0.05)
  
  top_lupus <- topTable(DE_results[["vfit_lupus"]], coef = "age", sort.by = "none",
                        number = Inf, p.value = 1)
  
  sig_lupus <- topTable(DE_results[["vfit_lupus"]], coef = "age", sort.by = "none",
                        number = Inf, p.value = 0.05)
  
  sig2_lupus <- topTable(DE_results[["vfit_lupus"]], coef = "age", sort.by = "none",
                         number = Inf, p.value = 0.1)
  
  sig2_control <- topTable(DE_results[["vfit_control"]], coef = "age", sort.by = "none",
                           number = Inf, p.value = 0.1)
  
  tables <- list(top_control, sig_control, top_lupus, sig_lupus, nrow(sig_control), nrow(sig_lupus), nrow(sig2_control), nrow(sig2_lupus))
  names(tables) <- c("top_control", "sig_control", "top_lupus", "sig_lupus", "nsigcont", "nsiglupus", "nsig2cont", "nsig2lupus")
  return(tables)
}

#Generate gsea results

## Inputs: topTables is a list with the following items: top_control is the result of toptable on the DE analysis for the cell-specific
## pseudobulk data for controls. sig_control is the same data but only for significantly DE genes. sig2_control
## is similar to sig_control, but with a less stringent threshold for DE genes nsigcont
## is the number of rows in Sig_control, which can be interepreted as the number of DE genes. nsig2cont
## is the number of DE genes with a less stringent p.value threshold
## The same naming scheme applies to _lupus.

## name is the cell type name
## number_lup is the total number of lupus samples passing QC
## number_control is the total number of control samples passing QC
## nsigl and nsigc are the total number of DE genes for the pseudobulk of a given cell type
## where the threshold is p < 0.05
## sig2l and sig2c are the total number of DE genes for the pseudobulk of a given cell type 
## where the threshold is p < 0.1

## Outputs: a list, "plots" with nsigl, nsigc, sig2l, sig2c, which are the same as inputs
## lupus.gsea, control.gsea which are the gsea results for the lupus and control groups
## for a given cell type's pseudobulk
## which is the same as the input topTables

plot <- function(topTables, name, number_lup, number_control, nsigl, nsigc, sig2l, sig2c) {
  set.seed(19) ## Reproducibility
  control_gene.ranks <- -log10(topTables[["top_control"]]$P.Value) * sign(topTables[["top_control"]]$logFC)
  names(control_gene.ranks) <- rownames(topTables[["top_control"]])
  control_gene.ranks <- sort(control_gene.ranks, decreasing = TRUE)
  
  lupus_gene.ranks <- -log10(topTables[["top_lupus"]]$P.Value) * sign(topTables[["top_lupus"]]$logFC)
  names(lupus_gene.ranks) <- rownames(topTables[["top_lupus"]])
  lupus_gene.ranks <- sort(lupus_gene.ranks, decreasing = TRUE)
  
  control.gsea <- fgseaMultilevel(
    pathways = reactome.list,
    stats = control_gene.ranks,
    minSize = 15,
    maxSize = 500,
    nproc = 1
  )
  lupus.gsea <- fgseaMultilevel(
    pathways = reactome.list,
    stats = lupus_gene.ranks,
    minSize = 15,
    maxSize = 500,
    nproc = 1
  )
  unique_pathways <- c(
    "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
    "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
    "REACTOME_INTERFERON_GAMMA_SIGNALING",
    "REACTOME_INTERFERON_SIGNALING",
    "REACTOME_DDX58_IFIH1_MEDIATED_INDUCTION_OF_INTERFERON_ALPHA_BETA",
    "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
    "REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING",
    "REACTOME_SIGNALING_BY_INTERLEUKINS",
    "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS",
    "REACTOME_TLR3_MEDIATED_TICAM1_DEPENDENT_PROGRAMMED_CELL_DEATH", "REACTOME_INTERLEUKIN_6_FAMILY_SIGNALING", "REACTOME_INTERLEUKIN_6_SIGNALING", "REACTOME_INTERLEUKIN_17_SIGNALING",
    "REACTOME_INTERLEUKIN_9_SIGNALING"
  )
  lupus.gsea <- lupus.gsea[lupus.gsea$pathway %in% unique_pathways, ]
  control.gsea <- control.gsea[control.gsea$pathway %in% unique_pathways, ]
  
  lupus.gsea$pathway <- sub('^REACTOME_', '', lupus.gsea$pathway)
  lupus.gsea$pathway <- gsub("_", " ", lupus.gsea$pathway)
  
  control.gsea$pathway <- sub('^REACTOME_', '', control.gsea$pathway)
  control.gsea$pathway <- gsub("_", " ", control.gsea$pathway)
  
  plots <- list(nsigl, nsigc, sig2l, sig2c, lupus.gsea, control.gsea, topTables)
  names(plots) <- c("0.05 significant lupus", "0.05 significant control" ,"0.1 significant_lupus", "0.1 significant_cont", paste0(name, "_lupus"), paste0(name, "_control"))
  return(plots)
  
}

# Call everything -----

## Input: The pseudobulk dataframe for a cell type along with its name
## Output: The output of the plot function (see above)
## When the function is complete with a cell type, it prints an output accordingly

complete <- function(original_df, name) {
  
  cleaned <- clean(original_df, name) ## Clean the dataframe
  selected <- select(cleaned) ## Select only the rows in counts which pass QC
  clupus <- selected[["nlup"]] ## Determine the number of lupus samples passing QC
  ccont <- selected[["ncont"]] ## Determine the number of control samples passing QC
  dma <- analysis(selected) ## Create design matrices for lupus and control
  top <- extract(dma) ## Extract top tables for lupus and control
  sigl <- top[["nsiglupus"]] ## Determine the number of DE genes for lupus patients
  sigc <- top[["nsigcont"]] ## Determine the number of DE genes for control patients
  sig2c <- top[["nsig2cont"]]
  sig2l <- top[["nsig2lupus"]]
  
  finish <- plot(top, name, clupus, ccont, sigl, sigc, sig2l, sig2c) ## Plot
  print(paste0("Done with ", name, "!"))
  return(finish)
  
}

# Apply to cells -----

T8_naive_plots <- complete(T8_temp, "T8_naive")
B_naive_plots <- complete(B_naive, "B_naive")
NK_plots <- complete(NK, "NK")
T4_em_plots <- complete(T4_em, "T4_em")
B_em_plots <- complete(B_mem, "B_mem")
T4_naive_plots <- complete(T4_naive, "T4_naive")
T4_reg_plots <- complete(T4_reg, "T4_reg")
ncM_plots <- complete(ncM, "ncM")
cM_plots <- complete(cM, "cM")

## To generate the cDC data, remove the sex variable from the formula in lines 165 and 168
cDC_plots <- complete(cDC, "cDC") #removed sex

# Plot for cell types of interest -----

## Note: the data from these output dataframes is used for the dotplots in the script dotplot.R

T8_CLUES.gsea <- T8_naive_plots[[5]]
T8_RS.gsea <- T8_naive_plots[[6]]
T8_CLUES.gsea$CT <- "T8_naive_CLUES"
T8_RS.gsea$CT <- "T8_naive_RS"

B_CLUES.gsea <- B_naive_plots[[5]]
B_RS.gsea <- B_naive_plots[[6]]
B_CLUES.gsea$CT <- "B_naive_CLUES"
B_RS.gsea$CT <- "B_naive_RS"

NK_CLUES.gsea <- NK_plots[[5]]
NK_RS.gsea <- NK_plots[[6]]
NK_CLUES.gsea$CT <- "NK_CLUES"
NK_RS.gsea$CT <- "NK_RS"

T4_em_CLUES.gsea <- T4_em_plots[[5]]
T4_em_RS.gsea <- T4_em_plots[[6]]
T4_em_CLUES.gsea$CT <- "T4_em_CLUES"
T4_em_RS.gsea$CT <- "T4_em_RS"

## Combine the subsets

three <- rbind(B_CLUES.gsea, B_RS.gsea, NK_CLUES.gsea, NK_RS.gsea, T4_em_CLUES.gsea, T4_em_RS.gsea)

pathway <- ggplot(data = three, aes(x=pathway, y=NES)) +
  geom_segment(aes(xend=reorder(pathway,NES)), yend=0) +
  geom_point(aes(fill = ifelse(padj<0.05, NES > 0, '6')),
             pch=21, size=2, stroke=1) +
  scale_fill_manual(name = "As Age Increases,
  the Pathway is", values = c('TRUE' = "#e66101", 'FALSE' = "#5e3c99", '6' = 'white'), drop = FALSE,
                    breaks = c('TRUE', 'FALSE'),
                    labels = c("TRUE" = "Upregulated","FALSE" = "Downregulated")) +
  labs(x = element_blank(), y="NES") +
  coord_flip() + my.theme + theme(panel.border = element_rect(color = "black", fill = NA), axis.text.y = element_text(size = 11), legend.position = "none", legend.box = "vertical",
                                  legend.spacing.x = unit(1, 'mm'), legend.spacing.y = unit(1, 'mm'), legend.margin=margin(0,0,0,0),
                                  legend.box.margin=margin(-10,-10,-10,-10), legend.text = element_text(size=10)) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) + ylim(c(-3.0, 3.0)) + geom_hline(yintercept = 0) + facet_wrap(~ CT, ncol = 6)
