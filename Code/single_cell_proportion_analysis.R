# Here I analyze the slopes of single cell data
# The slopes are change in cell proportion over an increase of age by one unit
# I will be comparing CLUES patients with Healthy Controls

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
library(janitor)
library(ggpubr)


# Read slope data -----
raw <- read.csv(
  "/Users/rithwikn/Documents/CZBiohub_2023/SC/Proportion_analysis/raw_slopes.csv",
  row.names = 1,
  check.names = FALSE
)

rownames(raw) <- raw$cell


# Set up ggplot2 theme -----

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

# Prepare for plotting -----

## Select relevant columns

raw_keep <- raw
raw_slopes <- raw_keep[,-c(3,4,6,7)]

## Format data for plot

ready_raw <- gather(raw_slopes, group, slope, 2:3)
ready_raw <- transform(
  ready_raw, pval= ifelse(group=="lupus_slope", raw_keep$lupus_corrected_pval, raw_keep$control_corrected_pval))
graph_raw <- ready_raw

graph_raw$slope <- as.numeric(as.character(graph_raw$slope))
graph_raw$pval <- as.numeric(as.character(graph_raw$pval))

## Change to change in % per ten years of age
graph_raw$adjust_slope = graph_raw$slope * 1000

## Rename

graph_raw <- graph_raw %>%
  mutate(expanded_names = case_when(
    cell == 'T4_naive' ~ 'CD4 Naive',
    cell == 'T8_naive' ~ 'CD8 Naive',
    cell == 'T4_em' ~ 'CD4 EM',
    cell == 'T4_reg' ~ 'CD4 Reg',
    cell == 'CytoT_GZMH+' ~ 'CD8 GZMH',
    cell == 'CytoT_GZMK+' ~ 'CD8 GZMK',
    cell == 'T_mait' ~ 'T MAIT',
    cell == 'NK_bright' ~ 'NK Bright',
    cell == 'NK_dim' ~ 'NK Dim',
    cell == 'B_naive' ~ 'B Naive',
    cell == 'B_mem' ~ 'B Mem',
    cell == 'B_plasma' ~ 'B Plasma',
    cell == 'B_atypical' ~ 'B Atypical',
    cell == 'Progen' ~ 'Progen'))

## Split into groups

raw_control <- graph_raw[graph_raw$group == "control_slope",]
raw_lupus <- graph_raw[graph_raw$group == "lupus_slope",]

# Plot -----

lupusg <- ggplot(data = raw_lupus, aes(x = reorder(expanded_names, slope), y = adjust_slope, fill = ifelse(pval<0.05, slope > 0, '6')))
lupusg <- lupusg + geom_bar(stat = "identity", position = 'dodge',  color = 'black')
lupusg <- lupusg + coord_flip() + labs(x="Cell type", y="Change in cell type percentage 
  per 10 years of age", 
                                       title="CLUES patients")
lupusg <- lupusg + my.theme
lupusg <- lupusg +  scale_fill_manual(name = "", values=c("TRUE" = "#e66101", "FALSE" = "#5e3c99", "6" = "white"), 
                                      labels = c("TRUE" = "Positive",
                                                 "FALSE" = "Negative",
                                                 "6" = substitute("Not significant (p"[adj]*" < 0.05)"))) + theme(panel.border = element_rect(color = "black", fill = NA), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),  
                                                                                                                  legend.text=element_text(size=10), legend.direction = "horizontal", legend.position = "none") + guides(fill=guide_legend(ncol=2)) + geom_hline(yintercept = 0) + scale_y_continuous(breaks = c(-1.5, 0.0, 1.5), limits = c(-1.5, 1.5))

controlg <- ggplot(data = raw_control, aes(x = reorder(raw_lupus$expanded_names, raw_lupus$slope), y = adjust_slope, fill = ifelse(pval<0.05, slope > 0, '6')))
controlg <- controlg + geom_bar(stat = "identity", position = 'dodge', color = 'black')
controlg <- controlg + coord_flip() + labs(x="Cell type", y="Change in cell type percentage 
  per 10 years of age", 
                                           title="Healthy individuals")
controlg <- controlg + my.theme
controlg <- controlg +  scale_fill_manual(name = "", values=c("TRUE" = "#ca0020", "FALSE" = "#0571b0", "6" = "white"), breaks = c("TRUE", "6", "FALSE"),
                                          labels = c("TRUE" = "Positive",
                                                     "FALSE" = "Negative",
                                                     "6" = substitute("Not significant (p"[adj]*" < 0.05)"), list(adj = "adj"))) + theme(panel.border = element_rect(color = "black", fill = NA), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), 
                                                                                                                                         legend.text=element_text(size=10), legend.direction = "horizontal", legend.position = "none") + guides(fill=guide_legend(ncol=2)) + geom_hline(yintercept = 0) + scale_y_continuous(breaks = c(-1.5, 0.0, 1.5), limits = c(-1.5, 1.5))
