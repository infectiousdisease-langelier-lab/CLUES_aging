# Here I create a dotplot based on the output of pseudobulk_analysis.R

# Are main idea comments

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
library(ggeffects)
library(purrr)
library(ggpubr)
library(svglite)
library("readxl")

# Read in data -----

third_ovr <- read_excel("./pathway_dotplot.xlsx")

# Separate data by CLUES and RS -----

CLUES <- subset(third_ovr, Group == 'CLUES')

RS <- subset(third_ovr, Group != 'CLUES')


CLUES$CT <- factor(CLUES$CT, levels = unique(CLUES$CT))

RS$CT <- factor(RS$CT, levels = unique(RS$CT))

# Create the dot plot

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm"),
    panel.grid = element_blank()
  )

## Define the desired order for Pathway and Cell type 

ct_order <- c("NK", "T4_naive", "T4_EM", "T4_reg", "T8_naive", "B_naive", "B_EM", "ncM", "cM", "cDC")
Pathway_Order <- c("Alpha_beta", "Gamma", "Cytokine_Signaling", "IL-1_Signaling")
PO <- rev(Pathway_Order)

## Convert Pathway and CT to factors with the desired order
CLUES$Pathway <- factor(CLUES$Pathway, levels = PO)
CLUES$CT <- factor(CLUES$CT, levels = ct_order)
RS$Pathway <- factor(RS$Pathway, levels = PO)
RS$CT <- factor(RS$CT, levels = ct_order)

pathway_CLUES <- ggplot(CLUES, aes(x = CT, y = Pathway, fill = NES, color = Padj < 0.05)) +
  geom_point(shape = 21, size = 6, aes(stroke = 2)) +
  scale_fill_gradientn(
    colors = c("#5e3c99", "white", "#e66101"),
    limits = c(-max(abs(CLUES$NES)), max(abs(CLUES$NES))),
    breaks = c(-1.3, 0, 1.3),
    guide = "colorbar"
  ) + 
  scale_color_manual(values = c("transparent", "black")) +
  ggtitle("CLUES") + my.theme +  
  geom_vline(xintercept=seq(1.5, length(unique(CLUES$CT)), 1), 
             lwd=0.5, colour="gray") + 
  geom_hline(yintercept=seq(1.5, length(unique(CLUES$Pathway)), 1), 
             lwd=0.5, colour="gray") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 2),
        plot.margin = margin(1, 1, 1, 1, "cm")) + coord_flip()                                                


pathway_RS <- ggplot(RS, aes(x = CT, y = Pathway, fill = NES, color = Padj < 0.05)) +
  geom_point(shape = 21, size = 6, aes(stroke = 2)) +
  scale_fill_gradientn(
    colors = c("#0571b0", "white", "#ca0020"),
    limits = c(-max(abs(RS$NES)), max(abs(RS$NES))),
    breaks = c(-1.3, 0, 1.3),
    guide = "colorbar"
  ) + 
  scale_color_manual(values = c("transparent", "black")) +
  ggtitle("Healthy controls") + my.theme + 
  geom_vline(xintercept=seq(1.5, length(unique(CLUES$CT)), 1), 
             lwd=0.5, colour="gray") + 
  geom_hline(yintercept=seq(1.5, length(unique(CLUES$Pathway)), 1), 
             lwd=0.5, colour="gray") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 2),
        plot.margin = margin(1, 1, 1, 1, "cm")) + coord_flip()



