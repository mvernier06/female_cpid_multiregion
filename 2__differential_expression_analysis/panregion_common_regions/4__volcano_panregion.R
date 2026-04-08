#### VOLCANO PLOT: DEG FOR ONE TP ####

rm(list=ls())

#### LIBRAIRIES ####
library(ggplot2)
library(tidyverse)
library(dplyr)

#### PATHS ####
project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

#### VOLCANO PLOT ####
for (tp in c("tp1", "tp2", "tp3")) {
  print(paste0("Executing volcano plot for: ", tp))
  deg_tp <- paste0("data/2__differential_expression_analysis/panregion_common_regions/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_", tp,".rds")
  volcano.path <- paste0("graphs_results/panregion_common_regions/deseq2_reg_group/volcano_panregion/volcano_panregion_deg_", tp,".png")
  
  deg_tp <- readRDS(deg_tp)
  
  volcano <- ggplot(deg_tp, aes(x = log2fc, y = -log10(pval), col=diffexpressed)) +
    geom_point() + theme_minimal() + 
    scale_color_manual(values=c("blue", "grey", "red"), name = "Differential expression") +
    #geom_vline(xintercept=c(-1, 1), col="grey", lty="dashed") +
    geom_hline(yintercept=-log10(0.05), col="grey", lty="dashed") +
    geom_hline(yintercept=-log10(0.05), col="grey", lty="dashed") +
    scale_y_continuous(limits = c(0, 8)) + 
    scale_x_continuous(limits = c(-0.5, 0.5)) +
    labs(title=paste0("Panregion DEG - ", tp), x ="log2FoldChange", y ="-log10 (p-value)")
  volcano
  
  #### SAVE ####
  ggsave(filename = volcano.path, plot = volcano, bg="white", width=1500, height=1000, units="px", scale=2)
  
  ## COUNTING VALUES ##
  counts_DEG <- table(deg_tp$diffexpressed)
  print(counts_DEG)
  
  # Filtrer les valeurs hors des limites définies
  outliers <- deg_tp %>%
    filter(log2fc < -0.5 | log2fc > 0.5 | -log10(pval) > 8 | is.na(log2fc) | is.na(pval))
  
  # Voir les lignes concernées
  print(outliers)
}
