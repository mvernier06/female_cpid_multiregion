#### VOLCANO PLOT : DEG PANREGION VS 2 UNIQUE REGIONS #### 

rm(list=ls())

#### LIBRARIES ####
library(tidyverse)
library(ggplot2)

#### PATHS ####
project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

## INPUT ##
panregion_tp1.path <- "data/2__differential_expression_analysis/panregion_common_regions/deseq2/design_reg_group/annotation_deg/panregion_deg_tp1.rds"
panregion_tp2.path <- "data/2__differential_expression_analysis/panregion_common_regions/deseq2/design_reg_group/annotation_deg/panregion_deg_tp2.rds"
panregion_tp3.path <- "data/2__differential_expression_analysis/panregion_common_regions/deseq2/design_reg_group/annotation_deg/panregion_deg_tp3.rds"
deglist.path <- "data/2__differential_expression_analysis/deglist.Rdata"

## OUTPUT ##
output_dir <- "graphs_results/panregion_common_regions/deseq2_reg_group/volcano_panregion_vs_3reg/"
dir.create(output_dir, recursive = TRUE)

#### DATA FORMATTING #### 
panregion_tp1 <- readRDS(panregion_tp1.path)
panregion_tp2 <- readRDS(panregion_tp2.path)
panregion_tp3 <- readRDS(panregion_tp3.path)
load(deglist.path)

deg_all <- list(
  tp1 = list(Hb = deg_Hb_tp1, Ins = deg_Ins_tp1, Nac = deg_Nac_tp1),
  tp2 = list(Hb = deg_Hb_tp2, Ins = deg_Ins_tp2, Nac = deg_Nac_tp2),
  tp3 = list(Hb = deg_Hb_tp3, Ins = deg_Ins_tp3, Nac = deg_Nac_tp3)
)

regions <- c("Hb", "Ins", "Nac")
timepoints <- names(deg_all)


create_volcano <- function(df, reg1, reg2, reg3, tp) {
  ggplot(df, aes(x = log2fc, y = -log10(pval), col = region)) +
    geom_point() +
    theme_minimal() +
    scale_colour_viridis_d(alpha = 0.6) +
    geom_vline(xintercept = c(-1, 1), col = "grey", lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "grey", lty = "dashed") +
    labs(
      title = paste0("Panregion vs ", reg1, ", ", reg2, " and ", reg3," (", tp, ")"),
      x = "log2 FoldChange",
      y = "-log10 (p_value)",
      col = "Regions"
    ) +
    theme(plot.title = element_text(size = 15))
}

for (tp in timepoints) {
  
  deg_list_tp <- deg_all[[tp]]
  
  # récupérer panregion correspondant
  panregion_df <- get(paste0("panregion_", tp))
  
  # ajouter label
  panregion_df$region <- "panregion"
  
  reg1 <- "Hb"
  reg2 <- "Ins"
  reg3 <- "Nac"
  
  df1 <- deg_list_tp[[reg1]]
  df2 <- deg_list_tp[[reg2]]
  df3 <- deg_list_tp[[reg3]]
  
  # Ajouter labels
  df1$region <- reg1
  df2$region <- reg2
  df3$region <- reg3
  
  # Merge
  df <- bind_rows(panregion_df, df1, df2, df3)
  
  # Verifier qu'il n'y a pas de Na
  df <- df[!is.na(df$log2fc) & !is.na(df$pval), ]
  
  # Plot
  p <- create_volcano(df, reg1, reg2, reg3, tp)
  
  print(p)
  
  # Nom fichier propre
  filename <- paste0(
    output_dir,
    "volcano_panregion_vs_", reg1, "_", reg2, "_", reg3, "_",  tp, ".png"
  )
  
  # Save
  ggsave(
    filename = filename,
    plot = p,
    bg = "white",
    width = 1500,
    height = 1000,
    units = "px",
    scale = 2
  )
  
}
