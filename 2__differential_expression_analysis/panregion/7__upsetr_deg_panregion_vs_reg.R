#### UPSETR COMPARAISON DES DEG PANREGIONS = DEG PAR REGION POUR VOIR LES INTERACTIONS ENTRE LES DEG ####

rm(list=ls())

#### LIBRARIES ####
library(dplyr)
library(UpSetR)

#### PATHS #### 
project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

deglist.path <- "data/2__differential_expression_analysis/deglist.Rdata"


load(deglist.path)

timepoints <- c("tp1", "tp2", "tp3")

for (tp in timepoints) {
  
  #### INPUT ####
  panregion_deg.path <- paste0(
    "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_", tp, ".rds"
  )
  
  #### OUTPUT ####
  upset.path <- paste0(
    "graphs_results/panregion/deseq2_reg_group/upset_plots/upset_deg_panregion_vs_sr_", tp, ".png"
  )
  
  #### LOAD PANREGION ####
  panregion_deg <- readRDS(panregion_deg.path)
  
  panregion_up   <- panregion_deg %>% filter(diffexpressed == "UP")
  panregion_down <- panregion_deg %>% filter(diffexpressed == "DOWN")
  
  #### LOAD DEG REGIONS ####
  ACC <- get(paste0("deg_ACC_", tp))
  Hb  <- get(paste0("deg_Hb_", tp))
  Ins <- get(paste0("deg_Ins_", tp))
  Nac <- get(paste0("deg_Nac_", tp))
  
  #### SPLIT UP / DOWN ####
  ACC_up   <- ACC %>% filter(diffexpressed == "UP")
  ACC_down <- ACC %>% filter(diffexpressed == "DOWN")
  
  Hb_up    <- Hb %>% filter(diffexpressed == "UP")
  Hb_down  <- Hb %>% filter(diffexpressed == "DOWN")
  
  Ins_up   <- Ins %>% filter(diffexpressed == "UP")
  Ins_down <- Ins %>% filter(diffexpressed == "DOWN")
  
  Nac_up   <- Nac %>% filter(diffexpressed == "UP")
  Nac_down <- Nac %>% filter(diffexpressed == "DOWN")
  
  #### LIST INPUT ####
  listInput <- list(
    panregion_up   = panregion_up$label,
    panregion_down = panregion_down$label,
    ACC_up   = ACC_up$label,
    ACC_down = ACC_down$label,
    Hb_up    = Hb_up$label,
    Hb_down  = Hb_down$label,
    Ins_up   = Ins_up$label,
    Ins_down = Ins_down$label,
    Nac_up   = Nac_up$label,
    Nac_down = Nac_down$label
  )
  
  #### PLOT ####
  png(upset.path, width = 3000, height = 2400, res = 300)
  
  print(
    upset(
      fromList(listInput),
      mainbar.y.label = paste("DEGs Intersections -", tp),
      sets.x.label = "Number of DEGs",
      order.by = "freq",
      nsets = length(listInput)
    )
  )
  
  dev.off()
}
