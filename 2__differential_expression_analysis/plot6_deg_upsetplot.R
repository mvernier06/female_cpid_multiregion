library(tidyverse)
library(UpSetR)

rm(list=ls())

#### PATHS ####
deglist.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/deglist.Rdata"
plot.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/deg_upsetplot"

dir.create(plot.path)
setwd(plot.path)
load(deglist.path)

regionlist <- c("ACC","Hb", "Ins", "Nac")

for(reg in regionlist){
  regtp1 <- get(paste0("deg_", reg, "_tp1"))
  regtp2 <- get(paste0("deg_", reg, "_tp2"))
  regtp3 <- get(paste0("deg_", reg, "_tp3"))
  
  listInput <- list(tp1 = regtp1$label,
                    tp2 = regtp2$label,
                    tp3 = regtp3$label)
  
  png(file=paste0("deg_upsetplot_", reg, ".png"),
      width = 1000, height = 700, 
      units = "px", res = 120)
  
  print(upset(fromList(listInput), 
              mainbar.y.label = paste0("DEGs Intersections in the ", reg), 
              sets.x.label = "Number of DEGs per timepoint", order.by="freq"))
  
  dev.off() 
}

for(reg in regionlist){
  regtp1_up <- get(paste0("deg_", reg, "_tp1")) %>%
    filter(diffexpressed == "UP")
  regtp1_down <- get(paste0("deg_", reg, "_tp1")) %>%
    filter(diffexpressed == "DOWN")
  regtp2_up <- get(paste0("deg_", reg, "_tp2")) %>%
    filter(diffexpressed == "UP")
  regtp2_down <- get(paste0("deg_", reg, "_tp2")) %>%
    filter(diffexpressed == "DOWN")
  regtp3_up <- get(paste0("deg_", reg, "_tp3")) %>%
    filter(diffexpressed == "UP")
  regtp3_down <- get(paste0("deg_", reg, "_tp3")) %>%
    filter(diffexpressed == "DOWN")
  
  listInput <- list(tp1_UP = regtp1_up$label,
                    tp1_DOWN = regtp1_down$label,
                    tp2_UP = regtp2_up$label,
                    tp2_DOWN = regtp2_down$label,
                    tp3_UP = regtp3_up$label,
                    tp3_DOWN = regtp3_down$label)
  
  png(file=paste0("deg_upsetplot_", reg, "_updown.png"),
      width = 1000, height = 700, 
      units = "px", res = 120)
  
  print(upset(fromList(listInput), 
              mainbar.y.label = paste0("DEGs Intersections in the ", reg), 
              sets.x.label = "Number of DEGs per timepoint", order.by="freq", nsets = 6))
  
  dev.off() 
}
