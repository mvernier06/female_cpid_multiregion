library(tidyverse)
library(UpSetR)

rm(list=ls())

project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

#### PATHS ####
deglist.path <- "data/2__differential_expression_analysis/deglist.Rdata"
plot.path <- "graphs_results/2__differential_expression_analysis/deg_upsetplot/multitp"

dir.create(plot.path)
load(deglist.path)
setwd(plot.path)

regionlist <- c("ACC","Hb", "Ins", "Nac")

for(reg in regionlist){
  regtp1 <- get(paste0("deg_", reg, "_tp1"))
  regtp2 <- get(paste0("deg_", reg, "_tp2"))
  regtp3 <- get(paste0("deg_", reg, "_tp3"))
  # print("TP1 : ")
  # print(length(regtp1$label))
  # print(length(unique(regtp1$label)))
  # print("TP2 : ")
  # print(length(regtp2$label))
  # print(length(unique(regtp2$label)))
  # print("TP3 : ")
  # print(length(regtp3$label))
  # print(length(unique(regtp3$label)))
  
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

### upset between region for all tp ###
setwd(project_path)
plot.path <- "graphs_results/2__differential_expression_analysis/deg_upsetplot/multireg"
dir.create(plot.path)
setwd(plot.path)

for(tp in c("tp1","tp2","tp3")){
  
  listInput <- list()
  print(tp)
  
  for(reg in regionlist){
    print(reg)
    
    deg <- get(paste0("deg_", reg, "_", tp))
    print(length(deg$label))
    print(length(unique(deg$label)))
    
    listInput[[reg]] <- deg$label
  }
  
  png(file=paste0("deg_upsetplot_", tp, "_regions.png"),
      width = 1000, height = 700, 
      units = "px", res = 120)
  
  print(
    upset(fromList(listInput),
          mainbar.y.label = paste0("DEGs intersections between regions (", tp, ")"),
          sets.x.label = "Number of DEGs per region",
          order.by = "freq")
  )
  
  dev.off()
}

for(tp in c("tp1","tp2","tp3")){
  
  listInput <- list()
  
  for(reg in regionlist){
    
    deg <- get(paste0("deg_", reg, "_", tp))
    
    listInput[[paste0(reg,"_UP")]] <- deg %>%
      filter(diffexpressed == "UP") %>%
      pull(label)
    
    listInput[[paste0(reg,"_DOWN")]] <- deg %>%
      filter(diffexpressed == "DOWN") %>%
      pull(label)
  }
  
  png(file=paste0("deg_upsetplot_", tp, "_regions_updown.png"),
      width = 1000, height = 700,
      units = "px", res = 120)
  
  print(
    upset(fromList(listInput),
          mainbar.y.label = paste0("DEGs intersections between regions (", tp, ")"),
          sets.x.label = "Number of DEGs",
          order.by = "freq")
  )
  
  dev.off()
}
