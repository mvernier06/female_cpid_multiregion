library(RRHO2)
library(tidyverse)

rm(list=ls())
setwd("/home/marinevernier/Documents/projets/female_cpid_multiregion/")

#### PATHS ####
rrho_obj.path <- "data/3__RRHO2/rrho_obj_multitp.Rdata"
output.path <- "graphs_results/3__RRHO2/multitp_scaled/"
dir.create(output.path, showWarnings = FALSE, recursive = TRUE)

#### scaled rrho plots ####
load(rrho_obj.path)
setwd(output.path)


### Regions et comparaisons
regions <- c("ACC", "Hb", "Ins", "Nac")
comparisons <- c("1vs2", "2vs3", "1vs3")

### Calcul du max hypermat par région ----

envlist <- ls()
max_list <- list()

for (region in regions) {
  
  rrholist <- envlist[grepl(paste0("RRHO_", region, "_"), envlist)]
  maxreg <- 0
  
  for (obj in rrholist) {
    rrho <- get(obj)
    maxvalue <- max(rrho$hypermat, na.rm = TRUE)
    
    if (maxvalue > maxreg) {
      maxreg <- maxvalue
    }
  }
  
  max_list[[region]] <- maxreg
  cat("\nMax value of", region, ":", maxreg)
}

for (region in regions) {
  setwd("/home/marinevernier/Documents/projets/female_cpid_multiregion/")
  setwd(output.path)
  directory <- paste0("./", region )
  dir.create(directory, showWarnings = FALSE)
  setwd(directory) 
  for (comp in comparisons) {
    
    obj_name <- paste0("RRHO_", region, "_", comp)
    
    if (!exists(obj_name)) next
    
    rrho_obj <- get(obj_name)
    maxval <- max_list[[region]]
    
    file_name <- paste0(region, "_", comp, ".png")
    plot_title <- paste(region, comp)
    
    png(file = file_name,
        width = 750,
        height = 500,
        units = "px",
        res = 100)
    
    RRHO2_heatmap(rrho_obj,
                  main = plot_title,
                  maximum = maxval)
    
    dev.off()
  }
}
