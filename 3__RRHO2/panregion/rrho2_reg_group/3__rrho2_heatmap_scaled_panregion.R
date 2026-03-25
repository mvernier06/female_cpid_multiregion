#### RRHO2: HEATMAPS (SCALED) ####

#### LIBRARIES ####
library(RRHO2)
library(tidyverse)

rm(list=ls())

#### PATHS #### 
project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

## INPUT ##
rrho_obj.path <- "data/3__RRHO2/panregion/design_reg_group/rrho2_genelist/rrho_obj_panregion_alltp.Rdata"

## OUTPUT ##
output.path <- "graphs_results/panregion/rrho2/design_reg_group/rrho2_scaled/"

#### RRHO2 HEATMAPS (SCALED) #### 
load(rrho_obj.path)
setwd(output.path)

#### MAX VALUE FOR SCALING ####
zmax <- max(
  RRHO_panregion_1vs2$hypermat,
  RRHO_panregion_2vs3$hypermat,
  RRHO_panregion_1vs3$hypermat,
  na.rm = TRUE
)

### HEATMAPS ### 
## TP1 VS TP2 ##
png(file="panregion_1vs2_scaled.png", width = 750, height = 500, units = "px", res = 100)
RRHO2_heatmap(RRHO_panregion_1vs2, main="panregion 1vs2 (design ~reg+group)", maximum = zmax)
dev.off()

## TP2 VS TP3 ##
png(file="panregion_2vs3_scaled.png", width = 750, height = 500, units = "px", res = 100)
RRHO2_heatmap(RRHO_panregion_2vs3, main="panregion 2vs3 (design ~reg+group)", maximum = zmax)
dev.off()

## TP1 VS TP3 ##
png(file="panregion_1vs3_scaled.png", width = 750, height = 500, units = "px", res = 100)
RRHO2_heatmap(RRHO_panregion_1vs3, main="panregion 1vs3 (design ~reg+group)", maximum = zmax)
dev.off()
