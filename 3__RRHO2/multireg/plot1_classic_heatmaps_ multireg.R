library(RRHO2)

rm(list=ls())
setwd("/home/marinevernier/Documents/projets/female_cpid_multiregion/")

#### PATHS ####
rrho_obj.path <- "data/3__RRHO2/rrho_obj_multireg.Rdata"
output.path <- "graphs_results/3__RRHO2/multireg_classic/"

load(rrho_obj.path)
dir.create(output.path, recursive=TRUE, showWarnings = FALSE)
setwd(output.path)

#### rrho2 heatmap for all tp ####
list1 <- c("ACC", "Hb", "Ins", "Nac")
list2 <- list1

# homemade loop that iterates each combination of regions
for(tp in c(1,2,3)){
  print(paste0("Computing RRHO2 heatmaps of the TP", tp))
  # setwd(output.path)
  directory <- paste0("./TP", tp)
  dir.create(directory, showWarnings = FALSE)
  setwd(directory)
  i=2
  for(reg in list1[1:length(list1)-1]){
    for(reg2 in list2[i:length(list2)]){
      rrho <- get(paste0("rrho", reg, "vs", reg2, tp))
      png(file=paste0(reg, "vs", reg2, tp, ".png"), width = 750, height = 500, units = "px", res = 100)
      RRHO2_heatmap(rrho, main = paste0(reg, " vs ", reg2, " tp", tp))
      dev.off()
    }
    i=i+1
  }
}

##### TP1 ####
#dir.create("./TP1")
#setwd("./TP1")
#
## DRN
#png(file="DRNvsBLA1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsBLA1, main="DRN vs BLA tp1")
#dev.off()
#png(file="DRNvsHb1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsHb1, main="DRN vs Hb tp1")
#dev.off()
#png(file="DRNvsIns1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsIns1, main="DRN vs Ins tp1")
#dev.off()
#png(file="DRNvsNAc1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsNAc1, main="DRN vs NAc tp1")
#dev.off()
#png(file="DRNvsVTA1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsVTA1, main="DRN vs VTA tp1")
#dev.off()
#
## BLA
#png(file="BLAvsHb1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsHb1, main = "BLA vs Hb tp1")
#dev.off()
#png(file="BLAvsIns1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsIns1, main = "BLA vs Ins tp1")
#dev.off()
#png(file="BLAvsNAc1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsNAc1, main = "BLA vs NAc tp1")
#dev.off()
#png(file="BLAvsVTA1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsVTA1, main = "BLA vs VTA tp1")
#dev.off()
#
## Hb
#png(file="HbvsIns1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoHbvsIns1, main = "Hb vs Ins tp1")
#dev.off()
#png(file="HbvsNAc1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoHbvsNAc1, main = "Hb vs NAc tp1")
#dev.off()
#png(file="HbvsVTA1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoHbvsVTA1, main = "Hb vs VTA tp1")
#dev.off()
#
## Ins
#png(file="InsvsNAc1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoInsvsNAc1, main = "Ins vs NAc tp1")
#dev.off()
#png(file="InsvsVTA1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoInsvsVTA1, main = "Ins vs VTA tp1")
#dev.off()
#
## NAc
#png(file="NAcvsVTA1.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoNAcvsVTA1, main = "NAc vs VTA tp1")
#dev.off()
#
##### TP2 ####
#setwd(output.path)
#dir.create("./TP2")
#setwd("./TP2")
#
## DRN
#png(file="DRNvsBLA2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsBLA2, main="DRN vs BLA tp2")
#dev.off()
#png(file="DRNvsHb2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsHb2, main="DRN vs Hb tp2")
#dev.off()
#png(file="DRNvsIns2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsIns2, main="DRN vs Ins tp2")
#dev.off()
#png(file="DRNvsNAc2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsNAc2, main="DRN vs NAc tp2")
#dev.off()
#png(file="DRNvsVTA2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsVTA2, main="DRN vs VTA tp2")
#dev.off()
#
## BLA
#png(file="BLAvsHb2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsHb2, main = "BLA vs Hb tp2")
#dev.off()
#png(file="BLAvsIns2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsIns2, main = "BLA vs Ins tp2")
#dev.off()
#png(file="BLAvsNAc2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsNAc2, main = "BLA vs NAc tp2")
#dev.off()
#png(file="BLAvsVTA2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsVTA2, main = "BLA vs VTA tp2")
#dev.off()
#
## Hb
#png(file="HbvsIns2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoHbvsIns2, main = "Hb vs Ins tp2")
#dev.off()
#png(file="HbvsNAc2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoHbvsNAc2, main = "Hb vs NAc tp2")
#dev.off()
#png(file="HbvsVTA2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoHbvsVTA2, main = "Hb vs VTA tp2")
#dev.off()
#
## Ins
#png(file="InsvsNAc2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoInsvsNAc2, main = "Ins vs NAc tp2")
#dev.off()
#png(file="InsvsVTA2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoInsvsVTA2, main = "Ins vs VTA tp2")
#dev.off()
#
## NAc
#png(file="NAcvsVTA2.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoNAcvsVTA2, main = "NAc vs VTA tp2")
#dev.off()
#
##### TP3 ####
#setwd(output.path)
#dir.create("./TP3")
#setwd("./TP3")
#
## DRN
#png(file="DRNvsBLA3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsBLA3, main="DRN vs BLA tp3")
#dev.off()
#png(file="DRNvsHb3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsHb3, main="DRN vs Hb tp3")
#dev.off()
#png(file="DRNvsIns3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsIns3, main="DRN vs Ins tp3")
#dev.off()
#png(file="DRNvsNAc3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsNAc3, main="DRN vs NAc tp3")
#dev.off()
#png(file="DRNvsVTA3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoDRNvsVTA3, main="DRN vs VTA tp3")
#dev.off()
#
## BLA
#png(file="BLAvsHb3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsHb3, main = "BLA vs Hb tp3")
#dev.off()
#png(file="BLAvsIns3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsIns3, main = "BLA vs Ins tp3")
#dev.off()
#png(file="BLAvsNAc3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsNAc3, main = "BLA vs NAc tp3")
#dev.off()
#png(file="BLAvsVTA3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoBLAvsVTA3, main = "BLA vs VTA tp3")
#dev.off()
#
## Hb
#png(file="HbvsIns3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoHbvsIns3, main = "Hb vs Ins tp3")
#dev.off()
#png(file="HbvsNAc3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoHbvsNAc3, main = "Hb vs NAc tp3")
#dev.off()
#png(file="HbvsVTA3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoHbvsVTA3, main = "Hb vs VTA tp3")
#dev.off()
#
## Ins
#png(file="InsvsNAc3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoInsvsNAc3, main = "Ins vs NAc tp3")
#dev.off()
#png(file="InsvsVTA3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoInsvsVTA3, main = "Ins vs VTA tp3")
#dev.off()
#
## NAc
#png(file="NAcvsVTA3.png", width = 750, height = 500, units = "px", res = 100)
#RRHO2_heatmap(rrhoNAcvsVTA3, main = "NAc vs VTA tp3")
#dev.off()
#