library(tidyverse)
library(UpSetR)

rm(list=ls())

#### PATH ####
gse_obj.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/gse/gse_obj.Rdata"
plot.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/go_gsea_analysis/gsea/"
# dir.create(plot.path, showWarnings = FALSE)



# upsetR résultats GSEA
setwd(plot.path)
load(gse_obj.path)

listInput <- list(tp1_up = gse_ACC_1$Description[which(sign(gse_ACC_1$NES)>0)],
                  tp1_down = gse_ACC_1$Description[which(sign(gse_ACC_1$NES)<0)],
                  tp2_up = gse_ACC_2$Description[which(sign(gse_ACC_2$NES)>0)],
                  tp2_down = gse_ACC_2$Description[which(sign(gse_ACC_2$NES)<0)],
                  tp3_up = gse_ACC_3$Description[which(sign(gse_ACC_3$NES)>0)],
                  tp3_down = gse_ACC_3$Description[which(sign(gse_ACC_3$NES)<0)])
png(filename = "upset_gse_acc.png", height = 1250, width = 1250, units = "px", res = 225)
upset(fromList(listInput), 
      mainbar.y.label = paste0("GO BP Intersections in the ACC"), 
      sets.x.label = "Number of BP per timepoint", order.by="freq")
dev.off()

listInput <- list(tp1_up = gse_Hb_1$Description[which(sign(gse_Hb_1$NES)>0)],
                  tp1_down = gse_Hb_1$Description[which(sign(gse_Hb_1$NES)<0)],
                  tp2_up = gse_Hb_2$Description[which(sign(gse_Hb_2$NES)>0)],
                  tp2_down = gse_Hb_2$Description[which(sign(gse_Hb_2$NES)<0)],
                  tp3_up = gse_Hb_3$Description[which(sign(gse_Hb_3$NES)>0)],
                  tp3_down = gse_Hb_3$Description[which(sign(gse_Hb_3$NES)<0)])
png(filename = "upset_gse_hb.png", height = 1250, width = 1250, units = "px", res = 225)
upset(fromList(listInput), 
      mainbar.y.label = paste0("GO BP Intersections in the Hb"), 
      sets.x.label = "Number of BP per timepoint", order.by="freq")
dev.off()

listInput <- list(tp1_up = gse_Ins_1$Description[which(sign(gse_Ins_1$NES)>0)],
                  tp1_down = gse_Ins_1$Description[which(sign(gse_Ins_1$NES)<0)],
                  tp2_up = gse_Ins_2$Description[which(sign(gse_Ins_2$NES)>0)],
                  tp2_down = gse_Ins_2$Description[which(sign(gse_Ins_2$NES)<0)],
                  tp3_up = gse_Ins_3$Description[which(sign(gse_Ins_3$NES)>0)],
                  tp3_down = gse_Ins_3$Description[which(sign(gse_Ins_3$NES)<0)])
png(filename = "upset_gse_ins.png", height = 1250, width = 1250, units = "px", res = 225)
upset(fromList(listInput), 
      mainbar.y.label = paste0("GO BP Intersections in the Ins"), 
      sets.x.label = "Number of BP per timepoint", order.by="freq")
dev.off()



listInput <- list(tp1_up = gse_Nac_1$Description[which(sign(gse_Nac_1$NES)>0)],
                  tp1_down = gse_Nac_1$Description[which(sign(gse_Nac_1$NES)<0)],
                  tp2_up = gse_Nac_2$Description[which(sign(gse_Nac_2$NES)>0)],
                  tp2_down = gse_Nac_2$Description[which(sign(gse_Nac_2$NES)<0)],
                  tp3_up = gse_Nac_3$Description[which(sign(gse_Nac_3$NES)>0)],
                  tp3_down = gse_Nac_3$Description[which(sign(gse_Nac_3$NES)<0)])
png(filename = "upset_gse_nac.png", height = 1250, width = 1250, units = "px", res = 225)
upset(fromList(listInput), 
      mainbar.y.label = paste0("GO BP Intersections in the NAc"), 
      sets.x.label = "Number of BP per timepoint", order.by="freq")
dev.off()
