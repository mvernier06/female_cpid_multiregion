### Problème pas de zScore avec enrichGO ?? ###

library(tidyverse)
library(UpSetR)

rm(list=ls())

#### PATH ####
go_obj.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/go_obj.Rdata"
plot.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/go_gsea_analysis/go/"
dir.create(plot.path, showWarnings = FALSE)



#### LOAD SAVED DATA ####
load(go_obj.path)



# intersect of GO pathways per tp for each region #
setwd(plot.path)
listInput <- list(tp1 = go_ACC_1$Description,
                  tp2 = go_ACC_2$Description,
                  tp3 = go_ACC_3$Description)
png(filename = "upset_go_acc.png", height = 1250, width = 1250, units = "px", res = 225)
upset(fromList(listInput), 
      mainbar.y.label = paste0("GO BP Intersections in the ACC"), 
      sets.x.label = "Number of BP per timepoint", order.by="freq")
dev.off()

listInput <- list(tp1 = go_Hb_1$Description,
                  tp2 = go_Hb_2$Description,
                  tp3 = go_Hb_3$Description)
png(filename = "upset_go_hb.png", height = 1250, width = 1250, units = "px", res = 225)
upset(fromList(listInput), 
      mainbar.y.label = paste0("GO BP Intersections in the Hb"), 
      sets.x.label = "Number of BP per timepoint", order.by="freq")
dev.off()

# Ne fonctionne pas car pas de GO enrichis 
# listInput <- list(tp1 = go_Ins_1$Description,
#                   tp2 = go_Ins_2$Description,
#                   tp3 = go_Ins_3$Description)
# png(filename = "upset_go_ins.png", height = 1250, width = 1250, units = "px", res = 225)
# upset(fromList(listInput),
#       mainbar.y.label = paste0("GO BP Intersections in the Ins"),
#       sets.x.label = "Number of BP per timepoint", order.by="freq")
# dev.off()

listInput <- list(tp1 = go_Nac_1$Description,
                  tp2 = go_Nac_2$Description,
                  tp3 = go_Nac_3$Description)
png(filename = "upset_go_nac.png", height = 1250, width = 1250, units = "px", res = 225)
upset(fromList(listInput), 
      mainbar.y.label = paste0("GO BP Intersections in the Nac"), 
      sets.x.label = "Number of BP per timepoint", order.by="freq")
dev.off()



#### GO TABLE ####
# big df with all results
acc1 <- data.frame(description=go_ACC_1$Description, acc1_pval=go_ACC_1$pvalue, acc1_zscore=go_ACC_1$q)
acc2 <- data.frame(description=go_ACC_2$Description, acc2_pval=go_ACC_2$pvalue, acc2_zscore=go_ACC_2$zScore)
acc3 <- data.frame(description=go_ACC_3$Description, acc3_pval=go_ACC_3$pvalue, acc3_zscore=go_ACC_3$zScore)

hb1 <- data.frame(description=go_Hb_1$Description, hb1_pval=go_Hb_1$pvalue, hb1_zscore=go_Hb_1$zScore)
hb2 <- data.frame(description=go_Hb_2$Description, hb2_pval=go_Hb_2$pvalue, hb2_zscore=go_Hb_2$zScore)
hb3 <- data.frame(description=go_Hb_3$Description, hb3_pval=go_Hb_3$pvalue, hb3_zscore=go_Hb_3$zScore)

ins1 <- data.frame(description=go_Ins_1$Description, ins1_pval=go_Ins_1$pvalue, ins1_zscore=go_Ins_1$zScore)
ins2 <- data.frame(description=go_Ins_2$Description, ins2_pval=go_Ins_2$pvalue, ins2_zscore=go_Ins_2$zScore)
ins3 <- data.frame(description=go_Ins_3$Description, ins3_pval=go_Ins_3$pvalue, ins3_zscore=go_Ins_3$zScore)

nac1 <- data.frame(description=go_Nac_1$Description, nac1_pval=go_Nac_1$pvalue, nac1_zscore=go_Nac_1$zScore)
nac2 <- data.frame(description=go_Nac_2$Description, nac2_pval=go_Nac_2$pvalue, nac2_zscore=go_Nac_2$zScore)
nac3 <- data.frame(description=go_Nac_3$Description, nac3_pval=go_Nac_3$pvalue, nac3_zscore=go_Nac_3$zScore)


go_allreg <- list(acc1, acc2, acc3,
                  hb1, hb2, hb3,
                  ins1, ins2, ins3,
                  nac1, nac2, nac3) %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="description"), .) %>% arrange(description)
output.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/go/"
dir.create(output.path)
setwd(output.path)
write_rds(go_allreg, "go_allreg.rds")

library(pheatmap)
library(tidyverse)
# heatmap function
top10=TRUE
go_heatmap <- function(query, top10 = FALSE){
  # format data
  test <- go_allreg %>% 
    filter(grepl(query, description)) %>%
    dplyr::select(description, contains("zscore"))
  if(top10==TRUE) test <- head(test, 10)
  test[,-1] <- lapply(test[,-1], function(x) round(x, digits=2))
  test2 <- test %>% column_to_rownames("description")
  
  # max absolute value to scale the color scale
  maxval <- test2 %>% 
    apply(2, function(x) replace_na(x, 0)) %>% 
    as.matrix %>% abs %>% max
  
  # make a table for all pathway found 
  p <- pheatmap(test2,
                display_numbers = TRUE,
                fontsize = 8,
                breaks = seq(-maxval, maxval, length.out = 100), # center the color scale on zero
                cluster_cols = FALSE,
                cluster_rows = FALSE,
                labels_col = c(paste("acc", c(1,2,3), sep = "_"),
                               paste("hb", c(1,2,3), sep = "_"),
                               paste("ins", c(1,2,3), sep = "_"),
                               paste("nac", c(1,2,3), sep = "_")),
                angle_col = 315,
                annotation_colors = "test")
  print(p)
}

# call to the heatmap function
go_heatmap("myelin")
go_heatmap("oligo")
go_heatmap("neuron")
go_heatmap("", top10=TRUE)

#### data pour ENORA ####
go_allreg_enora <- go_allreg %>% dplyr::select(contains("DRN"))

view(go_DRN_2)