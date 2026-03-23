#### PLOT CLUSTER PROFILER #### 
rm(list=ls())
#### LIBRARIES ####
library(dplyr)
library(tidyverse)
library(pheatmap)

#### PATHS ####
project.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project.path)
## INPUT ##
gsea_tp1 <- "data/2__differential_expression_analysis/panregion/go/gse/gse_panregion_tp1.rds"
gsea_tp2 <- "data/2__differential_expression_analysis/panregion/go/gse/gse_panregion_tp2.rds"
gsea_tp3 <- "data/2__differential_expression_analysis/panregion/go/gse/gse_panregion_tp3.rds"

## OUTPUT ##
gsea_alltp.path <- "data/2__differential_expression_analysis/panregion/go/gse/gsea_panregion_alltp.rds"


#### COMBINE ALL GSEA RESULTS FOR ALL TP  ####
gsea_tp1 <- as.data.frame(readRDS(gsea_tp1))
gsea_tp2 <- as.data.frame(readRDS(gsea_tp2))
gsea_tp3 <- as.data.frame(readRDS(gsea_tp3))

## TP1 ##
gsea_tp1 <-  gsea_tp1 %>% dplyr::select(Description, enrichmentScore, pvalue, p.adjust, qvalue)
colnames(gsea_tp1) <- c("Description", "enrichmentScore_tp1","pvalue_tp1",
                        "p.adjust_tp1", "qvalue_tp1")

## TP2 ##
gsea_tp2 <-  gsea_tp2 %>% dplyr::select(Description, enrichmentScore, pvalue, p.adjust, qvalue)
colnames(gsea_tp2) <- c("Description", "enrichmentScore_tp2","pvalue_tp2",
                        "p.adjust_tp2", "qvalue_tp2")

## TP3 ##
gsea_tp3 <-  gsea_tp3 %>% dplyr::select(Description, enrichmentScore, pvalue, p.adjust, qvalue)
colnames(gsea_tp3) <- c("Description", "enrichmentScore_tp3","pvalue_tp3",
                        "p.adjust_tp3", "qvalue_tp3")

# join df by "Description"
gsea_alltp <- full_join(gsea_tp1,gsea_tp2, by="Description")
gsea_alltp <- full_join(gsea_alltp,gsea_tp3, by = "Description" )

#### SAVE ####
saveRDS(gsea_alltp, gsea_alltp.path)



#### HEATMAP #### 

# heatmap function
top10=TRUE
go_heatmap <- function(query, top10 = FALSE){
  # format data
  test <- gsea_alltp %>% 
    filter(grepl(query, Description)) %>%
    dplyr::select(Description, contains("enrichmentScore"))
  if(top10==TRUE) test <- head(test, 10)
  test[,-1] <- lapply(test[,-1], function(x) round(x, digits=2))
  test2 <- test %>% column_to_rownames("Description")
  
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
                labels_col = c("TP1", "TP2", "TP3"),
                angle_col = 315,
                annotation_colors = "test")
  print(p)
}

#### QUERIES ####
go_heatmap("myelin")
# go_heatmap("calcium") # pas de calcium
go_heatmap("oligo")
# go_heatmap("neuron")
go_heatmap("axon")
go_heatmap("micro")
go_heatmap("", top10=TRUE)
