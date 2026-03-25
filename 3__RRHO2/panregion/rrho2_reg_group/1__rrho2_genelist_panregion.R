#### RRHO2 : PREPARATION OF GENE LISTS (ONE PER TIMEPOINT) ####

#### LIBRARIES ####
library(tidyverse)
library(dplyr)
library(devtools)
library(RRHO2)

rm(list=ls())

#### PATH #### 
project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

## INPUT ##
annotated_counts.path <- "data/3__RRHO2/panregion/design_reg_group/prep_rrho2/annotated_counts_filtered_panregion_rg_alltp.rds"

## OUTPUT ##
output.path <- "data/3__RRHO2/panregion/design_reg_group/rrho2_genelist/"


#### REMOVE NA VALUES ####
annotated_counts_NA <- readRDS(annotated_counts.path)
annotated_counts <- annotated_counts_NA %>% dplyr::filter(!is.na(log2fc_tp1) & 
                                                            !is.na(pval_tp1) &
                                                            !is.na(log2fc_tp2) & 
                                                            !is.na(pval_tp2) &
                                                            !is.na(log2fc_tp3) & 
                                                            !is.na(pval_tp3))


### GENELIST OF TP1 ###
list1_DDE <- c(-log10(annotated_counts[annotated_counts$log2fc_tp1 < 0,]$pval_tp1) * (-1), 
               -log10(annotated_counts[annotated_counts$log2fc_tp1 > 0,]$pval_tp1))
gene_list1 <- data.frame(Genes=c(annotated_counts[annotated_counts$log2fc_tp1 < 0,]$MGI.symbol, 
                                 annotated_counts[annotated_counts$log2fc_tp1 > 0,]$MGI.symbol),
                         DDE = list1_DDE,
                         stringsAsFactors = FALSE)


### GENELIST OF TP2 ###
list2_DDE <- c(-log10(annotated_counts[annotated_counts$log2fc_tp2 < 0,]$pval_tp2) * (-1), 
               -log10(annotated_counts[annotated_counts$log2fc_tp2 > 0,]$pval_tp2))
gene_list2 <- data.frame(Genes=c(annotated_counts[annotated_counts$log2fc_tp2 < 0,]$MGI.symbol, 
                                 annotated_counts[annotated_counts$log2fc_tp2 > 0,]$MGI.symbol),
                         DDE = list2_DDE,
                         stringsAsFactors = FALSE)


### GENELIST OF TP3 ### 
list3_DDE <- c(-log10(annotated_counts[annotated_counts$log2fc_tp3 < 0,]$pval_tp3) * (-1), 
               -log10(annotated_counts[annotated_counts$log2fc_tp3 > 0,]$pval_tp3))
gene_list3 <- data.frame(Genes=c(annotated_counts[annotated_counts$log2fc_tp3 < 0,]$MGI.symbol, 
                                 annotated_counts[annotated_counts$log2fc_tp3 > 0,]$MGI.symbol),
                         DDE = list3_DDE,
                         stringsAsFactors = FALSE)


#### ACTIVATE METHOD = "FISCHER" FOR ODDS RATIO HEATMAP ####
RRHO_panregion_1vs2 <- RRHO2_initialize(gene_list1, gene_list2, 
                                           labels = c("TP1", "TP2"), 
                                           log10.ind = TRUE)

RRHO_panregion_2vs3 <- RRHO2_initialize(gene_list2, gene_list3, 
                                           labels = c("TP2", "TP3"), 
                                           log10.ind = TRUE)

RRHO_panregion_1vs3 <- RRHO2_initialize(gene_list1, gene_list3,
                                           labels = c("TP1", "TP3"), 
                                           log10.ind = TRUE)

#view(RRHO_allreg_rg_1vs2$hypermat) # to see the matrix

#### SAVE ####
rrho_obj <- ls()[grepl("RRHO",ls())]
setwd(output.path)
save(list=rrho_obj, file="rrho_obj_panregion_alltp.Rdata")
