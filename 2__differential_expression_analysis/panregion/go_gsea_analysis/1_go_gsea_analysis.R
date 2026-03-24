#### ANALYSE PANREGION GO / GSEA ####
rm(list=ls())
#### LIBRARIES ####
dir.create("~/R/library", showWarnings = FALSE, recursive = TRUE)
.libPaths(c("~/R/library", .libPaths()))
Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always") 
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(fgsea)
library(msigdbr)

#### PATHS ####
project.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project.path)
## INPUT ## 
annotated_counts.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_alltp.rds"
deg_tp1 <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_tp1.rds"
deg_tp2 <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_tp2.rds"
deg_tp3 <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_tp3.rds"

## ANNOTATIONS DATABASE ##
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

## OUTPUT ##
output.path <- "data/2__differential_expression_analysis/panregion/go_gsea_analysis/"

#### ENRICHGO ####

## LOAD DATA ##
deg_panregion_1 <- readRDS(deg_tp1)
deg_panregion_2 <- readRDS(deg_tp2)
deg_panregion_3 <- readRDS(deg_tp3)
counts <- readRDS(annotated_counts.path)

## SELECT REGIONS AND TIMEPOINTS ##
tpList <- c("1", "2", "3")


## ENRICHGO ## 
for(tp in tpList) {
  message("Doing GO analysis of panregion for TP", tp )
  counts_temp <- counts %>% 
    dplyr::select(label, contains(paste0("tp", tp)) & contains("log2fc")) %>%
    na.omit
  
  deg_panreg_tp <- get(paste0("deg_panregion_", tp), envir = .GlobalEnv)
  
  go_panreg_tp <- enrichGO(gene = deg_panreg_tp$label, 
                           OrgDb = organism,
                           keyType = "SYMBOL", 
                           ont = "ALL", 
                           pvalueCutoff = 0.5,
                           qvalueCutoff = 0.5, 
                           pAdjustMethod = "BH", 
                           universe = counts_temp$label)
  
  assign(paste0("go_panregion_05_", tp), go_panreg_tp, envir = .GlobalEnv)
}

## SAVE RESULTS ##
rm(go_panreg_tp)
go_obj <- ls()[grepl("go_panregion_05_",ls())]
save(list=go_obj, file= paste0(output.path, "go_panregion_05_obj.Rdata"))

## COUNT RESULT ##
for(tp in tpList){
  go_obj <- paste0("go_panregion_05_", tp)
  print(paste0("Panregion tp", tp, ": ",nrow(get(go_obj))))
}


###  GSEA  ####
set.seed(123)
for(tp in tpList){
  print(tp)
  counts_temp <- counts %>%
    dplyr::select(label, 
                  matches(paste0("log2fc_tp", tp, "$")), 
                  matches(paste0("pval_tp", tp, "$"))) %>%
    na.omit()
  
  print(colnames(counts_temp))  # Vérifier les colonnes sélectionnées
  print(dim(counts_temp))  # Vérifier le nombre de lignes et colonnes
  
  if (ncol(counts_temp) == 3) {
    colnames(counts_temp) <- c("label", "log2fc", "pvalue")
  } else {
    print(" Problème de sélection des colonnes")
    print(colnames(counts_temp))
    next  # Passe à l'itération suivante si le problème persiste
  }
  
  # Création du score de ranking basé sur -log10(pvalue) * sign(log2fc)
  counts_temp <- counts_temp %>%
    mutate(rank_score = -log10(pvalue) * sign(log2fc))  # Création du score de ranking
  
  # Tri du tableau counts_temp par score de ranking
  counts_temp_sorted <- counts_temp %>%
    arrange(desc(rank_score))  # Trie les gènes en fonction du score de ranking
  
  # Création de la genelist avec le tri effectué
  genelist <- setNames(counts_temp_sorted$rank_score, counts_temp_sorted$label)
  
  
  # Exécution de la GSEA
  gse <- gseGO(geneList = genelist, 
               ont = "ALL", 
               keyType = "SYMBOL", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.1,
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "BH",
               by = "fgsea",
               seed = 123)
  
  assign(paste("fgsea_panregion_01", tp, sep = "_"), gse)
}

#### SAVE RESULTS ####
gse_obj <- ls()[grepl("fgsea_panregion_01_",ls())]
save(list=gse_obj, file= paste0(output.path, "fgsea_panregion_01_obj.Rdata"))


#### COUNTS RESULTS ####

for(tp in tpList){
  gse_obj <- paste("fgsea_panregion_01", tp, sep="_")
  print(paste0("Panregion tp", tp, ": ",nrow(get(gse_obj))))
}






