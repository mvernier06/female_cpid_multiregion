#### GO ENRICHMENT (ENRICHGO OR GSEA) ####

dir.create("~/R/library", showWarnings = FALSE, recursive = TRUE)
.libPaths(c("~/R/library", .libPaths()))

Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always") 

#### LIBRARIES ####
library(clusterProfiler)
library(tidyverse)

#### PATH ####
## INPUT ##
annotated_counts.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts_filtered.rds"
deglist.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/deglist.Rdata"

## ANNOTATIONS DATABASE ##
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

## OUTPUT ##
output.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/"


#### ENRICHGO ####

## LOAD DATA ##
load(deglist.path)
counts <- readRDS(annotated_counts.path)

## SELECT REGIONS AND TIMEPOINTS ##
regionList <-c("ACC", "Ins", "Hb", "Nac")
tpList <- c("1", "2", "3")

## ENRICHGO FUNCTION ##
for(reg in regionList){
  for(tp in tpList){
    message("Doing GO analysis of ", reg, " for TP", tp)
    counts_temp <- counts %>% 
      dplyr::select(MGI.symbol, contains(reg) & contains(tp) & !contains("padj")) %>% 
      na.omit
    
    deg_regtp <- get(paste0("deg_", reg, "_tp", tp))
    
    go_regtp <- enrichGO(gene = deg_regtp$label, 
                         OrgDb = organism, 
                         keyType = "SYMBOL", 
                         ont ="ALL",
                         qvalueCutoff = 1,
                         pAdjustMethod = "BH",
                         universe = counts_temp$MGI.symbol)
    
    assign(paste("go", reg, tp, sep="_"), go_regtp, envir = .GlobalEnv)
  }
}

## SAVE RESULTS ##
rm(go_regtp)
go_obj <- ls()[grepl("go_",ls())]
save(list=go_obj, file= paste0(output.path, "go_obj.Rdata"))

## COUNT RESULT ##
for(reg in regionList){
  for(tp in tpList){
    go_obj <- paste("go", reg, tp, sep="_")
    print(paste0(reg, " tp", tp, ": ",nrow(get(go_obj))))
  }
}


#### GSEA ####
set.seed(123)
for(reg in regionList){
  for(tp in tpList){
    message("Doing GSEA analysis of ", reg, " for TP", tp)
    
    counts_temp <- counts %>%
      dplyr::select(MGI.symbol, 
                    matches(paste0("^", reg, "_log2fc_tp", tp, "$")), 
                    matches(paste0("^", reg, "_pval_tp", tp, "$"))) %>%
      na.omit()
    
    #print(colnames(counts_temp))  # Vérifier les colonnes sélectionnées
    #print(dim(counts_temp))  # Vérifier le nombre de lignes et colonnes
    
    if (ncol(counts_temp) == 3) {
      colnames(counts_temp) <- c("MGI.symbol", "log2fc", "pvalue")
    } else {
      print(" Problème de sélection des colonnes")
      print(colnames(counts_temp))
      next
    }
    
    # Ranking
    counts_temp <- counts_temp %>%
      mutate(rank_score = -log10(pvalue) * sign(log2fc))  # score for ranking
    
    # Sort according to ranking score 
    counts_temp_sorted <- counts_temp %>%
      arrange(desc(rank_score))  # sorting 
    
    # gene list of the sorted gene set 
    genelist <- setNames(counts_temp_sorted$rank_score, counts_temp_sorted$MGI.symbol)
    
    
    # GSEA
    gse <- gseGO(geneList = genelist, 
                 ont = "ALL", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 eps = 1e-300,
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "BH", 
                 by = "fgsea",
                 seed = 123)
    
    assign(paste("fgsea", reg, tp, sep = "_"), gse)
  }
}

## SAVE ##
rm(gse)
gse_obj <- ls()[grepl("fgsea_",ls())]
save(list=gse_obj, file= paste0(output.path, "fgsea_obj.Rdata"))


## COUNT RESULT ##
for(reg in regionList){
  for(tp in tpList){
    gse_obj <- paste("fgsea", reg, tp, sep="_")
    print(paste0(reg, " tp", tp, ": ",nrow(get(gse_obj))))
  }
}

