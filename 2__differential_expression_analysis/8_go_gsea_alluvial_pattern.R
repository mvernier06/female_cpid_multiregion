.libPaths(c("~/R/library", .libPaths()))

Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always") 

#### LIBRARIES ####
library(clusterProfiler)
library(tidyverse)

rm(list=ls())

#### PATH ####
project.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project.path)
## INPUT ##
annotated_counts.path <- "data/2__differential_expression_analysis/annotated_counts_filtered.rds"
pattern.path <- "data/2__differential_expression_analysis/"

## ANNOTATIONS DATABASE ##
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

## OUTPUT ##
output.path <- "data/2__differential_expression_analysis/alluvial_go_gsea/"
dir.create(output.path)


## LOAD DATA ##

counts <- readRDS(annotated_counts.path)

## SELECT REGIONS AND TIMEPOINTS ##
regionList <-c("ACC", "Ins", "Hb", "Nac")
for (reg in regionList ){
  pattern_reg.path <- paste0(pattern.path, "alluvial_patterns_", reg, ".Rdata")
  load(pattern_reg.path)
  print(paste0("Processing : ", reg))
  
  for (pattern in unique(df_new$diffexpressed_alltp)) {
    print(paste0("Processing : ", pattern))
    df <- df_new %>% filter(diffexpressed_alltp == pattern)
    
    ## ENRICHGO FUNCTION ##
    go_pattern <- enrichGO(gene = df$MGI.symbol, 
                           OrgDb = organism, 
                           keyType = "SYMBOL", 
                           ont ="ALL",
                           qvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           universe = counts$MGI.symbol)
    
    if (!is.null(go_pattern)) {
      print(length(go_pattern@result$Description))
    }
    assign(paste("go", pattern, reg, sep="_"), go_pattern, envir = .GlobalEnv)
  }
  
  go_obj <- ls()[grepl(reg,ls())]
  save(list=go_obj, file= paste0(output.path, "go_patterns_", reg,".Rdata"))
  
}

## SAVE RESULTS ##
rm(go_pattern)
go_obj <- ls()[grepl("go_",ls())]
save(list=go_obj, file= paste0(output.path, "go_patterns_all.Rdata"))



