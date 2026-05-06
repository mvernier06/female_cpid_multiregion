.libPaths(c("~/R/library", .libPaths()))
Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always") 
library(clusterProfiler)
library(rrvgo)
library(org.Mm.eg.db)
library(dplyr)
library(openxlsx)
library(pheatmap)
library(tidyverse)

rm(list=ls())

project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

base.dir <- "data/2__differential_expression_analysis/alluvial_go_gsea/"
out_dir  <- "data/2__differential_expression_analysis/rrvgo_alluvial_pattern"
out_dir_graph <-  "graphs_results/2__differential_expression_analysis/rrvgo_alluvial_pattern"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_graph)

regions <- c("ACC", "Hb", "Ins", "Nac")
all_results <- list()

for (reg in regions) {
  
  file <- paste0(base.dir,"go_patterns_", reg, ".Rdata")
  go_patterns <- load(file)
  
  
  for (pat in go_patterns) {
    data <- get(pat)
    if (is.null(data) || nrow(data@result) == 0) next
    
    df <- data@result
    
    # garder termes significatifs
    df <- df %>% filter(p.adjust < 0.05)
    
    if (nrow(df) < 2) next
    
    # metadata depuis nom fichier
    pattern <- gsub(paste0("_",reg), "", pat) 
    pattern <- gsub("go_", "", pattern)
    
    df$pattern <- pattern
    df$name <- pat
    df$region <- reg 
    all_results[[pat]] <- df
  }
  
}

# combiner
go_all <- bind_rows(all_results)

# sauver tableau global

write.xlsx(go_all,
           file = file.path(out_dir, "GO_summary.xlsx"))

# boucle par comparaison

for (pattern_reg in unique(go_all$name)) {
  
  message("Processing ", pattern_reg)
  
  df <- go_all %>% filter(name == pattern_reg)
  
  ## créer dossier
  reg <- unique(df$reg)
  pattern <- unique(df$pattern)
  full_outdir <- file.path(out_dir_graph, reg, pattern)
  dir.create(full_outdir, recursive = TRUE, showWarnings = FALSE)
  
  ## -------------------------
  ## GO PAR ONTOLOGY
  ## -------------------------
  for (ont in c("BP","MF","CC")) {
    
    df_ont <- df %>% filter(ONTOLOGY == ont)
    
    if (nrow(df_ont) < 2) next
    
    terms  <- df_ont$ID
    scores <- setNames(-log10(df_ont$p.adjust), df_ont$ID)
    
    simMatrix <- calculateSimMatrix(
      terms,
      orgdb = org.Mm.eg.db,
      ont = ont,
      method = "Rel"
    )
    
    if (nrow(simMatrix) < 2) next
    
    reducedTerms <- reduceSimMatrix(
      simMatrix,
      scores,
      threshold = 0.7,
      orgdb = org.Mm.eg.db
    )
    
    ## -------------------------
    ## SAVE FILES
    ## -------------------------
    
    # TREEMAP
    png(file.path(full_outdir, paste0("treemap_", pattern_reg,"_", ont,  ".png")),
        width = 2000, height = 2000, res = 300)
    treemapPlot(reducedTerms)
    dev.off()
  }
}
