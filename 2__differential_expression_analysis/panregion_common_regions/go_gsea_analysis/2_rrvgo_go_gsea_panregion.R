#### SIMPLIFICATION ET VISUALISATION DES TERMES GO ENRICHIS : PANREGION ####

rm(list=ls())

#### LIBRARIES ####
library(tidyverse)
library(rrvgo)
library(GOSemSim)
library(AnnotationDbi)
library(GO.db)
library(org.Mm.eg.db)

#### PATHS ####
project.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project.path)

## INPUT ## 
gsea_results <- "data/2__differential_expression_analysis/panregion_common_regions/go_gsea_analysis/fgsea_panregion_01_obj.Rdata"
go_results <- "data/2__differential_expression_analysis/panregion_common_regions/go_gsea_analysis/go_panregion_05_obj.Rdata"

## OUTPUT ##
# GSEA # 
reduced_gsea <- "data/2__differential_expression_analysis/panregion_common_regions/go_gsea_analysis/gsea/reduced_terms/"
treemap_gsea <- "graphs_results/panregion_common_regions/go_gsea_analysis/gsea/treemap"

# GO # 
reduced_go <- "data/2__differential_expression_analysis/panregion_common_regions/go_gsea_analysis/go/reduced_terms/"
treemap_go <- "graphs_results/panregion_common_regions/go_gsea_analysis/go/treemap"


#### SELECT CRITERIA #### 
list_ontology <- c("BP", "CC", "MF")
list_tp <- c(1,2,3)
list_reg <- "panregion_"

#### LOAD DATA ####
load(gsea_results)
load(go_results)

#### REDUCED TERMS FUNCTION #### 

go_rrvgo <- function(regions, timepoints, ontologies) {
  
  ## GO ## 
  for (reg in regions) { 
    message(paste0("Launch ", reg, "..."))
    for (tp in timepoints) {
      message(paste0("TP", tp, "..."))
      for (ont in ontologies) {
        message(paste0("Ontology: ", ont))
        go_name <- paste0("go_", reg, "05_", tp)
        go_res <- get(go_name, envir = .GlobalEnv)
        go_ont <- go_res@result$ID[go_res@result$ONTOLOGY == ont]
        go_ont_qval <- go_res@result$qvalue[go_res@result$ONTOLOGY == ont]
        
        ## get simMAtrix ## 
        
        if (!is.null(go_ont) && length(go_ont) > 1) {
          simMatrix <- calculateSimMatrix(go_ont, 
                                          orgdb = "org.Mm.eg.db", 
                                          ont = ont, 
                                          method = "Rel")
          if (nrow(simMatrix) > 1) { 
            
            scores <- setNames(-log10(go_ont_qval), go_ont)
            
            reducedTerms <- reduceSimMatrix(simMatrix, 
                                            scores, 
                                            threshold = 0.7,
                                            orgdb = "org.Mm.eg.db")
            
            
            ## SAVE ##
            # PLOT # 
            plot_folder <- file.path(treemap_go)
            if (!dir.exists(plot_folder)) {
              dir.create(plot_folder, recursive = TRUE)
            }
            
            output_file <- file.path(plot_folder, paste0("treemap_", reg, "05_tp", 
                                                         tp , "_", ont, ".png"))
            
            png(output_file, width = 8, height = 8, units = "in", res = 600)
            treemapPlot(reducedTerms)  # Utilise grid.draw pour afficher
            dev.off()
            
            # DATA #
            assign(paste("red_go_05", reg, tp, ont, sep = "_"), reducedTerms, envir = .GlobalEnv)
          }
          
          
        }
      }
    }
    message(paste0(reg, " done"))
  }
  
  message(paste0("GO simplification done. Start GSEA simplification ..."))
  
  ## GSEA ##
  for (reg in regions) { 
    message(paste0("Launch ", reg, "..."))
    for (tp in timepoints) {
      message(paste0("TP", tp, "..."))
      for (ont in ontologies) {
        message(paste0("Ontology: ", ont))
        
        gsea_name <- paste0("fgsea_", reg, "01_", tp)
        gsea_res <- get(gsea_name, envir = .GlobalEnv)
        gsea_ont <- gsea_res@result$ID[gsea_res@result$ONTOLOGY == ont]
        gsea_ont_qval <- gsea_res@result$qvalue[gsea_res@result$ONTOLOGY == ont]
        
        ## get simMAtrix ## 
        if (!is.null(gsea_ont) && length(gsea_ont) > 1) {
          
          simMatrix <- calculateSimMatrix(gsea_ont, 
                                          orgdb = "org.Mm.eg.db", 
                                          ont = ont, 
                                          method = "Rel")
          if (nrow(simMatrix) > 1) {
            
            scores <- setNames(-log10(gsea_ont_qval), gsea_ont)
            
            reducedTerms <- reduceSimMatrix(simMatrix, 
                                            scores, 
                                            threshold = 0.7,
                                            orgdb = "org.Mm.eg.db")
            
            ## SAVE ##
            # PLOT # 
            plot_folder <- file.path(treemap_gsea)
            if (!dir.exists(plot_folder)) {
              dir.create(plot_folder, recursive = TRUE)
            }
            
            output_file <- file.path(plot_folder, paste0("treemap_", reg, "01_tp", 
                                                         tp, "_", ont, ".png"))
            
            png(output_file, width = 8, height = 8, units = "in", res = 600)
            treemapPlot(reducedTerms)  # Utilise grid.draw pour afficher
            dev.off()
            
            # DATA #
            assign(paste("red_gsea_01", reg, tp, ont, sep = "_"), reducedTerms, envir = .GlobalEnv)
          }
        }
        
      }
    }
  }
  message(paste0("GSEA simplification done."))
}

#### RUN FUNCTION ####
go_rrvgo(list_reg, list_tp, list_ontology)

red_go_obj <- ls()[grepl("red_go_panregion_05_",ls())]
save(list=red_go_obj, file= paste0(reduced_go, "rrvgo_go_panregion_05.Rdata"))

red_gsea_obj <- ls()[grepl("red_gsea_panregion_01_",ls())]
save(list=red_gsea_obj, file= paste0(reduced_gsea, "rrvgo_gsea_panregion_01_.Rdata"))
