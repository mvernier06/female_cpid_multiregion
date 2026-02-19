#### SIMPLIFICATION ET VISUALISATION DES TERMES GO ENRICHIS ####

rm(list=ls())

#### LIBRARIES ####
library(tidyverse)
library(rrvgo)
library(GOSemSim)
library(AnnotationDbi)
library(GO.db)
library(org.Mm.eg.db)

#### PATHS ####
## INPUT ## 
gsea_results <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/fgsea_obj.Rdata"
go_results <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/go_obj.Rdata"

## OUTPUT ##
# GSEA # 
reduced_gsea <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/go_gsea_analysis/gsea/reduced_terms/"
treemap_gsea <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/go_gsea_analysis/gsea/treemap"
# GO # 
reduced_go <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/go_gsea_analysis/go/reduced_terms/"
treemap_go <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/go_gsea_analysis/go/treemap"


#### SELECT CRITERIA #### 
list_ontology <- c("BP", "CC", "MF")
list_tp <- c(1,2,3)
list_reg <- c("ACC", "Hb", "Ins", "Nac")

#### LOAD DATA ####
load(gsea_results)
load(go_results)

reg <- "Hb"
tp <- "1"
ont <- "BP"

#### REDUCED TERMS FUNCTION #### 
go_rrvgo <- function(regions, timepoints, ontologies) {
  
  ## GO ## 
  for (reg in regions) { 
    message(paste0("Launch ", reg, "..."))
    for (tp in timepoints) {
      message(paste0("TP", tp, "..."))
      for (ont in ontologies) {
        message(paste0("Ontology: ", ont))
        go_name <- paste0("go_", reg, "_", tp)
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
            plot_folder <- file.path(treemap_go, reg)
            if (!dir.exists(plot_folder)) {
              dir.create(plot_folder, recursive = TRUE)
            }
            
            output_file <- file.path(plot_folder, paste0("treemap_", reg, "_tp", 
                                                         tp , "_", ont, ".png"))
            
            png(output_file, width = 8, height = 8, units = "in", res = 600)
            treemapPlot(reducedTerms)  # Utilise grid.draw pour afficher
            dev.off()
            
            # DATA #
            assign(paste("red_go", reg, tp, ont, sep = "_"), reducedTerms, envir = .GlobalEnv)
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
        
        gsea_name <- paste0("fgsea_", reg, "_", tp)
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
            plot_folder <- file.path(treemap_gsea, reg)
            if (!dir.exists(plot_folder)) {
              dir.create(plot_folder, recursive = TRUE)
            }
            
            output_file <- file.path(plot_folder, paste0("treemap_", reg, "_tp", 
                                                         tp, "_", ont, ".png"))
            
            png(output_file, width = 8, height = 8, units = "in", res = 600)
            treemapPlot(reducedTerms)  # Utilise grid.draw pour afficher
            dev.off()
            
            # DATA #
            assign(paste("red_gsea", reg, tp, ont, sep = "_"), reducedTerms, envir = .GlobalEnv)
            
            gsea_df <- gsea_res@result
            gsea_df <- gsea_df[gsea_df$ONTOLOGY == ont, ]
            
            if (nrow(gsea_df) > 1) {
              
              ## ======================
              ## Séparation UP / DOWN
              ## ======================
              
              gsea_up   <- gsea_df[gsea_df$NES > 0, ]
              gsea_down <- gsea_df[gsea_df$NES < 0, ]
              
              plot_list <- list()
              
              ## -------- UP --------
              if (nrow(gsea_up) > 1) {
                
                simMatrix_up <- calculateSimMatrix(gsea_up$ID, 
                                                   orgdb = "org.Mm.eg.db", 
                                                   ont = ont, 
                                                   method = "Rel")
                
                scores_up <- setNames(-log10(gsea_up$qvalue), gsea_up$ID)
                
                reduced_up <- reduceSimMatrix(simMatrix_up,
                                              scores_up,
                                              threshold = 0.7,
                                              orgdb = "org.Mm.eg.db")
                
                treemapPlot(reduced_up) 
                
                plot_list$up <- p_up
                
                assign(paste("red_gsea_up", reg, tp, ont, sep = "_"),
                       reduced_up, envir = .GlobalEnv)
              }
              
              ## -------- DOWN --------
              if (nrow(gsea_down) > 1) {
                
                simMatrix_down <- calculateSimMatrix(gsea_down$ID, 
                                                     orgdb = "org.Mm.eg.db", 
                                                     ont = ont, 
                                                     method = "Rel")
                
                scores_down <- setNames(-log10(gsea_down$qvalue), gsea_down$ID)
                
                reduced_down <- reduceSimMatrix(simMatrix_down,
                                                scores_down,
                                                threshold = 0.7,
                                                orgdb = "org.Mm.eg.db")
                
                p_down <- treemapPlot(reduced_down) + 
                  ggplot2::ggtitle("DOWN (NES < 0)")
                
                plot_list$down <- p_down
                
                assign(paste("red_gsea_down", reg, tp, ont, sep = "_"),
                       reduced_down, envir = .GlobalEnv)
              }
              
              ## ======================
              ## Sauvegarde PNG combiné
              ## ======================
              if (length(plot_list) > 0) {
                
                plot_folder <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/go_gsea_analysis/gsea/treemap/up_down"
                
                if (!dir.exists(plot_folder)) {
                  dir.create(plot_folder, recursive = TRUE)
                }
                
                output_file <- file.path(plot_folder,
                                         paste0("treemap_", reg, "_tp", tp, "_", ont, "_UP_DOWN.png"))
                
                png(output_file, width = 14, height = 7, units = "in", res = 600)
                
                if (length(plot_list) == 1) {
                  # Un seul plot → pas de grid.arrange
                  print(plot_list[[1]])
                } else {
                  gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
                }
                
                dev.off()
              }
                
            }
          }
        }
        
      }
    }
  }
  message(paste0("GSEA simplification done."))
}

#### RUN FUNCTION ####
go_rrvgo(list_reg, list_tp, list_ontology)


red_go_obj <- ls()[grepl("red_go_",ls())]
save(list=red_go_obj, file= paste0(reduced_go, "rrvgo_go.Rdata"))

red_gsea_obj <- ls()[grepl("red_gsea_",ls())]
save(list=red_gsea_obj, file= paste0(reduced_gsea, "rrvgo_gsea.Rdata"))