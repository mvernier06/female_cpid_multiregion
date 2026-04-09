#### VISUALIZATION GSEA #### 

rm(list=ls())

#### LIBRARIES ####
library(tidyverse)
library(dplyr)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(pheatmap)
library(grid)  
library(gtable) 
library(UpSetR)

#### PATHS ####
project.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project.path)

## INPUT ##
go_05 <- "data/2__differential_expression_analysis/panregion_common_regions/go_gsea_analysis/go_panregion_05_obj.Rdata"
# go_1 <- "data/2__differential_expression_analysis/panregion_common_regions/go_gsea_analysis/go_panregion_1_obj.Rdata"

## OUTPUT ##
data <- "data/2__differential_expression_analysis/panregion_common_regions/go_gsea_analysis/"
plot.file <- "graphs_results/panregion_common_regions/go_gsea_analysis/go/"

#### LOAD DATA ####
load(go_05)
# load(go_1)

#### CRITERIA ####
regionList <- "panregion"
tpList <- c(1, 2, 3)

################################################################################
######                               DOTPLOTS                             ######
################################################################################

get_dotplot <- function(reg, tp, x_axis, nb_terms, plot_path) { 
  
  # Récupère les résultats GO pour la région et le timepoint spécifiés
  res_go <- get(paste0("go_", reg, "_05_", tp))
  
  # Vérifie s'il y a des résultats GO
  if (nrow(res_go) == 0) {
    message(paste("Aucun résultat GO pour", reg, "TP", tp))
    return(NULL)  # Si aucun résultat, on quitte la fonction sans créer le graphique
  }
  
  # Limite le nombre de termes si nécessaire
  if (nrow(res_go) < nb_terms) {
    nb_terms <- nrow(res_go)  # Adapte nb_terms à la taille du jeu de données disponible
  }
  
  # Crée le dotplot avec les résultats GO
  dotplot_go <- dotplot(res_go, showCategory = nb_terms, x = x_axis, 
                        title = paste0("Dotplot of GO results 0.5 - ", reg, " TP", tp))
  
  # Crée le répertoire de sortie s'il n'existe pas
  output_dir <- file.path(plot_path, "dotplots",  x_axis)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Sauvegarde le graphique en PNG
  ggsave(filename = paste0(output_dir, "/dotplot_", x_axis, "_", reg, "_05_tp", tp, 
                           ".png"), 
         plot = dotplot_go, 
         bg = "white", 
         width = 1900, height = 1200, units="px", scale=2)
}

#### PARAMETERS ####
params <- expand.grid(reg = regionList, tp = tpList)

#### APPLY FUNCTION #### 
apply(params, 1, function(row) {
  get_dotplot(
    reg = row["reg"], 
    tp = row["tp"], 
    x_axis = "GeneRatio", 
    nb_terms = 15,
    plot_path = plot.file
  )
})


################################################################################
######                              HEATMAPS                              ######
################################################################################

terms_list <- c("myelin", "calcium", "oligo", "neuron", "axon", "micro", 
                "neurofil", "sero", "opio", "dendri", "astrocyte", "GABA", "sheath")

#### LOOP TO HAVE HEATMAP FILES #### 
heatmap_files <- function(region_list, tp_list) {
  
  for (reg in region_list) {
    for (tp in tp_list) {
      
      # Créer dynamiquement le nom de l'objet à partir de la région et du time point
      obj_name <- paste0("go_", reg, "_05_", tp)
      
      # Vérifier si l'objet existe dans l'environnement
      if (exists(obj_name)) {
        
        # Récupérer l'objet correspondant à partir de l'environnement
        file_tp <- get(obj_name, envir = .GlobalEnv)
        file_tp <- as.data.frame(file_tp)
        # Sélectionner les colonnes d'intérêt
        file_tp <- file_tp %>% 
          dplyr::select(Description, GeneRatio, BgRatio,pvalue, p.adjust, qvalue)
        
        # Renommer les colonnes
        colnames(file_tp) <- c("Description", 
                               paste0("GeneRatio_tp", tp),
                               paste0("BgRatio_tp", tp),
                               paste0("pvalue_tp", tp),
                               paste0("p.adjust_tp", tp),
                               paste0("qvalue_tp", tp))
        
        # Sauvegarder le tableau dans l'environnement global
        assign(paste0(reg, "_05_tp", tp), file_tp, envir = .GlobalEnv)
      } else {
        message(paste("L'objet", obj_name, "n'existe pas dans l'environnement."))
      }
    }
    
    # Fusionner les différents time points (tp1, tp2, tp3) pour la région donnée
    go_alltp <- Reduce(function(x, y) full_join(x, y, by = "Description"),
                       lapply(tp_list, function(tp) get(paste0(reg, "_05_tp", tp), envir = .GlobalEnv)))
    
    # Sauvegarder le résultat final pour la région
    assign(paste0("go_", reg, "_05_alltp"), go_alltp, envir = .GlobalEnv)
  }
}


#### LUNCH LOOP ####
heatmap_files(regionList, tpList)



#### FUNCTION TO GENERATE HEATMAPS ####

go_heatmap <- function(region_list, terms_list, plot.path, top10 = FALSE, max_line_length = 30) {
  
  for (reg in region_list) {  
    message(paste(reg, "..."))
    
    for (term in terms_list) {
      message(paste("Searching for", term, "in enriched terms..."))
      
      obj_name <- paste0("go_", reg, "_05_alltp")
      
      if (!exists(obj_name, envir = .GlobalEnv)) {
        message(paste("The object", obj_name, "doesn't exist in this environment."))
        next  
      }
      
      go_alltp <- get(obj_name, envir = .GlobalEnv)
      
      if (!"Description" %in% colnames(go_alltp)) {
        message(paste("The column 'Description' is missing in the object", obj_name))
        next
      }
      
      # Filtrer les termes contenant le mot-clé
      test <- go_alltp %>% 
        filter(grepl(term, Description, ignore.case = TRUE)) %>%
        dplyr::select(Description, contains("GeneRatio"))   # Modifier si une autre variable est souhaitée
      
      # Log transformation
      #test[,-1] <- lapply(test[,-1], function(x) -log10(x))
      
      if (nrow(test) == 0) {
        message(paste("No terms corresponding to", term, "in", obj_name))
        next
      }
      
      if (top10) test <- head(test, 10)
      
      # Enjoliver les noms longs
      test$Description <- sapply(test$Description, 
                                 function(x) paste(strwrap(x, width = max_line_length), collapse = "\n"))
      
      test <- test %>%
        mutate(across(-Description, ~ {
          if (is.character(.)) {
            sapply(strsplit(., "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
          } else {
            as.numeric(.)
          }
        })) %>%
        mutate(across(-Description, ~ round(., 2)))
      test2 <- test %>% column_to_rownames("Description")
      
      # Déterminer la valeur maximale pour la mise à l’échelle
      maxval <- test2 %>% 
        apply(2, function(x) replace_na(x, 0)) %>%
        as.matrix() %>% abs() %>% max()
      
      # Protéger contre les erreurs dues à maxval = 0
      if (maxval == 0 || is.na(maxval)) {
        message(paste("Skipping heatmap for", reg, "and term", term, "because maxval == 0 or NA."))
        next
      }
      
      # Créer les dossiers si besoin
      region_folder <- file.path(plot.path, "heatmaps", "GeneRatio")
      if (!dir.exists(region_folder)) {
        dir.create(region_folder, recursive = TRUE)
      }
      
      output_file <- file.path(region_folder, paste0("heatmap_", reg, "_05_", term, ".png"))
      
      png(output_file, width = 1400, height = 1200, units = "px", res = 300)
      
      breaks <- seq(-maxval, maxval, length.out = 100)
      
      pheatmap(test2,
               display_numbers = TRUE,
               fontsize = 5,
               fontsize_number = 6,
               breaks = breaks,
               cluster_cols = FALSE,
               cluster_rows = FALSE,
               labels_col = c("TP1", "TP2", "TP3"),
               main = paste0("Enrichment of GO terms in ", reg, " 05 (using the term ", term, ", -- GeneRatio)"),
               legend_breaks = seq(-maxval, maxval, length.out = 3),
               angle_col = 45
      )
      
      dev.off()
      message(paste("Heatmap saved:", output_file))
    }
  }
}


go_heatmap(regionList, terms_list, plot.file)
