#### VISUALIZATION GSEA #### 
#### ATTENTION, beaucoup de modifications par raport au script d'origine ####

rm(list=ls())

dir.create("~/R/library", showWarnings = FALSE, recursive = TRUE)
.libPaths(c("~/R/library", .libPaths()))

Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always") 

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
library(GOSemSim)

#### PATHS ####
## INPUT ##
go <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/go_obj.Rdata"
# go_1 <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/go_1_obj.Rdata"

## OUTPUT ##
data <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/"
plot.file <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/go_gsea_analysis/go/"

#### LOAD DATA ####
load(go)
# load(go_1)

#### CRITERIA ####
regionList <- c("ACC", "Ins", "Hb", "Nac")
tpList <- c(1, 2, 3)

################################################################################
######                               DOTPLOTS                             ######
################################################################################

get_dotplot <- function(reg, tp, x_axis, nb_terms, plot_path) { 
  
  obj_name <- paste0("go_", reg, "_", tp)
  
  if (!exists(obj_name)) {
    message("Objet ", obj_name, " inexistant")
    return(NULL)
  }
  
  res_go <- get(obj_name)
  
  if (is.null(res_go) || nrow(res_go@result) == 0) {
    message("Aucun résultat GO pour ", reg, " TP", tp)
    return(NULL)
  }
  
  # Ajout FoldEnrichment 
  if (x_axis == "FoldEnrichment") {
    
    df <- res_go@result
    
    geneRatio <- sapply(strsplit(df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
    bgRatio   <- sapply(strsplit(df$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
    
    res_go@result$FoldEnrichment <- geneRatio / bgRatio
  }
  
  nb_terms <- min(nb_terms, nrow(res_go@result))
  
  dotplot_go <- dotplot(
    res_go,
    showCategory = nb_terms,
    x = x_axis,
    title = paste0("Dotplot of GO results - ", reg, " TP", tp)
  )
  
  output_dir <- file.path(plot_path, "dotplots", reg, x_axis)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(
    filename = paste0(output_dir, "/dotplot_", x_axis, "_", reg, "_tp", tp, ".png"),
    plot = dotplot_go,
    bg = "white",
    width = 1900,
    height = 1200,
    units = "px",
    scale = 2
  )
}


#### PARAMETERS ####
params <- expand.grid(reg = regionList, tp = tpList)

#### APPLY FUNCTION #### 
apply(params, 1, function(row) {
  get_dotplot(
    reg = row["reg"], 
    tp = row["tp"], 
    x_axis = "geneRatio", 
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
    
    tp_tables <- list()
    
    for (tp in tp_list) {
      
      obj_name <- paste0("go_", reg, "_", tp)
      
      if (!exists(obj_name, envir = .GlobalEnv)) {
        message("L'objet ", obj_name, " n'existe pas. Création d'un TP vide.")
        df <- data.frame(Description = character(0),
                         RichFactor = numeric(0),
                         FoldEnrichment = numeric(0),
                         pvalue = numeric(0),
                         p.adjust = numeric(0),
                         qvalue = numeric(0))
      } else {
        res_go <- get(obj_name, envir = .GlobalEnv)
        if (is.null(res_go) || nrow(res_go@result) == 0) {
          message("Aucun résultat pour ", obj_name, ". Création d'un TP vide.")
          df <- data.frame(Description = character(0),
                           RichFactor = numeric(0),
                           FoldEnrichment = numeric(0),
                           pvalue = numeric(0),
                           p.adjust = numeric(0),
                           qvalue = numeric(0))
        } else {
          df <- as.data.frame(res_go)
          
          # Calcul FoldEnrichment
          geneRatio <- sapply(strsplit(df$GeneRatio, "/"),
                              function(x) as.numeric(x[1])/as.numeric(x[2]))
          bgRatio <- sapply(strsplit(df$BgRatio, "/"),
                            function(x) as.numeric(x[1])/as.numeric(x[2]))
          df$FoldEnrichment <- geneRatio / bgRatio
          
          # Calcul RichFactor
          bg_size <- sapply(strsplit(df$BgRatio, "/"),
                            function(x) as.numeric(x[1]))
          df$RichFactor <- df$Count / bg_size
          
          # Sélection
          df <- df %>%
            dplyr::select(Description, RichFactor, FoldEnrichment, pvalue, p.adjust, qvalue)
        }
      }
      
      # Rename colonnes avec TP
      colnames(df) <- c("Description",
                        paste0("RichFactor_tp", tp),
                        paste0("FoldEnrichment_tp", tp),
                        paste0("pvalue_tp", tp),
                        paste0("p.adjust_tp", tp),
                        paste0("qvalue_tp", tp))
      
      tp_tables[[as.character(tp)]] <- df
    }
    
    # Fusionner tous les TP
    go_alltp <- Reduce(function(x, y)
      dplyr::full_join(x, y, by = "Description"),
      tp_tables)
    
    assign(paste0("go_", reg, "_alltp"),
           go_alltp,
           envir = .GlobalEnv)
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
      
      obj_name <- paste0("go_", reg, "_alltp")
      
      if (!exists(obj_name, envir = .GlobalEnv)) {
        message(paste("The object", obj_name, "doesn't exist in this environment."))
        next  
      }
      
      gsea_alltp <- get(obj_name, envir = .GlobalEnv)
      
      if (!"Description" %in% colnames(gsea_alltp)) {
        message(paste("The column 'Description' is missing in the object", obj_name))
        next
      }
      
      # Filtrer les termes contenant le mot-clé
      test <- gsea_alltp %>% 
        filter(grepl(term, Description, ignore.case = TRUE)) %>%
        dplyr::select(Description, contains("RichFactor"))   # Modifier si une autre variable est souhaitée
      
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
      
      test[,-1] <- lapply(test[,-1], function(x) round(x, digits = 2))
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
      region_folder <- file.path(plot.path, "heatmaps", "RichFactor", reg)
      if (!dir.exists(region_folder)) {
        dir.create(region_folder, recursive = TRUE)
      }
      
      output_file <- file.path(region_folder, paste0("heatmap_", reg, "_", term, ".png"))
      
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
               main = paste0("Enrichment of GO terms in ", reg, " (using the term ", term, ", -- RichFactor)"),
               legend_breaks = seq(-maxval, maxval, length.out = 3),
               angle_col = 45
      )
      
      dev.off()
      message(paste("Heatmap saved:", output_file))
    }
  }
}


go_heatmap(regionList, terms_list, plot.file)



################################################################################
######                              UPSETR                                ######
################################################################################

reg <- "Hb"
plot.path <- plot.file
upset_plot <- function(reg, plot.path) { 
  # Chargement des packages nécessaires
  library(UpSetR)
  library(ggplot2)
  library(dplyr)
  
  # Fonction pour calculer RichFactor à partir d'un objet enrichGO
  calcRichFactor <- function(go_obj) {
    if (is.null(go_obj) || nrow(go_obj@result) == 0) return(data.frame())
    df <- as.data.frame(go_obj)
    
    # Calcul RichFactor
    bg_size <- sapply(strsplit(df$BgRatio, "/"), function(x) as.numeric(x[1]))
    df$RichFactor <- df$Count / bg_size
    df
  }
  
  # Récupération et préparation des TP
  tp_list <- 1:3
  df_list <- list()
  
  for (tp in tp_list) {
    obj_name <- paste0("go_", reg, "_", tp)
    
    if (!exists(obj_name, envir = .GlobalEnv)) {
      message(paste("L'objet", obj_name, "n'existe pas. TP", tp, "sera ignoré."))
      df_list[[tp]] <- data.frame(Description = character(0), RichFactor = numeric(0))
    } else {
      go_obj <- get(obj_name, envir = .GlobalEnv)
      df_list[[tp]] <- calcRichFactor(go_obj)
    }
  }
  
  # Préparer la liste pour UpSetR
  listInput <- list(
    tp1 = if(nrow(df_list[[1]])>0) df_list[[1]]$Description else character(0),
    tp2 = if(nrow(df_list[[2]])>0) df_list[[2]]$Description else character(0),
    tp3 = if(nrow(df_list[[3]])>0) df_list[[3]]$Description else character(0)
  )
  
  # Vérifier que la liste contient au moins un terme
  if (all(sapply(listInput, length) == 0)) {
    message(paste("Aucun terme enrichi pour la région", reg, ". Plot non généré."))
    return(NULL)
  }
  
  # Créer le dossier de sortie
  region_folder <- file.path(plot.path, "upset_plot")
  if (!dir.exists(region_folder)) dir.create(region_folder, recursive = TRUE)
  output_file <- file.path(region_folder, paste0("upsetr_go_", reg, ".png"))
  
  # Générer le plot
  png(filename = output_file, height = 1250, width = 1250, units = "px", res = 225)
  
  print(
    upset(
      fromList(listInput),
      mainbar.y.label = paste0("GO terms enriched - Intersections in ", reg),
      sets.x.label = "Enriched terms per TP",
      order.by = "freq"
    )
  )
  
  dev.off()
  
  
  # Vérifier le fichier généré
  if (file.exists(output_file) && file.info(output_file)$size > 0) {
    message(paste("UpSet plot sauvegardé pour", reg, ":", output_file))
  } else {
    message(paste("Erreur : le fichier est vide pour", reg))
  }
  
  return(output_file)
}


upset_plot(reg = "ACC", plot.path = plot.file) 
upset_plot(reg = "Ins", plot.path = plot.file)
upset_plot(reg = "Nac", plot.path = plot.file)
upset_plot(reg = "Hb", plot.path = plot.file)

# # Exécution de la fonction pour les six régions 
# results <- lapply(regionList, function(reg) {
#   tryCatch({
#     upset_plot(reg = reg, plot.path = plot.file)
#   }, error = function(e) {
#     cat("Erreur pour la région", reg, ":", conditionMessage(e), "\n")
#     NULL
#   })
# })

################################################################################
######                              GO Network                                ######
################################################################################

## test 
ego <- go_Hb_1

# Calcul similarité sémantique
ego <- pairwise_termsim(ego)

emapplot(
  ego,
  showCategory = 15,   # nombre de termes affichés
  layout = "kk"        # kamada-kawai (plus joli)
)

cnetplot(
  ego,
  showCategory = 10,
  foldChange = NULL,  # ou un vecteur log2FC nommé
  circular = FALSE,
  colorEdge = TRUE
)



