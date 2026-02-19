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
## INPUT ##
gsea <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/fgsea_obj.Rdata"
# gsea_1 <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/fgsea_1_obj.Rdata"
go <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/go_obj.Rdata"

## OUTPUT ##
data <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/"
plot.file <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/go_gsea_analysis/gsea/"

#### LOAD DATA ####
load(gsea)
# load(gsea_1)
load(go)

#### CRITERIA ####
regionList <- c("ACC",  "Ins", "Hb", "Nac")
tpList <- c(1, 2, 3)

################################################################################
######                               DOTPLOTS                             ######
################################################################################

get_dotplot <- function(reg, tp, x_axis, nb_terms, plot_path) { 
  
  obj_name <- paste0("fgsea_", reg, "_", tp)
  
  if (!exists(obj_name)) {
    message("Object ", obj_name, " does not exist")
    return(NULL)
  }
  
  res_fgsea <- get(obj_name)
  
  # Vérifie si résultats non vides
  if (is.null(res_fgsea) || nrow(res_fgsea@result) == 0) {
    message("No enriched terms for ", reg, " TP", tp)
    return(NULL)
  }
  
  dotplot_gsea <- dotplot(res_fgsea, 
                          showCategory = nb_terms, 
                          x = x_axis,
                          title = paste0("Dotplot of FGSEA results - ", reg, " TP", tp))
  
  output_dir <- file.path(plot_path, "dotplots", reg, x_axis)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  ggsave(
    filename = paste0(output_dir, "/dotplot_", x_axis, "_", reg, "_tp", tp, ".png"),
    plot = dotplot_gsea,
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
    x_axis = "GeneRatio", # ca fonctionne car dotplot calcul GeneRatio = Count / setSize en interne 
    nb_terms = 15,
    plot_path = plot.file
  )
})

#### APPLY FUNCTION #### 
apply(params, 1, function(row) {
  get_dotplot(
    reg = row["reg"], 
    tp = row["tp"], 
    x_axis = "NES", 
    nb_terms = 15,
    plot_path = plot.file
  )
})


################################################################################
######                          DOTPLOT ODD RATIO                         ######
################################################################################

### Pourquoi faire du enrichGO dans la partie gsea ??

# odds_ratio <- function(reg, tp, x_axis, plot_path) {
#   # S'assurer que le dossier existe
#   output_dir <- file.path(plot_path, "dotplots", reg, x_axis)
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir, recursive = TRUE)
#   }
#   
#   # Récupérer l'objet FGSEA 
#   enr_name <- paste0("go_", reg, "_", tp)
#   res_enrich <- get(enr_name)
#   
#   # Vérifier si l'objet existe
#   if (!exists(enr_name)) {
#     message(paste("Objet", enr_name, "n'existe pas."))
#     return(NULL)
#   }
#   
#   res_enrich <- get(enr_name)
#   
#   # Vérifier s'il y a des résultats significatifs
#   df <- as.data.frame(res_enrich)
#   if (nrow(df) == 0) {
#     message(paste("Aucun résultat dans", enr_name, "- on skip."))
#     return(NULL)
#   }
#   
#   # Transformer en data frame pour calculer les odds ratio
#   df <- as.data.frame(res_enrich)
#   
#   # Évaluer les ratios en nombre
#   eval_ratio <- function(ratio_str) {
#     parts <- strsplit(ratio_str, "/")[[1]]
#     as.numeric(parts[1]) / as.numeric(parts[2])
#   }
#   
#   df$GeneRatio_num <- unlist(sapply(strsplit(df$GeneRatio, "/"), 
#                                     function(x) as.numeric(x[1]) / as.numeric(x[2])))
#   df$BgRatio_num   <- unlist(sapply(strsplit(df$BgRatio, "/"),   
#                                     function(x) as.numeric(x[1]) / as.numeric(x[2])))
#   
#   # Calcul de l’odds ratio
#   df$odds_ratio <- (df$GeneRatio_num / (1 - df$GeneRatio_num)) /
#     (df$BgRatio_num / (1 - df$BgRatio_num))
#   
#   # Remettre les données dans l'objet enrichResult avec les nouvelles colonnes
#   res_enrich@result <- df
#   
#   # Faire le dotplot avec l'odds_ratio sur l’axe X
#   dotplot_gsea <- dotplot(
#     res_enrich,
#     showCategory = 15,
#     x = "odds_ratio",
#     title = paste0("Dotplot (Odds Ratio) - ", reg, " TP", tp)
#   )
#   
#   # Sauvegarde
#   ggsave(
#     filename = file.path(output_dir, paste0("dotplot_odds_ratio_", reg, "_tp", 
#                                             tp, ".png")),
#     plot = dotplot_gsea,
#     bg = "white",
#     width = 1900,
#     height = 1200,
#     units = "px",
#     scale = 2
#   )
# }
# 
# 
# #### PARAMETERS ####
# regionList <- c("ACC", "Ins", "Hb", "Nac")
# tpList <- c(1, 2, 3)
# params <- expand.grid(reg = regionList, tp = tpList)
# 
# #### APPLY FUNCTION #### 
# apply(params, 1, function(row) {
#   odds_ratio(
#     reg = row[["reg"]], 
#     tp = as.numeric(row[["tp"]]),
#     x_axis = "odds ratio", 
#     plot_path = plot.file
#   )
# })



################################################################################
######                            ENRICHPLOT                              ######
################################################################################

# get_enrichplot <- function(region, tp, term, plot.path) {
# 
#   fgsea_res <- get(paste0("fgsea_", region, "_", tp), envir = .GlobalEnv)
#   fgsea_df <- as.data.frame(fgsea_res)
#   # get term index
#   term_index <- which(fgsea_df$Description == term)
# 
#   # plot
#   plot <- gseaplot2(fgsea_res, geneSetID = term_index,
#                     title = paste0(term, " - ", region, " TP", tp ))
#   # save
#   region_dir <- file.path(plot.path, region, "enrichplot")
#   if (!dir.exists(region_dir)) dir.create(region_dir, recursive = TRUE)
# 
#   ggsave(file= paste0(region_dir, term, "_", region, "_tp", tp, ".png"), plot=plot,
#          bg="white", width=1900, height=1200, units="px", scale=2)
# 
# }
# 
# get_enrichplot(region = "Hb", tp = 1, term = "compact myelin", plot.path = plot.file)
# 

# #### GET GENES (EXCEL FILE) ####


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
      obj_name <- paste0("fgsea_", reg, "_", tp)
      
      # Vérifier si l'objet existe dans l'environnement
      if (exists(obj_name)) {
        
        # Récupérer l'objet correspondant à partir de l'environnement
        file_tp <- get(obj_name, envir = .GlobalEnv)
        file_tp <- as.data.frame(file_tp)
        # Sélectionner les colonnes d'intérêt
        file_tp <- file_tp %>% 
          dplyr::select(Description, NES, pvalue, p.adjust, qvalue)
        
        # Renommer les colonnes
        colnames(file_tp) <- c("Description", 
                               paste0("NES_tp", tp),
                               paste0("pvalue_tp", tp),
                               paste0("p.adjust_tp", tp),
                               paste0("qvalue_tp", tp))
        
        # Sauvegarder le tableau dans l'environnement global
        assign(paste0(reg, "_tp", tp), file_tp, envir = .GlobalEnv)
      } else {
        message(paste("L'objet", obj_name, "n'existe pas dans l'environnement."))
      }
    }
    
    # Fusionner les différents time points (tp1, tp2, tp3) pour la région donnée
    gsea_alltp <- Reduce(function(x, y) full_join(x, y, by = "Description"),
                         lapply(tp_list, function(tp) get(paste0(reg, "_tp", tp), envir = .GlobalEnv)))
    
    # Sauvegarder le résultat final pour la région
    assign(paste0("gsea_", reg, "_alltp"), gsea_alltp, envir = .GlobalEnv)
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
      
      obj_name <- paste0("gsea_", reg, "_alltp")
      
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
        dplyr::select(Description, contains("pval"))   # Modifier si une autre variable est souhaitée
      
      # Log transformation
      test[,-1] <- lapply(test[,-1], function(x) -log10(x))  # to comment if not pval or padj
      
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
      region_folder <- file.path(plot.path, "heatmaps", "pval", reg)
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
               main = paste0("Enrichment of GO terms in ", reg, " (using the term ", term, ", -- -log10(pval))"),
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


upset_plot <- function(reg, plot.path) { 
  library(UpSetR)
  library(ggplot2)
  
  fgsea_reg_1 <- get(paste0("fgsea_", reg, "_1"), envir = .GlobalEnv)
  fgsea_reg_2 <- get(paste0("fgsea_", reg, "_2"), envir = .GlobalEnv)
  fgsea_reg_3 <- get(paste0("fgsea_", reg, "_3"), envir = .GlobalEnv)
  
  listInput <- list(
    tp1_up = fgsea_reg_1@result$Description[which(sign(fgsea_reg_1@result$NES)>0)],
    tp1_down = fgsea_reg_1@result$Description[which(sign(fgsea_reg_1@result$NES)<0)],
    tp2_up = fgsea_reg_2@result$Description[which(sign(fgsea_reg_2@result$NES)>0)],
    tp2_down = fgsea_reg_2@result$Description[which(sign(fgsea_reg_2@result$NES)<0)],
    tp3_up = fgsea_reg_3@result$Description[which(sign(fgsea_reg_3@result$NES)>0)],
    tp3_down = fgsea_reg_3@result$Description[which(sign(fgsea_reg_3@result$NES)<0)]
  )
  
  region_folder <- file.path(plot.path, "upster_plot")
  if (!dir.exists(region_folder)) {
    dir.create(region_folder, recursive = TRUE)
  }
  
  output_file <- file.path(region_folder, paste0("upsetr_fgsea_", reg, ".png"))
  
  # Vérifier le contenu de listInput
  print("Contenu de listInput :")
  print(lapply(listInput, length))
  
  # Utiliser quartz() pour MacOS ou windows() pour Windows si nécessaire
  png(filename = output_file, height = 1250, width = 1250, units = "px", res = 225)
  
  # Ajouter des informations de débogage
  tryCatch({
    upset_plot <- upset(fromList(listInput), 
                        mainbar.y.label = paste0("GSEA : GO termes enriched - Intersections in the ", reg), 
                        sets.x.label = "Number of enriched terms per timepoint", 
                        order.by="freq")
    
    # Imprimer le plot
    print(upset_plot)
  }, error = function(e) {
    cat("Erreur lors de la création du plot :", conditionMessage(e), "\n")
  })
  
  dev.off()
  
  # Vérifier si le fichier existe et n'est pas vide
  file_info <- file.info(output_file)
  if (file_info$size > 0) {
    cat("Le plot a été sauvegardé avec succès à l'emplacement :", output_file, "\n")
  } else {
    cat("Erreur : Le fichier est vide\n")
  }
  
  return(output_file)
}

# Exécution de la fonction
upset_plot(reg = "ACC", plot.path = plot.file)
upset_plot(reg = "Ins", plot.path = plot.file)
upset_plot(reg = "Hb", plot.path = plot.file)
upset_plot(reg = "Nac", plot.path = plot.file)

# # Exécution de la fonction pour les six régions 
# results <- lapply(regionList, function(reg) {
#   tryCatch({
#     upset_plot(reg = reg, plot.path = plot.file)
#   }, error = function(e) {
#     cat("Erreur pour la région", reg, ":", conditionMessage(e), "\n")
#     NULL
#   })
# })

