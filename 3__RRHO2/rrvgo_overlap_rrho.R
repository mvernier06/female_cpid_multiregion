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

base_dir <- "data/3__RRHO2/overlap_GO_results"
out_dir  <- "data/3__RRHO2/rrvgo_overlaps_rrho"
out_dir_graph <-  "graphs_results/3__RRHO2/gene_overlap_rrho/"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_graph)

# récupérer tous les fichiers GO
files <- list.files(base_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

all_results <- list()

for (f in files) {
  
  ego <- readRDS(f)
  
  if (is.null(ego) || nrow(ego@result) == 0) next
  
  df <- ego@result
  
  # garder termes significatifs (tu peux ajuster)
  df <- df %>% filter(p.adjust < 0.05)
  
  if (nrow(df) < 2) next
  
  # metadata depuis nom fichier
  name <- basename(f)
  name <- gsub("_GO.rds", "", name)
  
  df$comparison <- name
  
  all_results[[name]] <- df
}

# combiner
go_all <- bind_rows(all_results)

# sauver tableau global

write.xlsx(go_all,
           file = file.path(out_dir, "GO_summary.xlsx"))
save(all_results, file=file.path(out_dir,"all_results.Rdata"))

# boucle par comparaison

for (comp in unique(go_all$comparison)) {
  
  message("Processing ", comp)
  
  df <- go_all %>% filter(comparison == comp)
  
  if (grepl("^RRHO_", comp)) {
    
    # multitp
    
    parts <- str_split(comp, "_")[[1]]
    
    type      <- "multitp"
    region    <- parts[2]
    comp_tp   <- parts[3]
    quadrant  <- paste(parts[4], parts[5], sep = "_")
    
    subdir <- file.path(type, region, comp_tp, quadrant)
    
  } else if (grepl("^rrho", comp)) {
    
    # multireg
    
    tmp <- str_remove(comp, "^rrho")
    
    # ACCvsHb1
    main <- str_split(tmp, "_")[[1]][1]
    
    quadrant <- str_split(tmp, "_")[[1]][2:3] %>% paste(collapse = "_")
    
    # séparer régions et TP
    reg_part <- str_extract(main, "^[A-Za-z]+vs[A-Za-z]+")
    tp       <- str_extract(main, "[0-9]+$")
    
    subdir <- file.path("multireg",
                        paste0("TP", tp),
                        reg_part,
                        quadrant)
  }
  
  ## créer dossier
  full_outdir <- file.path(out_dir_graph, subdir)
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
    png(file.path(full_outdir, paste0("treemap_", ont, ".png")),
        width = 2000, height = 2000, res = 300)
    treemapPlot(reducedTerms)
    dev.off()
  }
}
