.libPaths(c("~/R/library", .libPaths()))
Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always") 
library(tidyverse)
library(clusterProfiler)
library(dplyr)

rm(list=ls())

project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
counts_path <- "data/2__differential_expression_analysis/annotated_counts_filtered.rds"
setwd(project_path)

load("data/3__RRHO2/rrho_obj_multireg.Rdata")
load("data/3__RRHO2/rrho_obj_multitp.Rdata")
counts <- readRDS(counts_path)
organism <- "org.Mm.eg.db"
library(organism, character.only = TRUE)

### Extraire les meilleurs overlaps par quadrant
regionList <- c("ACC", "Hb", "Ins", "Nac")
tp_comparisons <- c("1vs2", "1vs3", "2vs3")
tpList <- c(1, 2, 3)

write_rrho_quadrant_overlaps <- function(rrho, outdir, prefix) {
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  overlaps <- list(
    up_up     = rrho$genelist_uu$gene_list_overlap_uu,
    down_down = rrho$genelist_dd$gene_list_overlap_dd,
    down_up   = rrho$genelist_du$gene_list_overlap_du,
    up_down   = rrho$genelist_ud$gene_list_overlap_ud
  )
  
  for (q in names(overlaps)) {
    outfile <- file.path(
      outdir,
      paste0(prefix, "_", q, ".txt")
    )
    
    writeLines(overlaps[[q]], outfile)
  }
}

base_dir <- "data/3__RRHO2/gene_overlap/multitp"

for (reg in regionList) {
  for (tp in tp_comparisons) {
    
    obj_name <- paste0("RRHO_", reg, "_", tp)
    rrho <- get(obj_name)
    
    outdir <- file.path(base_dir, reg, tp)
    
    write_rrho_quadrant_overlaps(
      rrho   = rrho,
      outdir = outdir,
      prefix = obj_name
    )
  }
}

base_dir <- "data/3__RRHO2/gene_overlap/multireg"


for (tp in tpList) {
  for (i in seq_len(length(regionList) - 1)) {
    for (j in (i + 1):length(regionList)) {
      
      reg1 <- regionList[i]
      reg2 <- regionList[j]
      
      obj_name <- paste0("rrho", reg1, "vs", reg2, tp)
      rrho <- get(obj_name)
      
      outdir <- file.path(
        base_dir,
        paste0("TP", tp),
        paste0(reg1, "_vs_", reg2)
      )
      
      write_rrho_quadrant_overlaps(
        rrho   = rrho,
        outdir = outdir,
        prefix = obj_name
      )
    }
  }
}

### Tester l'enrichissement de ces overalps 
#########################################
### Attention ce script est très long ###
#########################################
# Dossier racine
base_dirs <- c(
  "data/3__RRHO2/gene_overlap/multitp",
  "data/3__RRHO2/gene_overlap/multireg"
)

# Fonction GO
run_go <- function(genes, rrho_universe ) {
  
  ego <- enrichGO(gene = genes, 
                  OrgDb = organism, 
                  keyType = "SYMBOL", 
                  ont ="ALL",
                  qvalueCutoff = 1,
                  pAdjustMethod = "BH",
                  universe = rrho_universe)
  
  return(ego)
}

# Parcourir tous les fichiers
for (base_dir in base_dirs) {
  
  files <- list.files(
    base_dir,
    pattern = "\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  for (f in files) {
    
    print(f)
    # Lire gènes
    genes <- readLines(f)
    print(length(genes))
    
    # skip si trop peu de gènes
    if (length(genes) < 2) next
    
    # GO
    ego <- run_go(genes, counts$MGI.symbol )
    print(length(ego@result$Description))
    
    # skip si pas de résultat
    if (is.null(ego) || nrow(ego) == 0) next
    
    # output file
    outfile <- sub(
      "gene_overlap",
      "overlap_GO_results",
      f
    )
    
    outfile <- sub("\\.txt$", "_GO.rds", outfile)
    
    print(outfile)
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
     
    saveRDS(ego, outfile)
  }
}
