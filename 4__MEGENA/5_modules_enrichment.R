rm(list = ls()) # rm R working space

library(tidyverse)
library(GeneOverlap)
library(MEGENA)

setwd("/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/")

# Choose a region (case sensitive: ACC - Hb - Ins - Nac)
reg <- "Hb"


#### PATHS ####
alluvial_patterns.path <- paste0("data/2__differential_expression_analysis/alluvial_patterns_", reg, ".Rdata")
megena_results.path <- paste0("data/4__MEGENA/MEGENA_results/MEGENA.Results_", reg, ".Rdata")
modtable.path <- paste0("data/4__MEGENA/modtable/modtable_", reg, ".Rdata")
deglist.path <- "data/2__differential_expression_analysis/deglist.Rdata"
brain_cell_markers.path <- "data/TR_Cell_markers_for_MEGENA_annotation/Barres_lab_Cell-specific_genes_MG_21Nov2023.xlsx"
output.path <- "data/4__MEGENA/"


#### alluvial patterns enrichment in modules ####
load(alluvial_patterns.path)
load(megena_results.path)
load(modtable.path)

patterns <- df_new %>% select(MGI.symbol, diffexpressed_alltp)
genes_per_patterns <- aggregate(MGI.symbol~diffexpressed_alltp, data=patterns, FUN = function(x) paste0(x,collapse = '; '))
genes_per_patterns

data <- genes_per_patterns
colnames(genes_per_patterns)
output <- list()
for (i in 1:nrow(data)) {
  diffexpressed <- as.character(data$diffexpressed_alltp[i])
  print(diffexpressed)
  gene_list <- strsplit(data$MGI.symbol[i], "; ")[[1]]
  output[[diffexpressed]] <- gene_list
}

# Create individual df for each patterns
for (diffexpressed in names(output)) {
  assign(diffexpressed, output[[diffexpressed]])
}

pattern_list <- unique(patterns$diffexpressed_alltp)
module_list <- summary.output$modules

total.genes <- vcount(g) # Make this the total number of genes across all modules (size of the biggest module)
length(MEGENA.output$module.output)
modtable

for(pattern in pattern_list){
  overlapgenes <- get(pattern)
  overlaplist <- list()
  overlaplist[["genes"]] <- overlapgenes
  
  Object <- newGOM(overlaplist, module_list, total.genes)
  
  overlap.intersect <- getMatrix(Object, name="intersection")
  overlap.intersect.list <- getNestedList(Object, name="intersection")
  overlap.pval <- getMatrix(Object, name="pval")
  overlap.OR <- getMatrix(Object, name="odds.ratio")
  Overlap.sum <- cbind(overlap.intersect,overlap.pval,overlap.OR, overlap.intersect.list)
  Overlap.sum <- data.frame(Overlap.sum)
  Overlap.sum <- Overlap.sum %>% arrange(overlap.pval)
  
  assign(paste0("alluvial_enrichment_", pattern), Overlap.sum, envir = .GlobalEnv)
  # write.csv(Overlap.sum, file=paste0("~/Rcode/megena/enrichment_", 
  #                                    overlap, "_drn.csv"))
}
save(list=ls()[grepl("alluvial_enrichment",ls())],
     file=paste0(output.path, "enrichment_alluvial/enrichment_alluvial_", reg, ".Rdata"))




#### DEG enrichment in modules ####
load(deglist.path)
tplist <- c("tp1", "tp2", "tp3")
for(tp in tplist){
  degtp <- get(paste0("deg_", reg, "_", tp))
  degs <- degtp$label
  DEGs <- list()
  DEGs[["DEGs"]] <- degs
  
  Object <- newGOM(DEGs,module_list,total.genes)
  
  overlap.intersect <- getMatrix(Object, name="intersection")
  overlap.intersect.list <- getNestedList(Object, name="intersection")
  overlap.pval <- getMatrix(Object, name="pval")
  overlap.OR <- getMatrix(Object, name="odds.ratio")
  Overlap.sum <- cbind(overlap.intersect,overlap.pval,overlap.OR, overlap.intersect.list)
  Overlap.sum <- data.frame(Overlap.sum)
  Overlap.sum <- Overlap.sum %>% arrange(overlap.pval)
  
  assign(paste0("deg_enrichment_", reg, "_", tp), Overlap.sum, envir = .GlobalEnv)
}
save(list=ls()[grepl("deg_enrichment",ls())], 
     file=paste0(output.path, "enrichment_DEGs/enrichment_degs_", reg, ".Rdata"))




#### cellular function in modules using mithil association table ####
library(readxl)
brain_cell_markers <- read_excel(brain_cell_markers.path,
                                 sheet = "genes_>5fpkm_enrichment")

for(cell_type in colnames(brain_cell_markers)){
  genelist <- list()
  genelist[["genes"]] <- get(cell_type, brain_cell_markers)
  
  Object <- newGOM(genelist, module_list, total.genes)
  
  overlap.intersect <- getMatrix(Object, name="intersection")
  overlap.intersect.list <- getNestedList(Object, name="intersection")
  overlap.pval <- getMatrix(Object, name="pval")
  overlap.OR <- getMatrix(Object, name="odds.ratio")
  Overlap.sum <- cbind(overlap.intersect,overlap.pval,overlap.OR, overlap.intersect.list)
  Overlap.sum <- data.frame(Overlap.sum)
  Overlap.sum <- Overlap.sum %>% arrange(overlap.pval)
  
  assign(cell_type, Overlap.sum, envir = .GlobalEnv)
  # write.csv(Overlap.sum, file=paste0("~/Rcode/megena/enrichment_", 
  #                                    cell_type, "_drn.csv"))
}
save(list=c("Astrocytes_more_than_5_fpkm", "Neuron_more_than_5_fpkm", "Myelinating.Oligodendrocytes_more_than_5_fpkm", 
            "Microglia_more_than_5_fpkm", "Endothelial.Cells_more_than_5_fpkm"),
     file=paste0(output.path, "enrichment_cell_types/enrichment_cell_types_", reg, ".Rdata"))