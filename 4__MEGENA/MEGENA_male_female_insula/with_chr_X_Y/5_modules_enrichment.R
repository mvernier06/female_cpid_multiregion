### Le script a été éxécuté dans un environnement conda R4.4 pour faire fonctionner GeneOverlap ###

rm(list = ls()) # rm R working space
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", version = '3.17', repos = "https://cloud.r-project.org")

# GeneOverlap 1.40
BiocManager::install("GeneOverlap") 

# MEGENA
install.packages("MEGENA", repos = "https://cloud.r-project.org")  # ou devtools::install_github("jasonbian97/MEGENA")
library(tidyverse)
library(GeneOverlap)
library(MEGENA)

setwd("/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/")

# Choose a region : Ins
reg <- "Ins"
print(paste0("Running module enrichment for region: ", reg))


#### PATHS ####
# alluvial_patterns.path <- paste0("data/2__differential_expression_analysis/alluvial_patterns_", reg, ".Rdata")
megena_results.path <- paste0("data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/MEGENA.Results_", reg, ".Rdata")
modtable.path <- paste0("data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/modtable_", reg, ".Rdata")
# deglist.path <- "data/2__differential_expression_analysis/deglist.Rdata"
brain_cell_markers.path <- "data/4__MEGENA/TR_Cell_markers_for_MEGENA_annotation/Barres_lab_Cell-specific_genes_MG_21Nov2023.xlsx"
genes_chrXY.path <- "data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/genes_chrXY_MGI_symbols.txt"
output.path <- "data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/"

print(paste0("Running module enrichment for region: ", reg))
#### alluvial patterns enrichment in modules ####
# load(alluvial_patterns.path) # df_new
load(megena_results.path)
load(modtable.path)

genes_chrXY <- read.table(
  genes_chrXY.path,
  header = TRUE,
  stringsAsFactors = FALSE
)
genelist_chrXY <- list()
genelist_chrXY[["genes"]] <- genes_chrXY$MGI.symbol

# patterns <- df_new %>% select(MGI.symbol, diffexpressed_alltp)
# genes_per_patterns <- aggregate(MGI.symbol~diffexpressed_alltp, data=patterns, FUN = function(x) paste0(x,collapse = '; '))
# genes_per_patterns

# data <- genes_per_patterns
# colnames(genes_per_patterns)
# output <- list()
# for (i in 1:nrow(data)) {
#   diffexpressed <- as.character(data$diffexpressed_alltp[i])
#   print(diffexpressed)
#   gene_list <- strsplit(data$MGI.symbol[i], "; ")[[1]]
#   output[[diffexpressed]] <- gene_list
# }

# # Create individual df for each patterns
# for (diffexpressed in names(output)) {
#   assign(diffexpressed, output[[diffexpressed]])
# }

# pattern_list <- unique(patterns$diffexpressed_alltp)
module_list <- summary.output$modules

total.genes <- vcount(g) # Make this the total number of genes across all modules (size of the biggest module)
length(MEGENA.output$module.output)
modtable

# for(pattern in pattern_list){
#   overlapgenes <- get(pattern)
#   overlaplist <- list()
#   overlaplist[["genes"]] <- overlapgenes
  
#   Object <- newGOM(overlaplist, module_list, total.genes)
  
#   overlap.intersect <- getMatrix(Object, name="intersection")
#   overlap.intersect.list <- getNestedList(Object, name="intersection")
#   overlap.pval <- getMatrix(Object, name="pval")
#   overlap.OR <- getMatrix(Object, name="odds.ratio")
#   Overlap.sum <- cbind(overlap.intersect,overlap.pval,overlap.OR, overlap.intersect.list)
#   Overlap.sum <- data.frame(Overlap.sum)
#   Overlap.sum <- Overlap.sum %>% arrange(overlap.pval)
  
#   assign(paste0("alluvial_enrichment_", pattern), Overlap.sum, envir = .GlobalEnv)
#   # write.csv(Overlap.sum, file=paste0("~/Rcode/megena/enrichment_", 
#   #                                    overlap, "_drn.csv"))
# }
# save(list=ls()[grepl("alluvial_enrichment",ls())],
#      file=paste0(output.path, "enrichment_alluvial/enrichment_alluvial_", reg, ".Rdata"))
# print("Finished alluvial patterns enrichment in modules.")



# #### DEG enrichment in modules ####
# print("Starting DEG enrichment in modules.")
# load(deglist.path)
# tplist <- c("tp1", "tp2", "tp3")
# for(tp in tplist){
#   print(paste0("Processing time point: ", tp))
#   degtp <- get(paste0("deg_", reg, "_", tp))
#   degs <- degtp$label
#   DEGs <- list()
#   DEGs[["DEGs"]] <- degs
  
#   Object <- newGOM(DEGs,module_list,total.genes)
  
#   overlap.intersect <- getMatrix(Object, name="intersection")
#   overlap.intersect.list <- getNestedList(Object, name="intersection")
#   overlap.pval <- getMatrix(Object, name="pval")
#   overlap.OR <- getMatrix(Object, name="odds.ratio")
#   Overlap.sum <- cbind(overlap.intersect,overlap.pval,overlap.OR, overlap.intersect.list)
#   Overlap.sum <- data.frame(Overlap.sum)
#   Overlap.sum <- Overlap.sum %>% arrange(overlap.pval)
  
#   assign(paste0("deg_enrichment_", reg, "_", tp), Overlap.sum, envir = .GlobalEnv)
# }
# save(list=ls()[grepl("deg_enrichment",ls())], 
#      file=paste0(output.path, "enrichment_DEGs/enrichment_degs_", reg, ".Rdata"))
# print("Finished DEG enrichment in modules.")



#### cellular function in modules using mithil association table ####
library(readxl)
brain_cell_markers <- read_excel(brain_cell_markers.path,
                                 sheet = "genes_>5fpkm_enrichment")
print("Loaded brain cell markers.")
for(cell_type in colnames(brain_cell_markers)){
  genelist <- list()
  genelist[["genes"]] <- na.omit(brain_cell_markers[[cell_type]]) # rajouter na.omit au lieu de get 
  length(genelist[[1]])
  head(genelist[[1]])
  
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
print("Finished cell type enrichment in modules.")


### Enrichment of genes on chrX and chrY in modules ###
Object_chrXY <- newGOM(genelist_chrXY, module_list, total.genes)

overlap.intersect_chrXY <- getMatrix(Object_chrXY, name="intersection")
overlap.intersect.list_chrXY <- getNestedList(Object_chrXY, name="intersection")
overlap.pval_chrXY <- getMatrix(Object_chrXY, name="pval")
overlap.OR_chrXY <- getMatrix(Object_chrXY, name="odds.ratio")

Overlap_chrXY <- cbind(
  overlap.intersect_chrXY,
  overlap.pval_chrXY,
  overlap.OR_chrXY,
  overlap.intersect.list_chrXY
)

Overlap_chrXY <- data.frame(Overlap_chrXY)
Overlap_chrXY <- Overlap_chrXY %>% arrange(overlap.pval_chrXY)

save(Overlap_chrXY,
     file = paste0(output.path, "enrichment_chrXY/enrichment_chrXY_", reg, ".Rdata"))
print("Finished chrX and chrY genes enrichment in modules.")

