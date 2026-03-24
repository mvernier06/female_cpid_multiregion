rm(list=ls())

#### LIBRARIES ####
library(tidyverse)
library(ggplot2)

#### PATHS ####
project_path <- "/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/"
setwd(project_path)


overlap_deg_cell <- function(reg, tplist){
  print(paste0("Running overlap between DEGs and cell type markers for region: ", reg))
  deg_enrichment_reg.path <- paste0("data/4__MEGENA/MEGENA_without_RIN_correction/enrichment_DEGs/enrichment_degs_", reg, ".Rdata")
  enrichement_reg_cell_type.path <- paste0("data/4__MEGENA/MEGENA_without_RIN_correction/enrichment_cell_types/enrichment_cell_types_", reg, ".Rdata")
  
  load(deg_enrichment_reg.path)
  load(enrichement_reg_cell_type.path)
  
  deg_modules <- list()
  
  for(tp in tplist){
    
    ##recupère les module enrichis en DEG
    df <- get(paste0("deg_enrichment_", reg, "_", tp)) %>%
      rownames_to_column("module") %>%
      mutate(
        module = sub("\\..*", "", module),   # garder seulement "module"
        overlap.pval = as.numeric(overlap.pval),
        overlap.OR = as.numeric(overlap.OR)
      ) %>%
      filter(overlap.pval < 0.05) %>%
      mutate(timepoint = tp)
    
    deg_modules[[tp]] <- df
  }
  
  deg_modules_df <- bind_rows(deg_modules)
  
  
  celltypes <- c("Astrocytes_more_than_5_fpkm",
                 "Neuron_more_than_5_fpkm",
                 "Myelinating.Oligodendrocytes_more_than_5_fpkm",
                 "Microglia_more_than_5_fpkm",
                 "Endothelial.Cells_more_than_5_fpkm")
  
  cell_modules <- list()
  
  ## recupère les modules enrichis en cell type 
  for(ct in celltypes){
    
    df <- get(ct) %>%
      rownames_to_column("module") %>%
      mutate(
        module = sub("\\..*", "", module),
        overlap.pval = as.numeric(overlap.pval),
        overlap.OR = as.numeric(overlap.OR)
      ) %>%
      filter(overlap.pval < 0.05) %>%
      mutate(cell_type = ct)
    
    cell_modules[[ct]] <- df
  }
  
  cell_modules_df <- bind_rows(cell_modules)
  
  
  
  deg_cell_overlap <- inner_join(
    deg_modules_df,
    cell_modules_df,
    by = "module",
    suffix = c("_DEG", "_cell")
  )
  
  deg_cell_overlap %>%
    select(module, timepoint, cell_type,
           overlap.pval_DEG,
           overlap.pval_cell,
           overlap.OR_DEG,
           overlap.OR_cell)
  
  head(deg_cell_overlap)
  
  return(deg_cell_overlap)
}

module_enriched_cell_deg_ACC <- overlap_deg_cell("ACC", c("tp1","tp2","tp3"))
save(module_enriched_cell_deg_ACC, file="data/4__MEGENA/MEGENA_without_RIN_correction/enrichment_cell_types_DEGs/ACC_enrichment_cell_types_deg.Rdata")
module_enriched_cell_deg_Ins <- overlap_deg_cell("Ins", c("tp1","tp2","tp3"))
save(module_enriched_cell_deg_Ins, file="data/4__MEGENA/MEGENA_without_RIN_correction/enrichment_cell_types_DEGs/Ins_enrichment_cell_types_deg.Rdata")
module_enriched_cell_deg_Hb <- overlap_deg_cell("Hb", c("tp1","tp2","tp3"))
save(module_enriched_cell_deg_Hb, file="data/4__MEGENA/MEGENA_without_RIN_correction/enrichment_cell_types_DEGs/Hb_enrichment_cell_types_deg.Rdata")
module_enriched_cell_deg_Nac <- overlap_deg_cell("Nac", c("tp1","tp2","tp3"))
save(module_enriched_cell_deg_Nac, file="data/4__MEGENA/MEGENA_without_RIN_correction/enrichment_cell_types_DEGs/Nac_enrichment_cell_types_deg.Rdata")

print("Finished overlap between DEGs and cell type markers for all regions.")