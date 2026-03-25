#### PREP OF THE FILE USE BY RRHO20 ####
# Modification du fichier counts_final contenant counts + résultats deseq2 (à la fin) 

#### LIBRARIES ####
library(dplyr)

rm(list=ls())

#### PATHS ####
project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

input <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/deseq2_panregion/panregion_betaprior_rg_"
output <- "data/3__RRHO2/panregion/design_reg_group/prep_rrho2/annotated_counts_filtered_panregion_rg_"


list_tp <- as.list(c("tp1", "tp2", "tp3"))
combined_table <- NULL   # list of generated dataframes to create a file containing everything 


for (tp in list_tp) {
  input_file <- readRDS(paste0(input, tp, ".rds"))  #input file with rawcounts and DESeq2 results
  output_tp <- paste0(output, tp, ".rds") # output path
  
  genes_names <- input_file %>% select("MGI.symbol")  # save gene names
  
  # Generate df for each tp 
  table <- input_file %>% select("log2FoldChange", "pvalue", "padj")   # select columns of interest
  colnames(table) <- tolower(c(paste0("log2fc_", tp),           # change colnames
                                  paste0("pval_", tp), 
                                  paste0("padj_", tp)))
  
  table_combined <- cbind(genes_names, table)  # save df for each tp 
  saveRDS(table_combined, output_tp)
  
  # Combine all tp in one df
  if (is.null(combined_table)) {
    combined_table <- table_combined  # initialize 
  } else {
    # merge df 
    combined_table <- full_join(combined_table, table_combined, by = "MGI.symbol")
  }
}

#### SAVE ####
saveRDS(combined_table, paste0(output, "alltp.rds"))
