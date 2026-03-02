library(tidyverse)
library(edgeR)

setwd("/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/")

regionlist <- c("ACC", "Hb", "Ins", "Nac")

# PATHS ####
annotated_counts.path <- "data/2__differential_expression_analysis/annotated_counts.csv" # raw counts
output.path <- "data/4__MEGENA/" # folder to save results
dir.create(output.path)

for(reg in regionlist){
  # load unnormalised counts as dataframe 
  annotated_counts <- read_csv(annotated_counts.path)
  counts <- annotated_counts %>%
    distinct(MGI.symbol, .keep_all = TRUE) %>%
    column_to_rownames("MGI.symbol") %>%
    select(contains(paste0(reg, ".")))
  
  # keep genes with at least 1 cpm in >= 20% of samples in a single project
  cpm <- apply(counts, 2, function(x) (x/sum(x))*1000000)
  counts_filtered <- counts[rowSums(cpm < 1) < (dim(cpm)[2]/5),]
  
  # CTF normalization
  lib_size <- base::colSums(counts_filtered)
  norm_factors <- calcNormFactors(object = counts_filtered, lib.size = lib_size, method = "TMM")
  CTF_normalized <- sweep(counts_filtered, 2, norm_factors, "/")
  
  # save R object of normalized counts
  filename <- paste0("CTF_normalized_counts_", reg, ".Rdata")
  save(CTF_normalized, file=paste0(output.path, filename))
}