# install.packages("readODS", repos="https://cloud.r-project.org")
library(tidyverse)
library(edgeR)
library(limma)
library(readODS)

setwd("/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/")


regionlist <- c("ACC", "Hb", "Ins", "Nac")

# PATHS ####
annotated_counts.path <- "data/2__differential_expression_analysis/annotated_counts.csv" # raw counts
coldata.path <- "data/counts_m39_M32/coldata_without_outliers.ods" # pour récpérer les RIN et les groupes
output.path <- "data/4__MEGENA/" # folder to save results
dir.create(output.path)

for(reg in regionlist){
  # load unnormalised counts as dataframe 
  annotated_counts <- read_csv(annotated_counts.path)
  counts <- annotated_counts %>%
    distinct(MGI.symbol, .keep_all = TRUE) %>%
    column_to_rownames("MGI.symbol") %>%
    select(contains(paste0(reg, ".")))
  
  # load coldata
  cldata <- read_ods(coldata.path)
  
  # keep genes with at least 1 cpm in >= 20% of samples in a single project
  cpm <- apply(counts, 2, function(x) (x/sum(x))*1000000)
  counts_filtered <- counts[rowSums(cpm < 1) < (dim(cpm)[2]/5),]
  
  # TMM normalization
  dge <- DGEList(counts = counts_filtered)
  dge <- calcNormFactors(dge, method = "TMM")
  
  # LOG2 CPM 
  logCPM <- cpm(dge, log = TRUE)

  # Correction du RIN (en protegeant effet cuff/sham)
  samples_region <- colnames(logCPM)
  cldata_df <- as.data.frame(cldata)
  rownames(cldata_df) <- cldata_df$sample
  cldata_region <- cldata_df[samples_region, ]
  
  design <- model.matrix(~ group, data = cldata_region)
  
  logCPM_corrected <- removeBatchEffect(logCPM,
                                        covariates = cldata_region$RIN,
                                        design = design)
  
  # save R object of normalized counts
  filename <- paste0("logCPM_RINcorrected_", reg, ".Rdata")
  save(logCPM_corrected, file = paste0(output.path, filename))
}