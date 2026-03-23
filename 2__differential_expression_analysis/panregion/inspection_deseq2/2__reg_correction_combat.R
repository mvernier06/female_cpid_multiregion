####  VARIANCE PARTITION  &  COMBAT_SEQ   (LOOPED)  #### 

rm(list=ls())

#### LIBRAIRIES #### 
library(tidyverse)
library(dplyr)
library(variancePartition)
library(sva)
library(readODS)
library(reformulas)

#### PATHS ####
project.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project.path)
## INPUT ##
raw_counts.path <- "data/counts_m39_M32/cpid_multireg_counts.txt"
metadata.path <- "data/counts_m39_M32/coldata_without_outliers.ods"
annot_table.path <- "data/counts_m39_M32/annotation_final.csv"
filtered_counts.path <- "data/2__differential_expression_analysis/annotated_counts_filtered.rds"

## OUTPUT  ## 
output_vp.path <- "graphs_results/panregion/inspection_deseq2/variance_partition/"
output_combat.path <- "data/2__differential_expression_analysis/panregion/inspection_deseq2/combat_seq/"

#### TP LIST #### 
tp_list <- list(1, 2, 3)

#### MODIFICATIONS RAW COUNTS DATAFRAME #### 
raw_counts <- read.table(raw_counts.path, header = TRUE,
                         sep = "\t",
                         comment.char = "#",
                         stringsAsFactors = FALSE,
                         check.names = FALSE)
raw_counts <- raw_counts[, !(colnames(raw_counts) %in% c(
  "Chr", "Start", "End", "Strand", "Length"
))]
raw_counts$Geneid <- str_replace(raw_counts$Geneid, "\\..*", "")
rownames(raw_counts) <- raw_counts$Geneid
colnames(raw_counts) <- sub(
  "_R2\\.dedup\\.bam$",
  "",
  basename(colnames(raw_counts))
)
colnames(raw_counts) <- sub(
  "^([A-Za-z]+)([0-9]+)$",
  "\\1.\\2",
  colnames(raw_counts)
)
colnames(raw_counts) <- gsub("-", ".", colnames(raw_counts))


#### MODIFICATIONS METADATA ####
metadata <- read_ods(metadata.path)
rownames(metadata) <- metadata$sample

## Excluding outliers 
raw_counts <- raw_counts[, metadata$sample, drop = FALSE]
raw_counts$Geneid <- rownames(raw_counts)

#### GENES FILTERING ####
## ADD GENE NAMES ##
ens2symbol <- read.csv(annot_table.path, sep=";")
ens2symbol <- ens2symbol %>% dplyr::select(!X)
annotated_counts <- inner_join(raw_counts, ens2symbol, by=c("Geneid"="Gene.stable.ID"))
annotated_counts <- annotated_counts %>% relocate(MGI.symbol, .before=Geneid)
raw_counts_filtered <- annotated_counts[!duplicated(annotated_counts$MGI.symbol),]

## KEEP ONLY THE 14112 GENES FROM THE FILTERED GENES LIST ##  
filtered_counts <- readRDS(filtered_counts.path)
raw_counts_filtered <- raw_counts_filtered %>% filter(MGI.symbol %in% filtered_counts$MGI.symbol)

## FINAL DATAFRAME ## 
rownames(raw_counts_filtered) <- raw_counts_filtered$MGI.symbol
colonnes_a_exclure <- c("MGI.symbol", "Geneid")
genes_list <- raw_counts_filtered %>% select(all_of(colonnes_a_exclure))
raw_counts_filtered <- raw_counts_filtered %>% select(-all_of(colonnes_a_exclure))

#### LOOP : VARIANCE PARTITION PLOTS & COMBAT_SEQ CORRECTIONS ####
for (tp in tp_list) {
  print(paste("tp:", tp))
  
  ## SELECT CHOSEN TP CORRELECTED WITH METADATA ## 
  print("formatting metadata")
  metadata_tp <- metadata[which(metadata$timepoint == tp),] %>% arrange(sample)   
  list <- metadata_tp$sample[which(metadata_tp$timepoint == tp)] 
  rownames(metadata_tp) <- metadata_tp$sample
  
  ## GET CORRESPONDING RAW COUNTS MATRIX ## 
  print("raw count matrix")
  df_rawcounts_tp <- raw_counts_filtered[which(colnames(raw_counts_filtered) %in% list,)]
  
  all(colnames(df_rawcounts_tp) == rownames(metadata_tp))
  
  ##### VARIANCE PARTITION ##### 
  print("variance partition")
  formula_tp <- ~ RIN + (1 | reg) + (1 | group) 
  varPart_tp <- fitExtractVarPartModel(df_rawcounts_tp, formula_tp, metadata_tp)
  print("plot")
  
  # SAVE VARIATION PARTITION DATA #
  output_vp <- paste0(output_combat.path, "variance_partition_data_tp", tp, ".rds")
  saveRDS(varPart_tp, output_vp)
  
  ## SAVE VARIATION PARTITION PLOT ## 
  plot_tp <- plotVarPart(varPart_tp) + ggtitle(paste0("Variance Partition - TP", tp))
  output_plot_tp <- paste0(output_vp.path, "variance_partition_tp", tp, ".png")
  ggsave(output_plot_tp, plot = plot_tp, bg="white", width=1500, height=1000, units="px", scale=2)
  
  
  ##### COMBAT_SEQ ####
  print("combat correction")
  reg_batch <- as.factor(metadata_tp$reg)
  group <-  as.factor(metadata_tp$group)
  mat <- df_rawcounts_tp %>% as.matrix
  covar_mod <- data.frame(RIN = metadata_tp$RIN)
  
  ## RUN COMBAT_SEQ ##
  adjusted_counts_reg <- ComBat_seq(
    counts = mat,
    batch = reg_batch,
    group = group, 
    covar_mod = covar_mod
  )
  
  df_combat <- as.data.frame(adjusted_counts_reg)
  df_combat_merged <- bind_cols(genes_list, df_combat)
  
  ## SAVE COMBAT_SEQ DATA ##
  print("save combat data")
  saveRDS(df_combat_merged, paste0(output_combat.path, "adjusted_counts_reg_corrected_tp", tp, ".rds"))
  
  
  #### VARIATION PARTITION : CORRECTION CONTROL #### 
  formula_combat_tp <- ~ RIN + (1 | reg) + (1 | group) 
  print("variance partition : corrected version")
  varPart_combat_tp <- fitExtractVarPartModel(df_combat, formula_combat_tp, metadata_tp)
  plot_combat_tp <- plotVarPart(varPart_combat_tp) + ggtitle(paste0("Variance Partition (Reg corrected) - TP", tp))
  
  ## SAVE VARIATION PARTITION DATA ##
  output_vp <- paste0(output_combat.path, "variance_partition_data_corrected_tp", tp, ".rds")
  saveRDS(varPart_combat_tp, output_vp)
  
  ## SAVE VARIATION PARTITION PLOT ##
  output_plot_corrected <- paste0(output_vp.path, "variance_partition_reg_corrected_tp", tp, ".png")
  ggsave(output_plot_corrected, plot = plot_combat_tp, bg="white", width=1500, height=1000, units="px", scale=2)
  
}
