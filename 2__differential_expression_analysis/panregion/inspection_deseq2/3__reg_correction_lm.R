#### APPLY LINEAR MODEL ON RAW COUNTS ####

rm(list=ls())

#### LIBRAIRIES ####
library(tidyverse)
library(DESeq2)
library(variancePartition)
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


#### LINEAR MODEL ####
tp_list <- list(1, 2, 3)
for (tp in tp_list) {
  ## OUTPUT ## 
  data.path <- paste0("data/2__differential_expression_analysis/panregion/inspection_deseq2/lm/lm_matrix_reg_corrected_tp", tp, ".rds")
  plot.path <- paste0("graphs_results/panregion/inspection_deseq2/variance_partition/variance_partition_lm_tp", tp, ".png")
  data_vp.path <- paste0("data/2__differential_expression_analysis/panregion/inspection_deseq2/lm/variance_partition_matrix_tp",tp, ".rds")
  
  print(paste0("tp : ", tp))
  metadata_tp <- metadata[which(metadata$timepoint == tp),] %>% arrange(sample)   
  list <- metadata_tp$sample[which(metadata_tp$timepoint == tp)] 
  
  ## GET CORRESPONDING RAW COUNTS MATRIX ##
  df_rawcounts_tp <- raw_counts_filtered[which(colnames(raw_counts_filtered) %in% list,)]
  
  #### LINEAR MODEL FUCNTION #### 
  correct_expression_data <- function(count_matrix, metadata) {
    
    # Vérifier correspondance samples
    if (!all(colnames(count_matrix) %in% metadata$sample)) {
      stop("Sample IDs in count matrix and metadata don't match")
    }
    
    # Réordonner metadata
    metadata <- metadata[match(colnames(count_matrix), metadata$sample), ]
    
    # DESeq2 object (pour VST uniquement)
    dds <- DESeqDataSetFromMatrix(
      countData = count_matrix,
      colData = metadata,
      design = ~1
    )
    
    # VST transformation
    vst_data <- vst(dds, blind = TRUE)
    vst_matrix <- assay(vst_data)
    
    # Correction batch + covariate en protégeant group
    corrected <- removeBatchEffect(
      vst_matrix,
      batch = metadata$reg,
      covariates = metadata$RIN,
      design = model.matrix(~ group, data = metadata)
    )
    
    return(corrected)
  }
  

  corrected_matrix <- correct_expression_data(df_rawcounts_tp, metadata_tp)
  
  saveRDS(corrected_matrix, data.path)
  
    
  info <- metadata_tp %>%tibble::column_to_rownames(var = "sample")
  
  # Vérifier alignement
  stopifnot(all(colnames(corrected_matrix) == rownames(info)))
  
  # Modèle cohérent avec données corrigées
  formula <- ~ RIN + (1 | reg) + (1 | group)
  
  varPart <- fitExtractVarPartModel(corrected_matrix, formula, info)
  
  #### SAVE ####
  title <- paste0("Variance Partition (corrected reg + RIN) tp", tp)
  
  plot <- plotVarPart(varPart) + 
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(plot.path, plot = plot, bg="white",
         width=1500, height=1000, units="px", scale=2)
  
  saveRDS(varPart, data_vp.path)
  
 
}
