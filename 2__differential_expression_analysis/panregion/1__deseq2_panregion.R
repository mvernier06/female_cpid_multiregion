#### DESeq2 - design ~ reg + RIN + group #### 

rm(list=ls())

#### LIRBAIRIES ####
library(DESeq2)
library(tidyverse)
library(dplyr)

project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

#### PATHS #### 
## INPUT ##
raw_counts.path <- "data/counts_m39_M32/cpid_multireg_counts.txt"
coldata.path <- "data/counts_m39_M32/coldata.ods"
annot_table.path <- "data/counts_m39_M32/annotation_final.csv"
filtered_counts.path <- "data/2__differential_expression_analysis/annotated_counts_filtered.rds"

#### RAW COUNTS MODIFICATIONS ####  
raw_counts <- read.table(raw_counts.path, header = TRUE,
                         sep = "\t",
                         comment.char = "#",
                         stringsAsFactors = FALSE,
                         check.names = FALSE)
raw_counts <- raw_counts[, !(colnames(raw_counts) %in% c(
  "Chr", "Start", "End", "Strand", "Length"
))]
raw_counts$Geneid <-  str_replace(raw_counts$Geneid, "\\..*", "")
rownames(raw_counts) <- raw_counts$Geneid
raw_counts$Geneid <- NULL
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


#### COLDATA MODIFICATIONS #### 
coldata <- read_ods(coldata.path)
rownames(coldata) <- coldata$sample


## FACTOR ##
coldata$timepoint <- as.factor(coldata$timepoint) 
coldata$group <-  as.factor(coldata$group)
coldata$reg <- as.factor(coldata$reg)
coldata$sample <- as.factor(coldata$sample)

coldata <- coldata %>% arrange(by = sample) # order rows by the values of "sample"

### Excluding outliers ###
outliers <- c("Ins.1837", "Nac.1837", "Hb.1839", "Hb.2049")
coldata <- coldata %>%
  filter(!sample %in% outliers)
raw_counts <- raw_counts[, coldata$sample, drop = FALSE]

## Load annotation table and filtered genes
ens2symbol <- read.csv(annot_table.path, sep=";")
ens2symbol <- ens2symbol %>% dplyr::select(!X)
filtered_counts <- readRDS(filtered_counts.path)

for (tp in c("1", "2", "3")){
  ### TP SELECTION IN COLDATA ### 
  
  ### SELECT THE TIMEPOINTS DATA ###  
  metadata <- coldata[which(coldata$timepoint == tp),]        # get metadata for tp 
  list <- metadata$sample[which(metadata$timepoint == tp)]    # get list of samples from tp 
  df_counts <- raw_counts[which(colnames(raw_counts) %in% list,)] # get rawcounts only with 1tp's samples
  
  ### Center RIN to avoid collinearity with intercept ###
  metadata$RIN_c <- scale(metadata$RIN, center = TRUE, scale = FALSE)
  
  #### RUN DESEQ2 ####
  dds <- DESeqDataSetFromMatrix(countData=df_counts,   # creation of the DESeq2 object
                                colData= metadata, 
                                design=~reg+RIN_c+group)
  
  dds <- DESeq(object = dds, betaPrior = TRUE)  # launch DESeq2 with shrinkage (using betaPrior)
  
  normalized_counts <- counts(dds, normalized = TRUE) # normalized counts matrix (optional)
  
  res_tp <- results(dds, contrast=c("group", "cuff", "sham"), # sham is the control 
                    cooksCutoff=TRUE,                             
                    independentFiltering=TRUE, 
                    alpha=0.05) # significance threshold
  
  res_tp <-  as.data.frame(res_tp[, c('log2FoldChange', 'pvalue', 'padj'), 
                                  drop = FALSE]) # get results as a dataframe
  
  df_counts <- df_counts %>% rownames_to_column(var = "X") # Ensembl ID in columns instead of rownames
  raw_counts_deseq <- cbind(df_counts, res_tp) # binds rawcounts with DESeq2 results
  
  #### ADD GENE NAMES #### 
  
  # Join MGI.symbol with Ensembl ID 
 
  annotated_counts <- inner_join(raw_counts_deseq, ens2symbol, 
                                 by=c("X"="Gene.stable.ID"))
  annotated_counts <- annotated_counts %>% relocate(MGI.symbol, .before=X)
  annotated_counts <- annotated_counts[!duplicated(annotated_counts$MGI.symbol),] # on passe de 57010 à 56921 gènes 
  
  ## KEEP ONLY THE 14112 GENES FROM THE FILTERED GENES LIST ## 
  
  counts_final <- annotated_counts %>% 
    filter(MGI.symbol %in% filtered_counts$MGI.symbol) # file with counts and DESeq2 results
  deseq_panregion_tp <- counts_final %>% 
    select(MGI.symbol, log2FoldChange, pvalue, padj) # file with DESeq2 results 
  
  ## NORMALIZED COUNTS : ADD GENES NAMES ##
  normalized_counts <- as.data.frame(normalized_counts)
  normalized_counts <- normalized_counts %>% rownames_to_column(var="X")
  
  normalized_counts <- inner_join(normalized_counts, ens2symbol, by=c("X" = "Gene.stable.ID"))
  normalized_counts <- normalized_counts %>% relocate(MGI.symbol, .before=X)
  normalized_counts <- normalized_counts[!duplicated(normalized_counts$MGI.symbol),]
  
  # Keep only the 14112 genes
  normalized_counts <- normalized_counts %>% 
    filter(MGI.symbol %in% filtered_counts$MGI.symbol)
  
  
  betaprior.path <- paste0("data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/deseq2_panregion/panregion_betaprior_rg_tp",tp,".rds")                 
  deseq2_results.path <- paste0("data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/deseq2_panregion/panregion_deseq2_rg_tp",tp,".rds")       
  normalized_counts.path <- paste0("data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/deseq2_panregion/panregion_normalized_counts_rg_tp",tp,".rds") 
  
  #### SAVE FILES #### 
  saveRDS(counts_final, file = betaprior.path)   #file with counts + deseq2
  saveRDS(deseq_panregion_tp, file = deseq2_results.path)  #file with deseq2 results 
  saveRDS(normalized_counts, file = normalized_counts.path) #file with normalized counts (by size factor only) 
}




