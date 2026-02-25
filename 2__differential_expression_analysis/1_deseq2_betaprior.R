# betaprior sans outliers

library(tidyverse)
library(DESeq2)
library(readODS)

rm(list=ls())

setwd("/home/marinevernier/Documents/projets/")

#### PATHS ####
raw_counts.path <- "female_cpid_multiregion/data/counts_m39_M32/cpid_multireg_counts.txt"
coldata.path <- "female_cpid_multiregion/data/counts_m39_M32/coldata.ods"
output.path <- "female_cpid_multiregion/data/2__differential_expression_analysis/CPID_sham_vs_cuff_betaprior.csv" # file to save results

#### Formatting raw counts ####
raw_counts <- read.table(raw_counts.path, header = TRUE,
                         sep = "\t",
                         comment.char = "#",
                         stringsAsFactors = FALSE,
                         check.names = FALSE)

raw_counts <- raw_counts[, !(colnames(raw_counts) %in% c(
  "Chr", "Start", "End", "Strand", "Length"
))]

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

#### Formatting metadata ####
coldata <- read_ods(coldata.path)
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample
coldata$timepoint <- as.factor(coldata$timepoint) 
coldata$reg <- as.factor(coldata$reg)
coldata$sample <- as.factor(coldata$sample)
coldata <- coldata %>% arrange(by = sample)
coldata %>% filter(timepoint == 1) %>% nrow


### Excluding outliers ###
outliers <- c("Ins.1837", "Nac.1837", "Hb.1839", "Hb.2049")
coldata <- coldata %>%
  filter(!sample %in% outliers)
raw_counts <- raw_counts[, coldata$sample, drop = FALSE]

#### Differential expression analysis of sham vs cuff ####
Deseq2MultiReg <- function(regionList, tpList){
  
  for (tp in tpList){
    
    # Select rows of the timepoint in coldata
    metadata <- coldata[which(coldata$timepoint == tp),]
    
    # Select columns of the timepoint in raw_counts
    list <- metadata$sample[which(metadata$timepoint == tp)]
    df <- raw_counts[which(colnames(raw_counts) %in% list,)]
    
    for (region in regionList){
      
      print(paste0("Computing stastistics about TP", tp, " of ", region))
      
      # Compute pval, qval and lfc related to the difference of sham vs cuff for each timepoint per region
      dds <- DESeqDataSetFromMatrix(countData=df[, grepl(region, names(df))], 
                                    colData=metadata[str_detect(metadata$reg, region),], 
                                    design=~ RIN + group)
      
      dds <- DESeq(object = dds, betaPrior = TRUE)
      
      res_tp <- results(dds, contrast=c("group", "cuff", "sham"), 
                        cooksCutoff=TRUE, independentFiltering=TRUE, alpha=0.05)
      
      # Rearrange data and save results at the end of raw_counts
      res_tp <-  as.data.frame(res_tp[, c('log2FoldChange', 'pvalue', 'padj'), drop = FALSE])
      colnames(res_tp) <- c(
        paste0(region, "_log2fc_tp", tp),
        paste0(region, "_pval_tp", tp),
        paste0(region, "_padj_tp", tp)
      )
      
      raw_counts <- cbind(raw_counts, res_tp)
    }
  }
  raw_counts_out <- cbind(Geneid = rownames(raw_counts), raw_counts)
  write.csv(raw_counts_out, output.path, row.names = FALSE)

}

#### Call to the function ####
Deseq2MultiReg(unique(coldata$reg), unique(coldata$timepoint))

all(colnames(df) == rownames(metadata))
