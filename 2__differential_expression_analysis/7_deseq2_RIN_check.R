# test deseq2 avec et sans RIN 

library(tidyverse)
library(DESeq2)
library(readODS)

rm(list=ls())

setwd("/home/marinevernier/Documents/projets/")

#### PATHS ####
raw_counts.path <- "female_cpid_multiregion/data/counts_m39_M32/cpid_multireg_counts.txt"
coldata.path <- "female_cpid_multiregion/data/counts_m39_M32/coldata.ods"
annot_table.path <- "female_cpid_multiregion/data/counts_m39_M32/annotation_final.csv"

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

#### Differential expression analysis of sham vs cuff  with RIN ####
Deseq2MultiReg_RIN <- function(regionList, tpList){
  
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
  
}

#### Differential expression analysis of sham vs cuff  with RIN ####
Deseq2MultiReg_basic <- function(regionList, tpList){
  
  for (tp in tpList){
    
    # Select rows of the timepoint in coldata
    metadata <- coldata[which(coldata$timepoint == tp),]
    
    # Select columns of the timepoint in raw_counts
    list <- metadata$sample[which(metadata$timepoint == tp)]
    df <- raw_counts[which(colnames(raw_counts) %in% list,)]
    
    for (region in regionList){
      
      print(paste0("Computing stastistics about TP", tp, " of ", region))
      counts$Geneid <-  str_replace(counts$Geneid, "\\..*", "") 
      
      ens2symbol <- read.csv(annot_table.path, sep=";")
      ens2symbol <- ens2symbol %>% dplyr::select(!X)
      
      annotated_counts <- inner_join(counts, ens2symbol, by=c("Geneid"="Gene.stable.ID"))
      annotated_counts <- annotated_counts %>% relocate(MGI.symbol, .before=Geneid) 
      
      # remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)
      annotated_counts <- annotated_counts[!duplicated(annotated_counts$MGI.symbol),]
      # Compute pval, qval and lfc related to the difference of sham vs cuff for each timepoint per region
      dds <- DESeqDataSetFromMatrix(countData=df[, grepl(region, names(df))], 
                                    colData=metadata[str_detect(metadata$reg, region),], 
                                    design=~ group)
      
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
  
}

#### Call to the function ####
counts_RIN <- Deseq2MultiReg_RIN(unique(coldata$reg), unique(coldata$timepoint))
counts_basic <-  Deseq2MultiReg_basic(unique(coldata$reg), unique(coldata$timepoint))

counts_RIN$Geneid <-  str_replace(counts_RIN$Geneid, "\\..*", "") 

ens2symbol <- read.csv(annot_table.path, sep=";")
ens2symbol <- ens2symbol %>% dplyr::select(!X)

counts_RIN <- inner_join(counts_RIN, ens2symbol, by=c("Geneid"="Gene.stable.ID"))
counts_RIN <- counts_RIN %>% relocate(MGI.symbol, .before=Geneid) 

# remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)
counts_RIN <- counts_RIN[!duplicated(counts_RIN$MGI.symbol),]

counts_basic$Geneid <-  str_replace(counts_basic$Geneid, "\\..*", "") 

counts_basic <- inner_join(counts_basic, ens2symbol, by=c("Geneid"="Gene.stable.ID"))
counts_basic <- counts_basic %>% relocate(MGI.symbol, .before=Geneid) 

# remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)
counts_basic <- counts_basic[!duplicated(counts_basic$MGI.symbol),]


# Colonnes log2FC
l2fc_basic_cols <- grep("_log2fc_tp", colnames(counts_basic), value = TRUE)
l2fc_RIN_cols   <- grep("_log2fc_tp", colnames(counts_RIN), value = TRUE)

stopifnot(l2fc_basic_cols == l2fc_RIN_cols)

x_all <- unlist(counts_basic[, l2fc_basic_cols])
y_all <- unlist(counts_RIN[, l2fc_RIN_cols])

keep <- complete.cases(x_all, y_all)

cor_global <- cor.test(x_all[keep], y_all[keep])

plot(x_all[keep], y_all[keep],
     pch = 16,
     cex = 0.4,
     xlab = "log2FC sans RIN",
     ylab = "log2FC avec RIN",
     main = paste0("Correlation globale r = ",
                   round(cor_global$estimate, 3)))

abline(0,1,col="red", lwd=2)

plot <- list()
for (i in seq_along(l2fc_basic_cols)) {
  
  x <- counts_basic[[ l2fc_basic_cols[i] ]]
  y <- counts_RIN[[ l2fc_RIN_cols[i] ]]
  
  keep <- complete.cases(x, y)
  
  r <- cor(x[keep], y[keep])
  
  p <- plot(x[keep], y[keep],
       pch = 16,
       cex = 0.4,
       main = paste0(l2fc_basic_cols[i],
                     "\nr = ", round(r,3)),
       xlab = "sans RIN",
       ylab = "avec RIN")
  
  abline(0,1,col="red", lwd=2)
  plot[[i]] <- p
}


delta_matrix <- counts_RIN[, l2fc_RIN_cols] - 
  counts_basic[, l2fc_basic_cols]

rownames(delta_matrix) <- counts_basic$MGI.symbol

mean_delta <- rowMeans(delta_matrix, na.rm = TRUE)

result_delta <- data.frame(
  Geneid = rownames(delta_matrix),
  mean_delta_LFC = mean_delta,
  abs_mean_delta_LFC = abs(mean_delta)
)
result_delta <- result_delta[order(-result_delta$abs_mean_delta_LFC), ]

head(result_delta, 20)

## Facteur confondant : 
boxplot(RIN ~ group, data = coldata,
        main = "Distribution du RIN par groupe")
wilcox.test(RIN ~ group, data = coldata)


boxplot(RIN ~ group + reg, data = coldata)

interaction.plot(coldata$group,
                 coldata$reg,
                 coldata$RIN)

anova(lm(RIN ~ group + reg + timepoint, data = coldata))
