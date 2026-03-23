#### CLUSTER_PROFILER #### 
rm(list=ls())
#### LIBRARIES ####
dir.create("~/R/library", showWarnings = FALSE, recursive = TRUE)
.libPaths(c("~/R/library", .libPaths()))
Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always") 
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(fgsea)
library(stats)
library(ggplot2)

organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)


#### PATHS ####
project.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project.path)
## INPUT ##
# GO #
deg_tp1.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_tp1.rds" 
deg_tp2.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_tp2.rds"
deg_tp3.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_tp3.rds"
genes_filtered.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_alltp.rds"

# GSEA #
gse_genes_tp1.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_tp1.rds"
gse_genes_tp2.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_tp2.rds"
gse_genes_tp3.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_tp3.rds"

## OUTPUT ## 
# GO #
go_res_tp1 <- "data/2__differential_expression_analysis/panregion/go/go_enrich/go_enrich_panregion_tp1.rds"
go_res_tp2 <- "data/2__differential_expression_analysis/panregion/go/go_enrich/go_enrich_panregion_tp2.rds"
go_res_tp3 <- "data/2__differential_expression_analysis/panregion/go/go_enrich/go_enrich_panregion_tp3.rds"

go_plot_tp1 <- "graphs_results/panregion/go/cluster_profiler/go_enrich_panregion_tp1.png"
go_plot_tp2 <- "graphs_results/panregion/go/cluster_profiler/go_enrich_panregion_tp2.png"
go_plot_tp3 <- "graphs_results/panregion/go/cluster_profiler/go_enrich_panregion_tp3.png"

# GSE # 
gse_data_tp1.path <- "data/2__differential_expression_analysis/panregion/go/gse/gse_panregion_tp1.rds"
gse_data_tp2.path <- "data/2__differential_expression_analysis/panregion/go/gse/gse_panregion_tp2.rds"
gse_data_tp3.path <- "data/2__differential_expression_analysis/panregion/go/gse/gse_panregion_tp3.rds"

gse_plot_tp1.path <- "graphs_results/panregion/go/cluster_profiler/gse_panregion_tp1.png"
gse_plot_tp2.path <- "graphs_results/panregion/go/cluster_profiler/gse_panregion_tp2.png"
gse_plot_tp3.path <- "graphs_results/panregion/go/cluster_profiler/gse_panregion_tp3.png"


#### CLUSTER PROFILER : GO ENRICH ####
deg_tp1 <- readRDS(deg_tp1.path)
deg_tp2 <- readRDS(deg_tp2.path)
deg_tp3 <- readRDS(deg_tp3.path)

genes_filtered <- readRDS(genes_filtered.path)
genes <- genes_filtered$label

## TP1 ##
go_panregion_tp1 <- enrichGO(gene = deg_tp1$label, 
                             OrgDb = organism,                # mouse Db
                             keyType = "SYMBOL", 
                             ont ="ALL",                      # all ontologies
                             qvalueCutoff = 0.1,
                             pAdjustMethod = "BH",            # Benjamini-Hochberg (FDR)
                             universe = genes)                # reference genes for analysis

# overview (may be empty = no results)
view(go_panregion_tp1)

barplot(go_panregion_tp1, showCategory = 20)
dotplot(go_panregion_tp1, showCategory = 10)

### SAVE PLOT ###
png(go_plot_tp1, width = 1500, height = 1000)
barplot(go_panregion_tp1, showCategory = 20)
dev.off()

## TP2 ##
go_panregion_tp2 <- enrichGO(gene = deg_tp2$label, 
                             OrgDb = organism, 
                             keyType = "SYMBOL", 
                             ont ="ALL", 
                             qvalueCutoff = 0.1,
                             pAdjustMethod = "BH",
                             universe = genes)

# overview (may be empty = no results)
view(go_panregion_tp2)

# examples of visualization 
barplot(go_panregion_tp2, showCategory = 20)
dotplot(go_panregion_tp2, showCategory = 10)

# sort by ascending order (padj values)
go_panregion_tp2@result <- go_panregion_tp2@result[order(go_panregion_tp2@result$p.adjust), ]

# overview
barplot(go_panregion_tp2, showCategory = 20)


### SAVE PLOT ###
png(go_plot_tp2, width = 1500, height = 1000)
barplot(go_panregion_tp2, showCategory = 20)
dev.off()

## TP3 ##
go_panregion_tp3 <- enrichGO(gene = deg_tp3$label, 
                             OrgDb = organism, 
                             keyType = "SYMBOL", 
                             ont ="ALL", 
                             qvalueCutoff = 0.1,
                             pAdjustMethod = "BH",
                             universe = genes)

# overview (may be empty = no results)
view(go_panregion_tp3)

barplot(go_panregion_tp3, showCategory = 20)
dotplot(go_panregion_tp3, showCategory = 10)

### SAVE PLOT ###
png(go_plot_tp3, width = 1500, height = 1000)
barplot(go_panregion_tp3, showCategory = 20)
dev.off()


### SAVE GO ENRICH DATA ###
saveRDS(go_panregion_tp1, go_res_tp1)
saveRDS(go_panregion_tp2, go_res_tp2)
saveRDS(go_panregion_tp3, go_res_tp3)


#### GSEA #### 
## SELECT DATA ##
gse_tp1 <-  readRDS(gse_genes_tp1.path) %>% dplyr::select(label, log2fc)
gse_tp2 <-  readRDS(gse_genes_tp2.path) %>% dplyr::select(label, log2fc)
gse_tp3 <-  readRDS(gse_genes_tp3.path) %>% dplyr::select(label, log2fc)

## TP1 ##
log2fc <- gse_tp1$log2fc
names(log2fc) <- gse_tp1$label
# descending ranking of log2fc 
genelist <- sort(log2fc, decreasing = TRUE)

# set the seed for reproductibility 
set.seed(123)  # number is arbitrary

## LAUNCH GSEA ##
gse_tp1 <- gseGO(geneList=genelist, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose =  TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "BH",
                 seed = 123)
#eps = 0)

# Visualiser avec cnetplot (test)
# cnetplot(gse_tp1, foldChange = genelist)  # Génère un graphique en réseau des termes GO
# overview
dotplot(gse_tp1, showCategory = 25)
view(gse_tp1)


### SAVE GSE PLOT ### 
## TP1 ##
png(gse_plot_tp1.path, width = 1500, height = 1000)
dotplot(gse_tp1, showCategory = 25)
dev.off()

## TP2 ##
log2fc <- gse_tp2$log2fc
names(log2fc) <- gse_tp2$label
# descending ranking of log2fc 
genelist <- sort(log2fc, decreasing = TRUE)

# set the seed for reproductibility 
set.seed(123)  # number should be the same for each TP

## LAUNCH GSEA ##
gse_tp2 <- gseGO(geneList=genelist, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose =  TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "BH",
                 seed = 123)

# overview
dotplot(gse_tp2, showCategory = 25)
view(gse_tp2)

### SAVE GSE PLOT ### 
## TP2 ## 
png(gse_plot_tp2.path, width = 1500, height = 1000)
dotplot(gse_tp2, showCategory = 25)
dev.off()


## TP3 ##
log2fc <- gse_tp3$log2fc
names(log2fc) <- gse_tp3$label
# descending ranking of log2fc 
genelist <- sort(log2fc, decreasing = TRUE)

# set the seed for reproductibility 
set.seed(123)  # number should be the same for each TP

## LAUNCH GSEA ##
gse_tp3 <- gseGO(geneList=genelist, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose =  TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "BH",
                 seed = 123)

# overview
dotplot(gse_tp3, showCategory = 25)
view(gse_tp3)

### SAVE GSE PLOT ### 
## TP3 ##
png(gse_plot_tp3.path, width = 1500, height = 1000)
dotplot(gse_tp3, showCategory = 25)
dev.off()


### SAVE GSE DATA ###
saveRDS(gse_tp1, gse_data_tp1.path)
saveRDS(gse_tp2, gse_data_tp2.path)
saveRDS(gse_tp3, gse_data_tp3.path)
