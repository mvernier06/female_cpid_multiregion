#### DEGs ANNOTATION PER TP ####

rm(list=ls())

#### LIBRARIES ####
library(dplyr)

#### PATHS #### 

project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

## INPUT ##
panregion_tp1.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/deseq2_panregion/panregion_deseq2_rg_tp1.rds" # file with just deseq2 results 
panregion_tp2.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/deseq2_panregion/panregion_deseq2_rg_tp2.rds"
panregion_tp3.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/deseq2_panregion/panregion_deseq2_rg_tp3.rds"

## OUTPUT ##
annotated_genes_tp1.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_tp1.rds"
annotated_genes_tp2.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_tp2.rds"
annotated_genes_tp3.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_tp3.rds"

annotated_genes_alltp.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_alltp.rds"

deg_tp1.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_tp1.rds"
deg_tp2.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg__tp2.rds"
deg_tp3.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_tp3.rds"

deg_alltp.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_alltp.rds"


#### DEGs ANNOTATION ####
## TP1 ##
panregion_tp1 <- readRDS(panregion_tp1.path)

panregion_tp1$diffexpressed <- "ns" 
panregion_tp1$diffexpressed[panregion_tp1$log2FoldChange > 0 & panregion_tp1$pvalue < 0.05] <- "UP"
panregion_tp1$diffexpressed[panregion_tp1$log2FoldChange < 0 & panregion_tp1$pvalue < 0.05] <- "DOWN"

summary_annotation_TP1 <- table(panregion_tp1$diffexpressed) # summary DEG 
print(summary_annotation_TP1)

## TP2 ## 
panregion_tp2 <- readRDS(panregion_tp2.path)

panregion_tp2$diffexpressed <- "ns" 
panregion_tp2$diffexpressed[panregion_tp2$log2FoldChange > 0 & panregion_tp2$pvalue < 0.05] <- "UP"
panregion_tp2$diffexpressed[panregion_tp2$log2FoldChange < 0 & panregion_tp2$pvalue < 0.05] <- "DOWN"

summary_annotation_TP2 <- table(panregion_tp2$diffexpressed) # summary DEG 
print(summary_annotation_TP2)


## TP3 ## 
panregion_tp3 <- readRDS(panregion_tp3.path)

panregion_tp3$diffexpressed <- "ns" 
panregion_tp3$diffexpressed[panregion_tp3$log2FoldChange > 0 & panregion_tp3$pvalue < 0.05] <- "UP"
panregion_tp3$diffexpressed[panregion_tp3$log2FoldChange < 0 & panregion_tp3$pvalue < 0.05] <- "DOWN"

summary_annotation_TP3 <- table(panregion_tp3$diffexpressed) # summary DEG 
print(summary_annotation_TP3)

### FORMATTING COLNAMES ###
colnames(panregion_tp1) <-  c("label", "log2fc", "pval", "padj", "diffexpressed")
colnames(panregion_tp2) <-  c("label", "log2fc", "pval", "padj", "diffexpressed")
colnames(panregion_tp3) <-  c("label", "log2fc", "pval", "padj", "diffexpressed")

### SAVE ###
## WITH PADJ ##
saveRDS(panregion_tp1, file = annotated_genes_tp1.path)
saveRDS(panregion_tp2, file = annotated_genes_tp2.path)
saveRDS(panregion_tp3, file = annotated_genes_tp3.path)

#### KEEP ONLY DEG, REMOVE PADJ ####
panregion_deg_tp1 <- panregion_tp1 %>% filter(diffexpressed != "ns") %>%
  select("label", "log2fc", "pval", "diffexpressed")
panregion_deg_tp2 <- panregion_tp2 %>% filter(diffexpressed != "ns") %>%
  select("label", "log2fc", "pval", "diffexpressed")
panregion_deg_tp3 <- panregion_tp3 %>% filter(diffexpressed != "ns") %>%
  select("label", "log2fc", "pval", "diffexpressed")

### SAVE ###
## WITHOUT PADJ and NS ##
saveRDS(panregion_deg_tp1, file = deg_tp1.path)
saveRDS(panregion_deg_tp2, file = deg_tp2.path)
saveRDS(panregion_deg_tp3, file = deg_tp3.path)


### ALLTP ###
colnames(panregion_tp1) <-  c("label", "log2fc_tp1", "pval_tp1", "padj_tp1", "diffexpressed_tp1")
colnames(panregion_tp2) <-  c("label", "log2fc_tp2", "pval_tp2", "padj_tp2", "diffexpressed_tp2")
colnames(panregion_tp3) <-  c("label", "log2fc_tp3", "pval_tp3", "padj_tp3", "diffexpressed_tp3")

# merge them all
panregion_alltp <- merge(panregion_tp1, panregion_tp2, by = "label", all = TRUE)
panregion_alltp <- merge(panregion_alltp, panregion_tp3, by = "label", all = TRUE)

### SAVE ###
saveRDS(panregion_alltp, file =annotated_genes_alltp.path)

### ALLTP: WITHOUT PADJ NOR NS ###
panregion_deg_tp1 <- panregion_tp1 %>% filter(diffexpressed_tp1 != "ns") %>%
  select("label", "log2fc_tp1", "pval_tp1", "diffexpressed_tp1")
panregion_deg_tp2 <- panregion_tp2 %>% filter(diffexpressed_tp2 != "ns") %>%
  select("label", "log2fc_tp2", "pval_tp2", "diffexpressed_tp2")
panregion_deg_tp3 <- panregion_tp3 %>% filter(diffexpressed_tp3 != "ns") %>%
  select("label", "log2fc_tp3", "pval_tp3", "diffexpressed_tp3")

# merge them all 
deg_alltp <- merge(panregion_deg_tp1, panregion_deg_tp2, by ="label", all = TRUE)
deg_alltp <- merge(deg_alltp, panregion_deg_tp3, by="label", all = TRUE)

### SAVE ###
saveRDS(deg_alltp, file =deg_alltp.path)
