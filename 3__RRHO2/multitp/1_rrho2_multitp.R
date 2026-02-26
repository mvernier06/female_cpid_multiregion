library(tidyverse)
library(RRHO2)

rm(list=ls())
setwd("/home/marinevernier/Documents/projets/female_cpid_multiregion/")

# PATHS ####
annotated_counts.path <- "data/2__differential_expression_analysis/annotated_counts_filtered.rds"
output.path <- "data/3__RRHO2/"
# dir.create(output.path, showWarnings = FALSE, recursive = TRUE)

# Import counts
annotated_counts <- readRDS(annotated_counts.path)

# ACC ####
print("Computing RRHO2 objects of the ACC")
# Filter out NA values
counts_noNA <- annotated_counts %>% dplyr::filter(!is.na(ACC_log2fc_tp1) & 
                                                    !is.na(ACC_pval_tp1) &
                                                    !is.na(ACC_log2fc_tp2) & 
                                                    !is.na(ACC_pval_tp2) &
                                                    !is.na(ACC_log2fc_tp3) & 
                                                    !is.na(ACC_pval_tp3) )

resACC <- annotated_counts %>% 
  select(contains("ACC_") & !contains("padj")) %>%
  na.omit()

# Genelist of TP1.
list1_DDE <- c(-log10(counts_noNA[counts_noNA$ACC_log2fc_tp1 < 0,]$ACC_pval_tp1) * (-1), 
               -log10(counts_noNA[counts_noNA$ACC_log2fc_tp1 > 0,]$ACC_pval_tp1))
gene_list1 <- data.frame(Genes=c(counts_noNA[counts_noNA$ACC_log2fc_tp1 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$ACC_log2fc_tp1 > 0,]$MGI.symbol),
                         DDE = list1_DDE,
                         stringsAsFactors = FALSE)
# Genelist of TP2.
list2_DDE <- c(-log10(counts_noNA[counts_noNA$ACC_log2fc_tp2 < 0,]$ACC_pval_tp2) * (-1),
               -log10(counts_noNA[counts_noNA$ACC_log2fc_tp2 > 0,]$ACC_pval_tp2)) 
gene_list2 <- data.frame(Genes=c(counts_noNA[counts_noNA$ACC_log2fc_tp2 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$ACC_log2fc_tp2> 0,]$MGI.symbol),
                         DDE = list2_DDE,
                         stringsAsFactors = FALSE)
# Genelist of TP3.
list3_DDE <- c(-log10(counts_noNA[counts_noNA$ACC_log2fc_tp3 < 0,]$ACC_pval_tp3) * (-1),
               -log10(counts_noNA[counts_noNA$ACC_log2fc_tp3 > 0,]$ACC_pval_tp3)) 
gene_list3 <- data.frame(Genes=c(counts_noNA[counts_noNA$ACC_log2fc_tp3 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$ACC_log2fc_tp3> 0,]$MGI.symbol),
                         DDE = list3_DDE,
                         stringsAsFactors = FALSE)

RRHO_ACC_1vs2 <- RRHO2_initialize(gene_list1, gene_list2, 
                                  labels = c("Tp1", "Tp2"), 
                                  log10.ind = TRUE)
# Activate method = "fischer" for odds ratio heatmap.
RRHO_ACC_2vs3 <- RRHO2_initialize(gene_list2, gene_list3, 
                                  labels = c("Tp2", "Tp3"), 
                                  log10.ind = TRUE)
RRHO_ACC_1vs3 <- RRHO2_initialize(gene_list1, gene_list3,
                                  labels = c("Tp1", "Tp3"), 
                                  log10.ind = TRUE)
 

# Hb ####
print("Computing RRHO2 objects of the Hb")
counts_noNA <- annotated_counts %>% filter(!is.na(Hb_log2fc_tp1) & 
                                             !is.na(Hb_pval_tp1) &
                                             !is.na(Hb_log2fc_tp2) & 
                                             !is.na(Hb_pval_tp2) &
                                             !is.na(Hb_log2fc_tp3) & 
                                             !is.na(Hb_pval_tp3) )
# Genelist of TP1.
list1_DDE <- c(-log10(counts_noNA[counts_noNA$Hb_log2fc_tp1 < 0,]$Hb_pval_tp1) * (-1), 
               -log10(counts_noNA[counts_noNA$Hb_log2fc_tp1 >= 0,]$Hb_pval_tp1)) 

gene_list1 <- data.frame(Genes=c(counts_noNA[counts_noNA$Hb_log2fc_tp1 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$Hb_log2fc_tp1 >= 0,]$MGI.symbol),
                         DDE = list1_DDE,
                         stringsAsFactors = FALSE)
# Genelist of TP2.
list2_DDE <- c(-log10(counts_noNA[counts_noNA$Hb_log2fc_tp2 < 0,]$Hb_pval_tp2) * (-1),
               -log10(counts_noNA[counts_noNA$Hb_log2fc_tp2 >= 0,]$Hb_pval_tp2)) 
gene_list2 <- data.frame(Genes=c(counts_noNA[counts_noNA$Hb_log2fc_tp2 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$Hb_log2fc_tp2 >= 0,]$MGI.symbol),
                         DDE = list2_DDE,
                         stringsAsFactors = FALSE)
# Genelist of TP3.
list3_DDE <- c(-log10(counts_noNA[counts_noNA$Hb_log2fc_tp3 < 0,]$Hb_pval_tp3) * (-1),
               -log10(counts_noNA[counts_noNA$Hb_log2fc_tp3 >= 0,]$Hb_pval_tp3)) 
gene_list3 <- data.frame(Genes=c(counts_noNA[counts_noNA$Hb_log2fc_tp3 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$Hb_log2fc_tp3 >= 0,]$MGI.symbol),
                         DDE = list3_DDE,
                         stringsAsFactors = FALSE)

RRHO_Hb_1vs2 <- RRHO2_initialize(gene_list1, gene_list2, 
                                 labels = c("Tp1", "Tp2"), 
                                 log10.ind = TRUE)
RRHO_Hb_2vs3 <- RRHO2_initialize(gene_list2, gene_list3, 
                                 labels = c("Tp2", "Tp3"), 
                                 log10.ind = TRUE)
RRHO_Hb_1vs3 <- RRHO2_initialize(gene_list1, gene_list3,
                                 labels = c("Tp1", "Tp3"), 
                                 log10.ind = TRUE)

# Ins ####
print("Computing RRHO2 objects of the Ins")
counts_noNA <- annotated_counts %>% filter(!is.na(Ins_log2fc_tp1) & 
                                             !is.na(Ins_pval_tp1) &
                                             !is.na(Ins_log2fc_tp2) & 
                                             !is.na(Ins_pval_tp2) &
                                             !is.na(Ins_log2fc_tp3) & 
                                             !is.na(Ins_pval_tp3) )
# Genelist of TP1.
list1_DDE <- c(-log10(counts_noNA[counts_noNA$Ins_log2fc_tp1 < 0,]$Ins_pval_tp1) * (-1), 
               -log10(counts_noNA[counts_noNA$Ins_log2fc_tp1 >= 0,]$Ins_pval_tp1))
gene_list1 <- data.frame(Genes=c(counts_noNA[counts_noNA$Ins_log2fc_tp1 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$Ins_log2fc_tp1 >= 0,]$MGI.symbol),
                         DDE = list1_DDE,
                         stringsAsFactors = FALSE)
# Genelist of TP2.
list2_DDE <- c(-log10(counts_noNA[counts_noNA$Ins_log2fc_tp2 < 0,]$Ins_pval_tp2) * (-1),
               -log10(counts_noNA[counts_noNA$Ins_log2fc_tp2 >= 0,]$Ins_pval_tp2)) 
gene_list2 <- data.frame(Genes=c(counts_noNA[counts_noNA$Ins_log2fc_tp2 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$Ins_log2fc_tp2 >= 0,]$MGI.symbol),
                         DDE = list2_DDE,
                         stringsAsFactors = FALSE)
# Genelist of TP3.
list3_DDE <- c(-log10(counts_noNA[counts_noNA$Ins_log2fc_tp3 < 0,]$Ins_pval_tp3) * (-1),
               -log10(counts_noNA[counts_noNA$Ins_log2fc_tp3 >= 0,]$Ins_pval_tp3)) 
gene_list3 <- data.frame(Genes=c(counts_noNA[counts_noNA$Ins_log2fc_tp3 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$Ins_log2fc_tp3 >= 0,]$MGI.symbol),
                         DDE = list3_DDE,
                         stringsAsFactors = FALSE)

RRHO_Ins_1vs2 <- RRHO2_initialize(gene_list1, gene_list2, 
                                  labels = c("Tp1", "Tp2"), 
                                  log10.ind = TRUE)
RRHO_Ins_2vs3 <- RRHO2_initialize(gene_list2, gene_list3, 
                                  labels = c("Tp2", "Tp3"), 
                                  log10.ind = TRUE)
RRHO_Ins_1vs3 <- RRHO2_initialize(gene_list1, gene_list3,
                                  labels = c("Tp1", "Tp3"), 
                                  log10.ind = TRUE)

# Nac ####
print("Computing RRHO2 objects of the Nac")
counts_noNA <- annotated_counts %>% filter(!is.na(Nac_log2fc_tp1) & 
                                             !is.na(Nac_pval_tp1) &
                                             !is.na(Nac_log2fc_tp2) & 
                                             !is.na(Nac_pval_tp2) &
                                             !is.na(Nac_log2fc_tp3) & 
                                             !is.na(Nac_pval_tp3) )
# Genelist of TP1.
list1_DDE <- c(-log10(counts_noNA[counts_noNA$Nac_log2fc_tp1 < 0,]$Nac_pval_tp1) * (-1), 
               -log10(counts_noNA[counts_noNA$Nac_log2fc_tp1 >= 0,]$Nac_pval_tp1)) 
gene_list1 <- data.frame(Genes=c(counts_noNA[counts_noNA$Nac_log2fc_tp1 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$Nac_log2fc_tp1 >= 0,]$MGI.symbol),
                         DDE = list1_DDE,
                         stringsAsFactors = FALSE)
# Genelist of TP2.
list2_DDE <- c(-log10(counts_noNA[counts_noNA$Nac_log2fc_tp2 < 0,]$Nac_pval_tp2) * (-1),
               -log10(counts_noNA[counts_noNA$Nac_log2fc_tp2 >= 0,]$Nac_pval_tp2)) 
gene_list2 <- data.frame(Genes=c(counts_noNA[counts_noNA$Nac_log2fc_tp2 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$Nac_log2fc_tp2 >= 0,]$MGI.symbol),
                         DDE = list2_DDE,
                         stringsAsFactors = FALSE)
# Genelist of TP3.
list3_DDE <- c(-log10(counts_noNA[counts_noNA$Nac_log2fc_tp3 < 0,]$Nac_pval_tp3) * (-1),
               -log10(counts_noNA[counts_noNA$Nac_log2fc_tp3 >= 0,]$Nac_pval_tp3)) 
gene_list3 <- data.frame(Genes=c(counts_noNA[counts_noNA$Nac_log2fc_tp3 < 0,]$MGI.symbol, 
                                 counts_noNA[counts_noNA$Nac_log2fc_tp3 >= 0,]$MGI.symbol),
                         DDE = list3_DDE,
                         stringsAsFactors = FALSE)

RRHO_Nac_1vs2 <- RRHO2_initialize(gene_list1, gene_list2, 
                                  labels = c("Tp1", "Tp2"), 
                                  log10.ind = TRUE)
RRHO_Nac_2vs3 <- RRHO2_initialize(gene_list2, gene_list3, 
                                  labels = c("Tp2", "Tp3"), 
                                  log10.ind = TRUE)
RRHO_Nac_1vs3 <- RRHO2_initialize(gene_list1, gene_list3,
                                  labels = c("Tp1", "Tp3"), 
                                  log10.ind = TRUE)

# save rrho objects
rrho_obj <- ls()[grepl("RRHO",ls())]
setwd(output.path)
save(list=rrho_obj, file="rrho_obj_multitp.Rdata")
