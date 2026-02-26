library(tidyverse)
library(RRHO2)

rm(list=ls())
setwd("/home/marinevernier/Documents/projets/female_cpid_multiregion/")

# PATHS ####
annotated_counts.path <- "data/2__differential_expression_analysis/annotated_counts_filtered.rds"
output.path <- "data/3__RRHO2/"
dir.create(output.path, showWarnings = FALSE, recursive = TRUE)

# Import counts
annotated_counts <- readRDS(annotated_counts.path)

# Filter out NA values and padj
raw_noNA_tp1 <- annotated_counts %>% 
  select(c(MGI.symbol, contains("tp1") & !contains("padj"))) %>%
  na.omit()

raw_noNA_tp2 <- annotated_counts %>% 
  select(c(MGI.symbol, contains("tp2") & !contains("padj"))) %>%
  na.omit()

raw_noNA_tp3 <- annotated_counts %>% 
  select(c(MGI.symbol, contains("tp3") & !contains("padj"))) %>%
  na.omit()

#### TP1 ####
print("Computing RRHO2 heatmaps of the TP1")

# Genelist of ACC
listDDE <- c(-log10(raw_noNA_tp1[raw_noNA_tp1$ACC_log2fc_tp1 < 0,]$ACC_pval_tp1) * (-1), 
             -log10(raw_noNA_tp1[raw_noNA_tp1$ACC_log2fc_tp1 > 0,]$ACC_pval_tp1))
geneListACC1 <- data.frame(Genes=c(raw_noNA_tp1[raw_noNA_tp1$ACC_log2fc_tp1 < 0,]$MGI.symbol, 
                                   raw_noNA_tp1[raw_noNA_tp1$ACC_log2fc_tp1 > 0,]$MGI.symbol),
                           DDE = listDDE,
                           stringsAsFactors = FALSE)

# Genelist of Hb.
listDDE <- c(-log10(raw_noNA_tp1[raw_noNA_tp1$Hb_log2fc_tp1 < 0,]$Hb_pval_tp1) * (-1), 
             -log10(raw_noNA_tp1[raw_noNA_tp1$Hb_log2fc_tp1 > 0,]$Hb_pval_tp1))
geneListHb1 <- data.frame(Genes=c(raw_noNA_tp1[raw_noNA_tp1$Hb_log2fc_tp1 <= 0,]$MGI.symbol, 
                                  raw_noNA_tp1[raw_noNA_tp1$Hb_log2fc_tp1 > 0,]$MGI.symbol),
                          DDE = listDDE,
                          stringsAsFactors = FALSE)

# Genelist of Ins.
listDDE <- c(-log10(raw_noNA_tp1[raw_noNA_tp1$Ins_log2fc_tp1 < 0,]$Ins_pval_tp1) * (-1), 
             -log10(raw_noNA_tp1[raw_noNA_tp1$Ins_log2fc_tp1 > 0,]$Ins_pval_tp1))
geneListIns1 <- data.frame(Genes=c(raw_noNA_tp1[raw_noNA_tp1$Ins_log2fc_tp1 <= 0,]$MGI.symbol, 
                                   raw_noNA_tp1[raw_noNA_tp1$Ins_log2fc_tp1 > 0,]$MGI.symbol),
                           DDE = listDDE,
                           stringsAsFactors = FALSE)

# Genelist of NAc.
listDDE <- c(-log10(raw_noNA_tp1[raw_noNA_tp1$Nac_log2fc_tp1 < 0,]$Nac_pval_tp1) * (-1), 
             -log10(raw_noNA_tp1[raw_noNA_tp1$Nac_log2fc_tp1 > 0,]$Nac_pval_tp1))
geneListNac1 <- data.frame(Genes=c(raw_noNA_tp1[raw_noNA_tp1$Nac_log2fc_tp1 <= 0,]$MGI.symbol, 
                                   raw_noNA_tp1[raw_noNA_tp1$Nac_log2fc_tp1 > 0,]$MGI.symbol),
                           DDE = listDDE,
                           stringsAsFactors = FALSE)


# ACC
rrhoACCvsHb1 <- RRHO2_initialize(geneListACC1, geneListHb1,
                                 labels = c("ACC1", "Hb1"),
                                 log10.ind = TRUE)
rrhoACCvsIns1 <- RRHO2_initialize(geneListACC1, geneListIns1,
                                  labels = c("ACC1", "Ins1"),
                                  log10.ind = TRUE)
rrhoACCvsNac1 <- RRHO2_initialize(geneListACC1, geneListNac1,
                                  labels = c("ACC1", "Nac1"),
                                  log10.ind = TRUE)

# Hb
rrhoHbvsIns1 <- RRHO2_initialize(geneListHb1, geneListIns1,
                                 labels = c("Hb1", "Ins1"),
                                 log10.ind = TRUE)
rrhoHbvsNac1 <- RRHO2_initialize(geneListHb1, geneListNac1,
                                 labels = c("Hb1", "Nac1"),
                                 log10.ind = TRUE)

# Ins
rrhoInsvsNac1 <- RRHO2_initialize(geneListIns1, geneListNac1,
                                  labels = c("Ins1", "Nac1"),
                                  log10.ind = TRUE)



#### TP2 ####
print("Computing RRHO2 heatmaps of the TP2")

# Genelist of ACC
listDDE <- c(-log10(raw_noNA_tp2[raw_noNA_tp2$ACC_log2fc_tp2 < 0,]$ACC_pval_tp2) * (-1), 
             -log10(raw_noNA_tp2[raw_noNA_tp2$ACC_log2fc_tp2 > 0,]$ACC_pval_tp2))
geneListACC2 <- data.frame(Genes=c(raw_noNA_tp2[raw_noNA_tp2$ACC_log2fc_tp2 <= 0,]$MGI.symbol, 
                                   raw_noNA_tp2[raw_noNA_tp2$ACC_log2fc_tp2 > 0,]$MGI.symbol),
                           DDE = listDDE,
                           stringsAsFactors = FALSE)

# Genelist of Hb.
listDDE <- c(-log10(raw_noNA_tp2[raw_noNA_tp2$Hb_log2fc_tp2 < 0,]$Hb_pval_tp2) * (-1), 
             -log10(raw_noNA_tp2[raw_noNA_tp2$Hb_log2fc_tp2 > 0,]$Hb_pval_tp2))
geneListHb2 <- data.frame(Genes=c(raw_noNA_tp2[raw_noNA_tp2$Hb_log2fc_tp2 <= 0,]$MGI.symbol, 
                                  raw_noNA_tp2[raw_noNA_tp2$Hb_log2fc_tp2 > 0,]$MGI.symbol),
                          DDE = listDDE,
                          stringsAsFactors = FALSE)

# Genelist of Ins.
listDDE <- c(-log10(raw_noNA_tp2[raw_noNA_tp2$Ins_log2fc_tp2 < 0,]$Ins_pval_tp2) * (-1), 
             -log10(raw_noNA_tp2[raw_noNA_tp2$Ins_log2fc_tp2 > 0,]$Ins_pval_tp2))
geneListIns2 <- data.frame(Genes=c(raw_noNA_tp2[raw_noNA_tp2$Ins_log2fc_tp2 <= 0,]$MGI.symbol, 
                                   raw_noNA_tp2[raw_noNA_tp2$Ins_log2fc_tp2 > 0,]$MGI.symbol),
                           DDE = listDDE,
                           stringsAsFactors = FALSE)

# Genelist of Nac.
listDDE <- c(-log10(raw_noNA_tp2[raw_noNA_tp2$Nac_log2fc_tp2 < 0,]$Nac_pval_tp2) * (-1), 
             -log10(raw_noNA_tp2[raw_noNA_tp2$Nac_log2fc_tp2 >= 0,]$Nac_pval_tp2))
geneListNac2 <- data.frame(Genes=c(raw_noNA_tp2[raw_noNA_tp2$Nac_log2fc_tp2 < 0,]$MGI.symbol, 
                                   raw_noNA_tp2[raw_noNA_tp2$Nac_log2fc_tp2 >= 0,]$MGI.symbol),
                           DDE = listDDE,
                           stringsAsFactors = FALSE)


# ACC
rrhoACCvsHb2 <- RRHO2_initialize(geneListACC2, geneListHb2,
                                 labels = c("ACC2", "Hb2"),
                                 log10.ind = TRUE)
rrhoACCvsIns2 <- RRHO2_initialize(geneListACC2, geneListIns2,
                                  labels = c("ACC2", "Ins2"),
                                  log10.ind = TRUE)
rrhoACCvsNac2 <- RRHO2_initialize(geneListACC2, geneListNac2,
                                  labels = c("ACC2", "Nac2"),
                                  log10.ind = TRUE)
# Hb
rrhoHbvsIns2 <- RRHO2_initialize(geneListHb2, geneListIns2,
                                 labels = c("Hb2", "Ins2"),
                                 log10.ind = TRUE)
rrhoHbvsNac2 <- RRHO2_initialize(geneListHb2, geneListNac2,
                                 labels = c("Hb2", "Nac2"),
                                 log10.ind = TRUE)
# Ins
rrhoInsvsNac2 <- RRHO2_initialize(geneListIns2, geneListNac2,
                                  labels = c("Ins2", "Nac2"),
                                  log10.ind = TRUE)


#### TP3 ####
print("Computing RRHO2 heatmaps of the TP3")

# Genelist of ACC.
listDDE <- c(-log10(raw_noNA_tp3[raw_noNA_tp3$ACC_log2fc_tp3 < 0,]$ACC_pval_tp3) * (-1), 
             -log10(raw_noNA_tp3[raw_noNA_tp3$ACC_log2fc_tp3 > 0,]$ACC_pval_tp3))
geneListACC3 <- data.frame(Genes=c(raw_noNA_tp3[raw_noNA_tp3$ACC_log2fc_tp3 <= 0,]$MGI.symbol, 
                                   raw_noNA_tp3[raw_noNA_tp3$ACC_log2fc_tp3 > 0,]$MGI.symbol),
                           DDE = listDDE,
                           stringsAsFactors = FALSE)


# Genelist of Hb.
listDDE <- c(-log10(raw_noNA_tp3[raw_noNA_tp3$Hb_log2fc_tp3 < 0,]$Hb_pval_tp3) * (-1), 
             -log10(raw_noNA_tp3[raw_noNA_tp3$Hb_log2fc_tp3 > 0,]$Hb_pval_tp3))
geneListHb3 <- data.frame(Genes=c(raw_noNA_tp3[raw_noNA_tp3$Hb_log2fc_tp3 <= 0,]$MGI.symbol, 
                                  raw_noNA_tp3[raw_noNA_tp3$Hb_log2fc_tp3 > 0,]$MGI.symbol),
                          DDE = listDDE,
                          stringsAsFactors = FALSE)

# Genelist of Ins.
listDDE <- c(-log10(raw_noNA_tp3[raw_noNA_tp3$Ins_log2fc_tp3 < 0,]$Ins_pval_tp3) * (-1), 
             -log10(raw_noNA_tp3[raw_noNA_tp3$Ins_log2fc_tp3 > 0,]$Ins_pval_tp3))
geneListIns3 <- data.frame(Genes=c(raw_noNA_tp3[raw_noNA_tp3$Ins_log2fc_tp3 <= 0,]$MGI.symbol, 
                                   raw_noNA_tp3[raw_noNA_tp3$Ins_log2fc_tp3 > 0,]$MGI.symbol),
                           DDE = listDDE,
                           stringsAsFactors = FALSE)

# Genelist of Nac.
listDDE <- c(-log10(raw_noNA_tp3[raw_noNA_tp3$Nac_log2fc_tp3 < 0,]$Nac_pval_tp3) * (-1), 
             -log10(raw_noNA_tp3[raw_noNA_tp3$Nac_log2fc_tp3 > 0,]$Nac_pval_tp3))
geneListNac3 <- data.frame(Genes=c(raw_noNA_tp3[raw_noNA_tp3$Nac_log2fc_tp3 <= 0,]$MGI.symbol, 
                                   raw_noNA_tp3[raw_noNA_tp3$Nac_log2fc_tp3 > 0,]$MGI.symbol),
                           DDE = listDDE,
                           stringsAsFactors = FALSE)



# ACC
rrhoACCvsHb3 <- RRHO2_initialize(geneListACC3, geneListHb3,
                                 labels = c("ACC3", "Hb3"),
                                 log10.ind = TRUE)
rrhoACCvsIns3 <- RRHO2_initialize(geneListACC3, geneListIns3,
                                  labels = c("ACC3", "Ins3"),
                                  log10.ind = TRUE)
rrhoACCvsNac3 <- RRHO2_initialize(geneListACC3, geneListNac3,
                                  labels = c("ACC3", "Nac3"),
                                  log10.ind = TRUE)

# Hb
rrhoHbvsIns3 <- RRHO2_initialize(geneListHb3, geneListIns3,
                                 labels = c("Hb3", "Ins3"),
                                 log10.ind = TRUE)
rrhoHbvsNac3 <- RRHO2_initialize(geneListHb3, geneListNac3,
                                 labels = c("Hb3", "Nac3"),
                                 log10.ind = TRUE)

# Ins
rrhoInsvsNac3 <- RRHO2_initialize(geneListIns3, geneListNac3,
                                  labels = c("Ins3", "Nac3"),
                                  log10.ind = TRUE)


# save rrho objects
rrho_obj <- ls()[grepl("rrho",ls())]
setwd(output.path)
save(list=rrho_obj, file="rrho_obj_multireg.Rdata")
