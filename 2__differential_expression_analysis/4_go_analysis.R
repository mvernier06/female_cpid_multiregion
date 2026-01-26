################# problème de biocmanager, faut R4.3 ####################





library(clusterProfiler)
library(tidyverse)

BiocManager::install("clusterProfiler")
BiocManager::install("ggtree")

rm(list=ls())

#### PATH ####
annotated_counts.path <- "~/CPID_multiregion/data/2__differential_expression_analysis/annotated_counts_filtered.rds"
deglist.path <- "~/CPID_multiregion/data/2__differential_expression_analysis/deglist.Rdata"
output.path <- "~/CPID_multiregion/data/2__differential_expression_analysis/"
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

load(deglist.path)
counts <- readRDS(annotated_counts.path)
counts

regionList <- c("BLA", "DRN", "Hb", "Ins", "NAc", "VTA")
tpList <- c("1", "2", "3")

# over-representation of DEGs vs background (hypergeometric)
for(reg in regionList){
  for(tp in tpList){
    print(paste0("Doing GO analysis of ", reg, " for tp", tp))
    counts_temp <- counts %>% 
      dplyr::select(MGI.symbol, contains(reg) & contains(tp) & !contains("padj")) %>% 
      na.omit
    
    deg_regtp <- get(paste0("deg_", reg, "_tp", tp))
    
    go_regtp <- enrichGO(gene = deg_regtp$label, 
                         OrgDb = organism, 
                         keyType = "SYMBOL", 
                         ont ="ALL", 
                         qvalueCutoff = 0.1,
                         pAdjustMethod = "BH",
                         universe = counts_temp$MGI.symbol)
    
    assign(paste("go", reg, tp, sep="_"), go_regtp, envir = .GlobalEnv)
  }
}


# save go objects #
rm(go_regtp)
go_obj <- ls()[grepl("go",ls())]
setwd(output.path)
save(list=go_obj, file="go_obj.Rdata")

load(paste0(output.path, "go_obj.Rdata"))
nrow(go_DRN_1)
for(reg in regionList){
  for(tp in tpList){
    go_obj <- paste("go", reg, tp, sep="_")
    print(paste0(reg, " tp", tp, ": ",nrow(get(go_obj))))
  }
}
# => Ins tp1+3 n'ont pas d'enrichissement GO
# TESTER SANS BACKGROUND PERSO
view(go_DRN_3)
# test:
counts_temp <- counts %>% 
  dplyr::select(MGI.symbol, contains("Ins") & contains("1") & !contains("padj")) %>% 
  na.omit

deg_regtp <- get(paste0("deg_", reg, "_tp", tp))

go_test_Ins_1 <- enrichGO(gene = deg_regtp$label, 
                          OrgDb = organism, 
                          keyType = "SYMBOL", 
                          ont ="BP", 
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          pAdjustMethod = "BH",
                          universe = counts_temp$MGI.symbol)

go_test_Ins_1_wt_background <- enrichGO(gene = deg_regtp$label, 
                                        OrgDb = organism, 
                                        keyType = "SYMBOL", 
                                        ont ="BP", 
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.2,
                                        pAdjustMethod = "BH")

go_test_Ins_1_wt_pval_cutoff <- enrichGO(gene = deg_regtp$label, 
                                         OrgDb = organism, 
                                         keyType = "SYMBOL", 
                                         ont ="BP", 
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "BH",
                                         universe = counts_temp$MGI.symbol)

go_test_Ins_1_wt_background_and_pval_cutoff <- enrichGO(gene = deg_regtp$label, 
                                                        OrgDb = organism, 
                                                        keyType = "SYMBOL", 
                                                        ont ="BP", 
                                                        pvalueCutoff = 0.05,
                                                        pAdjustMethod = "BH")

nrow(go_test_Ins_1)
nrow(go_test_Ins_1_wt_background)
nrow(go_test_Ins_1_wt_pval_cutoff)
nrow(go_test_Ins_1_wt_background_and_pval_cutoff)
# qval cutoff change rien, on perd le signal en utilisant notre background custom
# end test

#### GSEA ####
# All genes -> detect small consistent changes
for(reg in regionList){
  for(tp in tpList){
    print(c(reg, tp))
    counts_temp <- counts %>% 
      dplyr::select(MGI.symbol, contains(reg) & contains(paste0("tp", tp)) & contains("log2fc")) %>% 
      na.omit
    colnames(counts_temp) <- c("MGI.symbol", "log2fc")
    
    log2fc <- counts_temp$log2fc
    names(log2fc) <- counts_temp$MGI.symbol
    genelist <- sort(log2fc, decreasing = TRUE)
    
    gse <- gseGO(geneList=genelist, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose =  TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "BH",
                 seed = 123)
    
    assign(paste("gse", reg, tp, sep="_"), gse)
  }
}

# save gse objects
gse_obj <- ls()[grepl("gse",ls())]
setwd(output.path)
dir.create("./gse")
setwd("./gse")
save(list=gse_obj, file="gse_obj.Rdata")
load("~/CPID_multiregion/data/2__differential_expression_analysis/gse/gse_obj.Rdata")

for(reg in regionList){
  for(tp in tpList){
    gse_obj <- paste("gse", reg, tp, sep="_")
    print(paste0(reg, " tp", tp, ": ",nrow(get(gse_obj))))
  }
}


#### intersect gsea ####
length(intersect(gse_BLA_1$Description, gse_BLA_2$Description)) # 75   après correction #116
length(intersect(gse_BLA_2$Description, gse_BLA_3$Description)) # 2  # 3
length(intersect(gse_BLA_1$Description, gse_BLA_3$Description)) # 2  # 4
length(intersect(gse_BLA_1$Description, intersect(gse_BLA_2$Description, gse_BLA_3$Description))) # 0

length(intersect(gse_DRN_1$Description, gse_DRN_2$Description)) # 64   # 114
length(intersect(gse_DRN_2$Description, gse_DRN_3$Description)) # 19   # 42
length(intersect(gse_DRN_1$Description, gse_DRN_3$Description)) # 16   # 50
length(intersect(gse_DRN_1$Description, intersect(gse_DRN_2$Description, gse_DRN_3$Description))) # 10 # 28

length(intersect(gse_Hb_1$Description, gse_Hb_2$Description)) # 60 # 96
length(intersect(gse_Hb_2$Description, gse_Hb_3$Description)) # 117 # 231
length(intersect(gse_Hb_1$Description, gse_Hb_3$Description)) # 39 # 66
length(intersect(gse_Hb_1$Description, intersect(gse_Hb_2$Description, gse_Hb_3$Description))) # 33 # 56

length(intersect(gse_Ins_1$Description, gse_Ins_2$Description)) # 0 
length(intersect(gse_Ins_2$Description, gse_Ins_3$Description)) # 0
length(intersect(gse_Ins_1$Description, gse_Ins_3$Description)) # 0 # 1

length(intersect(gse_NAc_1$Description, gse_NAc_2$Description)) # 95 # 123
length(intersect(gse_NAc_2$Description, gse_NAc_3$Description)) # 37 # 58
length(intersect(gse_NAc_1$Description, gse_NAc_3$Description)) # 7  # 12 
length(intersect(gse_NAc_1$Description, intersect(gse_NAc_2$Description, gse_NAc_3$Description))) # 6 # 10

length(intersect(gse_VTA_1$Description, gse_VTA_2$Description)) # 51 # 97
length(intersect(gse_VTA_2$Description, gse_VTA_3$Description)) # 59 # 101
length(intersect(gse_VTA_1$Description, gse_VTA_3$Description)) # 283 # 394
length(intersect(gse_VTA_1$Description, intersect(gse_VTA_2$Description, gse_VTA_3$Description))) # 37 # 73

view(gse_BLA_1)
# big df with all results
bla1 <- data.frame(description=gse_BLA_1$Description, bla1_pval=gse_BLA_1$pvalue, bla1_ES=gse_BLA_1$enrichmentScore,
                   bla1_ranking=sign(gse_BLA_1$enrichmentScore)*-log10(gse_BLA_1$pvalue))
bla2 <- data.frame(description=gse_BLA_2$Description, bla2_pval=gse_BLA_2$pvalue, bla2_ES=gse_BLA_2$enrichmentScore,
                   bla2_ranking=sign(gse_BLA_2$enrichmentScore)*-log10(gse_BLA_2$pvalue))
bla3 <- data.frame(description=gse_BLA_3$Description, bla3_pval=gse_BLA_3$pvalue, bla3_ES=gse_BLA_3$enrichmentScore,
                   bla3_ranking=sign(gse_BLA_3$enrichmentScore)*-log10(gse_BLA_3$pvalue))

drn1 <- data.frame(description=gse_DRN_1$Description, drn1_pval=gse_DRN_1$pvalue, drn1_ES=gse_DRN_1$enrichmentScore,
                   drn1_ranking=sign(gse_DRN_1$enrichmentScore)*-log10(gse_DRN_1$pvalue))
drn2 <- data.frame(description=gse_DRN_2$Description, drn2_pval=gse_DRN_2$pvalue, drn2_ES=gse_DRN_2$enrichmentScore,
                   drn2_ranking=sign(gse_DRN_2$enrichmentScore)*-log10(gse_DRN_2$pvalue))
drn3 <- data.frame(description=gse_DRN_3$Description, drn3_pval=gse_DRN_3$pvalue, drn3_ES=gse_DRN_3$enrichmentScore,
                   drn3_ranking=sign(gse_DRN_3$enrichmentScore)*-log10(gse_DRN_3$pvalue))

hb1 <- data.frame(description=gse_Hb_1$Description, hb1_pval=gse_Hb_1$pvalue, hb1_ES=gse_Hb_1$enrichmentScore,
                  hb1_ranking=sign(gse_Hb_1$enrichmentScore)*-log10(gse_Hb_1$pvalue))
hb2 <- data.frame(description=gse_Hb_2$Description, hb2_pval=gse_Hb_2$pvalue, hb2_ES=gse_Hb_2$enrichmentScore,
                  hb2_ranking=sign(gse_Hb_2$enrichmentScore)*-log10(gse_Hb_2$pvalue))
hb3 <- data.frame(description=gse_Hb_3$Description, hb3_pval=gse_Hb_3$pvalue, hb3_ES=gse_Hb_3$enrichmentScore,
                  hb3_ranking=sign(gse_Hb_3$enrichmentScore)*-log10(gse_Hb_3$pvalue))

ins1 <- data.frame(description=gse_Ins_1$Description, ins1_pval=gse_Ins_1$pvalue, ins1_ES=gse_Ins_1$enrichmentScore,
                   ins1_ranking=sign(gse_Ins_1$enrichmentScore)*-log10(gse_Ins_1$pvalue))
ins2 <- data.frame(description=gse_Ins_2$Description, ins2_pval=gse_Ins_2$pvalue, ins2_ES=gse_Ins_2$enrichmentScore,
                   ins2_ranking=sign(gse_Ins_2$enrichmentScore)*-log10(gse_Ins_2$pvalue))
ins3 <- data.frame(description=gse_Ins_3$Description, ins3_pval=gse_Ins_3$pvalue, ins3_ES=gse_Ins_3$enrichmentScore,
                   ins3_ranking=sign(gse_Ins_3$enrichmentScore)*-log10(gse_Ins_3$pvalue))

nac1 <- data.frame(description=gse_NAc_1$Description, nac1_pval=gse_NAc_1$pvalue, nac1_ES=gse_NAc_1$enrichmentScore,
                   nac1_ranking=sign(gse_NAc_1$enrichmentScore)*-log10(gse_NAc_1$pvalue))
nac2 <- data.frame(description=gse_NAc_2$Description, nac2_pval=gse_NAc_2$pvalue, nac2_ES=gse_NAc_2$enrichmentScore,
                   nac2_ranking=sign(gse_NAc_2$enrichmentScore)*-log10(gse_NAc_2$pvalue))
nac3 <- data.frame(description=gse_NAc_3$Description, nac3_pval=gse_NAc_3$pvalue, nac3_ES=gse_NAc_3$enrichmentScore,
                   nac3_ranking=sign(gse_NAc_3$enrichmentScore)*-log10(gse_NAc_3$pvalue))

vta1 <- data.frame(description=gse_VTA_1$Description, vta1_pval=gse_VTA_1$pvalue, vta1_ES=gse_VTA_1$enrichmentScore,
                   vta1_ranking=sign(gse_VTA_1$enrichmentScore)*-log10(gse_VTA_1$pvalue))
vta2 <- data.frame(description=gse_VTA_2$Description, vta2_pval=gse_VTA_2$pvalue, vta2_ES=gse_VTA_2$enrichmentScore,
                   vta2_ranking=sign(gse_VTA_2$enrichmentScore)*-log10(gse_VTA_2$pvalue))
vta3 <- data.frame(description=gse_VTA_3$Description, vta3_pval=gse_VTA_3$pvalue, vta3_ES=gse_VTA_3$enrichmentScore,
                   vta3_ranking=sign(gse_VTA_3$enrichmentScore)*-log10(gse_VTA_3$pvalue))

gse_allreg <- list(bla1, bla2, bla3,
                   drn1, drn2, drn3,
                   hb1, hb2, hb3,
                   ins1, ins2, ins3,
                   nac1, nac2, nac3,
                   vta1, vta2, vta3) %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="description"), .) %>% arrange(description)
write_rds(gse_allreg, "gse_allreg.rds")




#### PRINT GSEA RESULT ####
# load data
gse_allreg <- readRDS("~/CPID_multiregion/data/2__differential_expression_analysis/gse/gse_allreg.rds")

gse_heatmap <- function(query){
  # format data
  test <- gse_allreg %>% 
    filter(grepl(query, description)) %>%
    dplyr::select(description, contains("ranking"))
  test[,-1] <- lapply(test[,-1], function(x) round(x, digits=2))
  test2 <- test %>% column_to_rownames("description")
  
  # max absolute value to scale the color scale
  maxval <- test2 %>% 
    apply(2, function(x) replace_na(x, 0)) %>% 
    as.matrix %>% abs %>% max
  
  # make a table for all pathway found 
  library(pheatmap)
  p <- pheatmap(test2,
                display_numbers = TRUE,
                fontsize = 8,
                breaks = seq(-maxval, maxval, length.out = 100), # center the color scale on zero
                cluster_cols = FALSE,
                cluster_rows = FALSE,
                labels_col = c(paste("bla", c(1,2,3), sep = "_"),
                               paste("drn", c(1,2,3), sep = "_"),
                               paste("hb", c(1,2,3), sep = "_"),
                               paste("ins", c(1,2,3), sep = "_"),
                               paste("nac", c(1,2,3), sep = "_"),
                               paste("vta", c(1,2,3), sep = "_")),
                angle_col = 315,
                annotation_colors = "test")
  print(p)
}

gse_heatmap("myelin")
gse_heatmap("calcium")