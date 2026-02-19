dir.create("~/R/library", showWarnings = FALSE, recursive = TRUE)
.libPaths(c("~/R/library", .libPaths()))

Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always") 

# Installer remotes si nécessaire
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Installer treeio depuis GitHub (corrige random_ref)
remotes::install_github("YuLab-SMU/treeio", lib="~/R/library")

# Installer ggtree depuis GitHub (corrige check_linewidth)
remotes::install_github("YuLab-SMU/ggtree", lib="~/R/library")

# Installer les autres dépendances de clusterProfiler
BiocManager::install(c("DOSE","enrichplot","clusterProfiler"), lib="~/R/library", ask=FALSE)


library(clusterProfiler)
library(tidyverse)

rm(list=ls())

setwd("/home/marinevernier/Documents/projets/")
#### PATH ####
annotated_counts.path <- "female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts_filtered.rds"
deglist.path <- "female_cpid_multiregion/data/2__differential_expression_analysis/deglist.Rdata"
output.path <- "female_cpid_multiregion/data/2__differential_expression_analysis/"
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

load(deglist.path)
counts <- readRDS(annotated_counts.path)
counts

regionList <- c("ACC", "Ins", "Hb", "Nac")
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
nrow(go_ACC_1)
for(reg in regionList){
  for(tp in tpList){
    go_obj <- paste("go", reg, tp, sep="_")
    print(paste0(reg, " tp", tp, ": ",nrow(get(go_obj))))
  }
}

# TESTER SANS BACKGROUND PERSO
view(go_Nac_1)
# test:
counts_temp <- counts %>% 
  dplyr::select(MGI.symbol, contains("Acc") & contains("2") & !contains("padj")) %>% 
  na.omit

deg_regtp <- get(paste0("deg_", reg, "_tp", tp))

go_test_Acc_2 <- enrichGO(gene = deg_regtp$label, 
                          OrgDb = organism, 
                          keyType = "SYMBOL", 
                          ont ="BP", 
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          pAdjustMethod = "BH",
                          universe = counts_temp$MGI.symbol)

go_test_Acc_2_wt_background <- enrichGO(gene = deg_regtp$label, 
                                        OrgDb = organism, 
                                        keyType = "SYMBOL", 
                                        ont ="BP", 
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.2,
                                        pAdjustMethod = "BH")

go_test_Acc_2_wt_pval_cutoff <- enrichGO(gene = deg_regtp$label, 
                                         OrgDb = organism, 
                                         keyType = "SYMBOL", 
                                         ont ="BP", 
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "BH",
                                         universe = counts_temp$MGI.symbol)

go_test_Acc_2_wt_background_and_pval_cutoff <- enrichGO(gene = deg_regtp$label, 
                                                        OrgDb = organism, 
                                                        keyType = "SYMBOL", 
                                                        ont ="BP", 
                                                        pvalueCutoff = 0.05,
                                                        pAdjustMethod = "BH")

nrow(go_test_Acc_2)
nrow(go_test_Acc_2_wt_background)
nrow(go_test_Acc_2_wt_pval_cutoff)
nrow(go_test_Acc_2_wt_background_and_pval_cutoff)
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
load("gse_obj.Rdata")

for(reg in regionList){
  for(tp in tpList){
    gse_obj <- paste("gse", reg, tp, sep="_")
    print(paste0(reg, " tp", tp, ": ",nrow(get(gse_obj))))
  }
}


#### intersect gsea ####
length(intersect(gse_ACC_1$Description, gse_ACC_2$Description)) 
length(intersect(gse_ACC_2$Description, gse_ACC_3$Description)) 
length(intersect(gse_ACC_1$Description, gse_ACC_3$Description))
length(intersect(gse_ACC_1$Description, intersect(gse_ACC_2$Description, gse_ACC_3$Description)))

length(intersect(gse_Hb_1$Description, gse_Hb_2$Description)) 
length(intersect(gse_Hb_2$Description, gse_Hb_3$Description)) 
length(intersect(gse_Hb_1$Description, gse_Hb_3$Description)) 
length(intersect(gse_Hb_1$Description, intersect(gse_Hb_2$Description, gse_Hb_3$Description))) 

length(intersect(gse_Ins_1$Description, gse_Ins_2$Description)) 
length(intersect(gse_Ins_2$Description, gse_Ins_3$Description)) 
length(intersect(gse_Ins_1$Description, gse_Ins_3$Description)) 

length(intersect(gse_Nac_1$Description, gse_Nac_2$Description)) 
length(intersect(gse_Nac_2$Description, gse_Nac_3$Description)) 
length(intersect(gse_Nac_1$Description, gse_Nac_3$Description)) 
length(intersect(gse_Nac_1$Description, intersect(gse_Nac_2$Description, gse_Nac_3$Description))) 


view(gse_Hb_1)
# big df with all results
acc1 <- data.frame(description=gse_ACC_1$Description, acc1_pval=gse_ACC_1$pvalue, acc1_ES=gse_ACC_1$enrichmentScore,
                   acc_ranking=sign(gse_ACC_1$enrichmentScore)*-log10(gse_ACC_1$pvalue))
acc2 <- data.frame(description=gse_ACC_2$Description, acc2_pval=gse_ACC_2$pvalue, acc2_ES=gse_ACC_2$enrichmentScore,
                   acc2_ranking=sign(gse_ACC_2$enrichmentScore)*-log10(gse_ACC_2$pvalue))
acc3 <- data.frame(description=gse_ACC_3$Description, acc3_pval=gse_ACC_3$pvalue, acc3_ES=gse_ACC_3$enrichmentScore,
                   acc3_ranking=sign(gse_ACC_3$enrichmentScore)*-log10(gse_ACC_3$pvalue))

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

nac1 <- data.frame(description=gse_Nac_1$Description, nac1_pval=gse_Nac_1$pvalue, nac1_ES=gse_Nac_1$enrichmentScore,
                   nac1_ranking=sign(gse_Nac_1$enrichmentScore)*-log10(gse_Nac_1$pvalue))
nac2 <- data.frame(description=gse_Nac_2$Description, nac2_pval=gse_Nac_2$pvalue, nac2_ES=gse_Nac_2$enrichmentScore,
                   nac2_ranking=sign(gse_Nac_2$enrichmentScore)*-log10(gse_Nac_2$pvalue))
nac3 <- data.frame(description=gse_Nac_3$Description, nac3_pval=gse_Nac_3$pvalue, nac3_ES=gse_Nac_3$enrichmentScore,
                   nac3_ranking=sign(gse_Nac_3$enrichmentScore)*-log10(gse_Nac_3$pvalue))

gse_allreg <- list(acc1, acc2, acc3,
                   hb1, hb2, hb3,
                   ins1, ins2, ins3,
                   nac1, nac2, nac3) %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="description"), .) %>% arrange(description)
write_rds(gse_allreg, "gse_allreg.rds")




#### PRINT GSEA RESULT ####
# load data
setwd("/home/marinevernier/Documents/projets/")
gse_allreg <- readRDS("female_cpid_multiregion/data/2__differential_expression_analysis/gse/gse_allreg.rds")

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

########################################################################################################################################

deg_summary <- data.frame()

for(reg in regionList){
  for(tp in tpList){
    
    # Sélection des colonnes correspondant à la région et tp
    df_temp <- counts %>%
      dplyr::select(
        MGI.symbol,
        contains(reg) & contains(paste0("tp", tp)) & 
          (contains("log2fc") | contains("pval") | contains("padj"))
      ) %>%
      na.omit()
    
    colnames(df_temp) <- c("gene", "log2fc", "pval", "padj")
    
    # Comptages
    n_pval  <- sum(df_temp$pval < 0.05)
    n_padj  <- sum(df_temp$padj < 0.05)
    
    deg_summary <- rbind(
      deg_summary,
      data.frame(
        region = reg,
        tp = tp,
        DEG_pvalue = n_pval,
        DEG_padj = n_padj
      )
    )
  }
}

deg_summary

