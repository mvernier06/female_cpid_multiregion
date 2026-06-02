library(tidyverse)
library(RRHO2)
library(openxlsx)

rm(list=ls())
setwd("/home/marinevernier/Documents/projets/")

# PATHS ####
counts_female.path <- "female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts_filtered.rds"
counts_male.path <- "cpid_multiregion/data/2__differential_expression_analysis/annotated_counts_filtered.rds"
output.path <- "female_cpid_multiregion/data/0_comparison_male_female/RRHO"
dir.create(output.path, showWarnings = FALSE, recursive = TRUE)
plot.path <- "female_cpid_multiregion/graphs_results/0_comparison_male_female/RRHO"
dir.create(plot.path, showWarnings = FALSE, recursive = TRUE)

# Import counts
counts_female <- readRDS(counts_female.path) %>% select("MGI.symbol", contains(c("log2fc","pval")) )
counts_male <- readRDS(counts_male.path) %>% select("MGI.symbol", contains(c("log2fc","pval")) )
colnames(counts_male) <- gsub("NAc", "Nac", colnames(counts_male))


common_regions <- c("Ins", "Nac", "Hb")
timepoints <- c("tp1", "tp2", "tp3")

for (reg in common_regions) {
  for (tp in timepoints) {
    counts_female_tmp <- counts_female %>% select("MGI.symbol", contains(reg)) %>% select("MGI.symbol", contains(tp))
    counts_male_tmp <- counts_male %>% select("MGI.symbol", contains(reg)) %>% select("MGI.symbol", contains(tp))
    
    counts_both <- inner_join(counts_female_tmp, counts_male_tmp, by = "MGI.symbol", suffix = c("_female", "_male")) 
    counts_noNA <- counts_both %>% 
      filter(!is.na(.data[[paste0(reg,"_log2fc_",tp,"_female")]]) &
               !is.na(.data[[paste0(reg,"_pval_",tp, "_female")]]) &
               !is.na(.data[[paste0(reg,"_log2fc_",tp,"_male")]]) &
               !is.na(.data[[paste0(reg,"_pval_",tp, "_male")]]))
    
    # Genelist of female.
    log2fc_col <- paste0(reg, "_log2fc_", tp, "_female")
    pval_col   <- paste0(reg, "_pval_", tp, "_female")

    list1_DDE <- c(
      -log10(counts_noNA[[pval_col]][counts_noNA[[log2fc_col]] < 0]) * (-1),
      -log10(counts_noNA[[pval_col]][counts_noNA[[log2fc_col]] > 0])
    )
    
    gene_list1 <- data.frame(Genes=c(counts_noNA$MGI.symbol[counts_noNA[[log2fc_col]] < 0], 
                       counts_noNA$MGI.symbol[counts_noNA[[log2fc_col]] > 0]),
               DDE = list1_DDE,
               stringsAsFactors = FALSE)
    
  
    
    # Genelist of male.
    log2fc_col <- paste0(reg, "_log2fc_", tp, "_male")
    pval_col   <- paste0(reg, "_pval_", tp, "_male")
    
    list2_DDE <- c(
      -log10(counts_noNA[[pval_col]][counts_noNA[[log2fc_col]] < 0]) * (-1),
      -log10(counts_noNA[[pval_col]][counts_noNA[[log2fc_col]] > 0])
    )
    
    gene_list2 <- data.frame(Genes=c(counts_noNA$MGI.symbol[counts_noNA[[log2fc_col]] < 0], 
                      counts_noNA$MGI.symbol[counts_noNA[[log2fc_col]] > 0]),
              DDE = list1_DDE,
              stringsAsFactors = FALSE)
    
    RRHO_femalevsmale_l2fc_pval <- RRHO2_initialize(gene_list1, gene_list2, 
                                                    labels = c("female", "male"), 
                                                    log10.ind = TRUE)
    filename <- file.path(output.path, paste0(reg,"_rrho_male_female_",tp,".Rdata"))
    save(RRHO_femalevsmale_l2fc_pval, file = filename)
    filename <- file.path(plot.path,paste0("/",reg,"/",reg, "_male_vs_female_", tp, ".png"))
    png(file=filename, width = 750, height = 500, units = "px", res = 100)
    RRHO2_heatmap(RRHO_femalevsmale_l2fc_pval, main = paste0(reg, "male vs female ", tp))
    dev.off()
    
    }
}
