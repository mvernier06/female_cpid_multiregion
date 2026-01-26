library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggalluvial)

rm(list=ls())

# PATHS ####
annotated_counts.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts_filtered.rds"
output.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/"
plot.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/deg_volcano"
dir.create(plot.path)

annotated_counts <- read_rds(annotated_counts.path)
setwd(plot.path)

regionlist <- c("ACC", "Hb", "Ins", "Nac")
tplist <- c("tp1", "tp2", "tp3")

# volcano plots for each tp of every reg
for(reg in regionlist){
  for(tp in tplist){
    temp_counts <- annotated_counts %>% select(c(MGI.symbol, contains(reg) & contains(tp)))
    
    # rename temp_counts colnames so that they are universal for all samples
    temp_counts <- temp_counts %>% 
      rename_at(vars(contains("log2fc")), ~ "log2fc") %>%
      rename_at(vars(contains("pval")), ~ "pval") %>%
      rename_at(vars(contains("padj")), ~ "padj")
    
    temp_counts$diffexpressed <- "ns"
    
    temp_counts$diffexpressed[temp_counts$log2fc > 0  & temp_counts$pval < 0.05] <- "UP"
    temp_counts$diffexpressed[temp_counts$log2fc < 0  & temp_counts$pval < 0.05] <- "DOWN"
    
    temp_counts$label <- NA
    temp_counts$label[temp_counts$diffexpressed != "ns"] <- temp_counts$MGI.symbol[temp_counts$diffexpressed != "ns"]
    
    assign(paste0("temp_counts_", reg, "_", tp), temp_counts, envir = .GlobalEnv)
    
    # volcano plot tp1
    p <- ggplot(data=temp_counts, aes(x=log2fc, y=-log10(pval), col=diffexpressed)) +
      geom_point() + 
      theme_minimal() +
      scale_color_manual(values=c("blue", "grey", "red")) +
      geom_vline(xintercept=c(-1, 1), col="grey", lty="dashed") +
      geom_hline(yintercept=-log10(0.05), col="grey", lty="dashed") +
      guides(color = guide_legend(override.aes = list(linetype = 1, size=2)), 
             label = "test") +
      labs(title=paste0("DEGs of ", reg, " at ", tp),
           x="log2 (Fold Change)", y="-log10 (p-value)",
           colour="DEG")
    plot.filename <- paste0("deg_volcano_", reg,"_" ,tp, ".png")
    ggsave(plot=p, filename=plot.filename, bg="white", width=1900, height=1200, units="px", scale=2)
    
    # return DEGs
    deg <- temp_counts %>% 
      dplyr::select(label, pval, log2fc, diffexpressed) %>%
      na.omit() %>%
      arrange(pval)
    assign(paste0("deg_", reg, "_", tp), deg, envir = .GlobalEnv)
  }
}

# save degs
setwd(output.path)
deglist <- ls()[grepl("deg_",ls())]
save(list=deglist, file="deglist.Rdata")
# save temp_counts
tempcountlist <- ls()[grepl("temp_counts_",ls())]
save(list=tempcountlist, file="tempcountlist.Rdata")
