.libPaths(c("~/R/library", .libPaths()))
Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always") 
library(clusterProfiler)
library(rrvgo)
library(org.Mm.eg.db)
library(dplyr)
library(openxlsx)
library(pheatmap)
library(tidyverse)
library(ggplot2)

rm(list=ls())

project_path <- "/home/marinevernier/Documents/projets/"
setwd(project_path)

intersections.path <- "female_cpid_multiregion/data/0_comparison_male_female/overlaps_rrho/intersection_terms_overlap_rrho_male_vs_female_summary.xlsx"
intersections_full.path <- "female_cpid_multiregion/data/0_comparison_male_female/overlaps_rrho/intersection_terms_overlap_rrho_male_vs_female.xlsx"
output.path <- "female_cpid_multiregion/graphs_results/0_comparison_male_female/overlaps_rrho/"
dir.create(output.path)

intersections <- read.xlsx(intersections.path)
intersection_full <- read.xlsx(intersections_full.path)


comp <- "rrhoHbvsNac1_down_down"

for (comp in unique(intersection_full$comparison)) {
  df <- intersection_full %>% filter(comparison == comp) %>%
    mutate(    
      log_male   = -log10(padj_male),
      log_female = -log10(padj_female)
    )
  
 p <- ggplot(df, aes(x = log_male, y = log_female, col = Description)) + 
   geom_point(alpha = 0.9) +
   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
   scale_y_continuous(limits = c(0, 30)) + 
   scale_x_continuous(limits = c(0, 30)) +
   labs(
     x = "-log10(padj) Male",
     y = "-log10(padj) Female",
     title = paste0("Enrichment overlap ", comp, " male vs female")
   ) +
   theme(lengend.text = element_text(size = 0.5) )
 p
 
 ggsave(p, filename=paste0(output.path, "/comparison_pvalue_male_female_overlaps_", comp,".png"))
}

