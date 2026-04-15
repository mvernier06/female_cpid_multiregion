.libPaths(c("~/R/library", .libPaths()))

Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always") 

library(tidyverse)
library(tibble)
library(ggplot2)
library(clusterProfiler)

rm(list=ls())

project.path <- "/home/marinevernier/Documents/projets/"
setwd(project.path)

## Path 
fgsea_female <- "female_cpid_multiregion/data/2__differential_expression_analysis/fgsea_obj.Rdata"
fgsea_male <- "cpid_supplementary/2__differential_expression_analysis/data/fgsea_obj.Rdata"
output_data_path <- "female_cpid_multiregion/data/0_comparison_male_female/"
output_path <- "female_cpid_multiregion/graphs_results/0_comparison_male_female/"

## Load 
female_env <- new.env()
load(fgsea_female, envir = female_env)
female_objs <- ls(female_env, pattern = "^fgsea_")
for (obj in female_objs) {
  assign(paste0(obj, "_female"), female_env[[obj]])
}

male_env <- new.env()
load(fgsea_male, envir = male_env)
male_objs <- ls(male_env, pattern = "^fgsea_")
for (obj in male_objs) {
  assign(paste0(obj, "_male"), male_env[[obj]])
}

# Homogénéiser les noms 
fgsea_NAc_1_female <- fgsea_Nac_1_female
fgsea_NAc_2_female <- fgsea_Nac_2_female
fgsea_NAc_3_female <- fgsea_Nac_3_female


timepoints <- c("1", "2", "3")
regions <- c("Ins", "Hb", "NAc")
intersections_list <- list()

fgsea <- fgsea_Hb_1_female

get_terms_description <- function(fgsea, pval_threshold = 0.05, padj_threshold = 0.05, qval_threshold = 0.05,  NES_threshold = 0){
  fgsea@result %>%
    filter(pvalue < pval_threshold & p.adjust < padj_threshold & qvalue < qval_threshold & abs(NES) > NES_threshold) %>%
    pull(Description)
}

get_terms_id <- function(fgsea, pval_threshold = 0.05, padj_threshold = 0.05, qval_threshold = 0.05,  NES_threshold = 0){
  fgsea@result %>%
    filter(pvalue < pval_threshold & p.adjust < padj_threshold & qvalue < qval_threshold & abs(NES) > NES_threshold) %>%
    pull(ID)
}

for(tp in timepoints){
  for(reg in regions){
    
    male_name <- paste0("fgsea_", reg, "_", tp, "_male")
    female_name <- paste0("fgsea_", reg, "_", tp, "_female")
    
    if(exists(male_name) & exists(female_name)){
      
      male_df <- get(male_name)
      female_df <- get(female_name)
      
      male_terms <- get_terms_description(male_df)
      female_terms <- get_terms_description(female_df)
      
      inter <- intersect(male_terms, female_terms)
      
      intersections_list[[paste(reg, tp, sep = "_")]] <- inter
      
    }
  }
}
intersection_sizes <- sapply(intersections_list, length)
intersection_sizes
save(intersections_list, file=paste0(output_data_path, "intersection_terms/intersection_list_fgsea_terms.Rdata"))


intersection_df <- tibble(
  comparison = names(intersections_list),
  n_terms = intersection_sizes
) %>%
  tidyr::separate(comparison, into = c("region", "timepoint"), sep = "_")

intersection_df
intersection_df$timepoint <- as.factor(intersection_df$timepoint)

ggplot(intersection_df, aes(x=timepoint, y=n_terms, color=region)) +
  geom_point(aes(pch=region)) +
  geom_line(aes(group=region)) +
  labs(title="FGSEA terms overlap between male and female",
       y = "Number of overlapping terms",
       color = "Region", shape = "Region") +
  theme_bw() +
  theme(axis.title.x=element_blank())
ggsave(plot=last_plot(), paste0(output_path, "intersection_terms/overlaping_terms_male_female.png"))
