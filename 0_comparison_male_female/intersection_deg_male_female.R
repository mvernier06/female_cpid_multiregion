library(tidyverse)
library(tibble)
library(ggplot2)

rm(list=ls())

project.path <- "/home/marinevernier/Documents/projets/"
setwd(project.path)

deg_male_path <- "cpid_multiregion/data/2__differential_expression_analysis/deglist.Rdata"
deg_female_path <- "female_cpid_multiregion/data/2__differential_expression_analysis/deglist.Rdata"

output_path <- "female_cpid_multiregion/graphs_results/0_comparison_male_female/intersection_deg/"
dir.create(output_path)

# -------- LOAD MALE --------
male_env <- new.env()
load(deg_male_path, envir = male_env)

male_objs <- ls(male_env, pattern = "^deg_")

for(obj in male_objs){
  assign(paste0(obj, "_male"), male_env[[obj]])
}


# -------- LOAD FEMALE --------
female_env <- new.env()
load(deg_female_path, envir = female_env)

female_objs <- ls(female_env, pattern = "^deg_")

for(obj in female_objs){
  assign(paste0(obj, "_female"), female_env[[obj]])
}

deg_NAc_tp1_female <- deg_Nac_tp1_female
deg_NAc_tp2_female <- deg_Nac_tp2_female
deg_NAc_tp3_female <- deg_Nac_tp3_female

regions <- c("Ins", "NAc", "Hb")
timepoints <- c("tp1", "tp2", "tp3")

get_deg_genes <- function(df, pval_threshold = 0.05, lfc_thresh = 0){
  df %>%
    filter(pval < pval_threshold & abs(log2fc) > lfc_thresh) %>%
    pull(label)
}

intersections_list <- list()

for(tp in timepoints){
  for(reg in regions){
    
    male_name <- paste0("deg_", reg, "_", tp, "_male")
    female_name <- paste0("deg_", reg, "_", tp, "_female")
    
    if(exists(male_name) & exists(female_name)){
      
      male_df <- get(male_name)
      female_df <- get(female_name)
      
      male_genes <- get_deg_genes(male_df)
      female_genes <- get_deg_genes(female_df)
      
      inter <- intersect(male_genes, female_genes)
      
      intersections_list[[paste(reg, tp, sep = "_")]] <- inter
      
    }
  }
}
intersection_sizes <- sapply(intersections_list, length)
intersection_sizes



intersection_df <- tibble(
  comparison = names(intersections_list),
  n_genes = intersection_sizes
) %>%
  tidyr::separate(comparison, into = c("region", "timepoint"), sep = "_")

intersection_df
intersection_df$timepoint <- as.factor(intersection_df$timepoint)
  
ggplot(intersection_df, aes(x=timepoint, y=n_genes, color=region)) +
  geom_point(aes(pch=region)) +
  geom_line(aes(group=region)) +
  labs(title="DEG overlap between male and female",
       y = "Number of overlapping DEG",
       color = "Region", shape = "Region") +
  theme_bw() +
  theme(axis.title.x=element_blank())
ggsave(plot=last_plot(), paste0(output_path, "overlaping_deg_male_female.png"))
