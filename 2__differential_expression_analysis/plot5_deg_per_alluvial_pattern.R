library(tidyverse)
library(ggplot2)

rm(list=ls())

#### PATHS ####
plot.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/deg_per_alluvial_pattern"
dir.create(plot.path)

setwd(plot.path)
regionlist <- c("ACC", "Hb", "Ins", "Nac")

plots <- list()
i=1
#### number of DEGs per alluvial pattern for each reg ####
for(reg in regionlist){
  print(reg)
  
  alluvial_patterns.path <- paste0("/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/alluvial_patterns_", reg, ".Rdata")
  load(alluvial_patterns.path)
  
  df <- df_new %>%
    dplyr::select(diffexpressed_alltp) %>%
    group_by(diffexpressed_alltp) %>%
    mutate(n = n())
  
  p <- df |>
    ggplot(aes(fct_infreq(diffexpressed_alltp))) +
    geom_bar() +
    geom_label(aes(y = n, label = n), data = df) +
    labs(x = "Alluvial patterns",
         y = "Number of genes",
         title = paste0("Number of genes in the ", reg, " per alluvial pattern")) +
    theme(axis.text.x = element_text(angle = 30, vjust=1, hjust=1))
  ggsave(plot=p, filename=paste0("deg_per_alluvial_pattern_", reg, ".png"),
         height = 2000, width = 3000, units = "px")
  plots[i] <- p
  i=i+1
}
plots[4]
