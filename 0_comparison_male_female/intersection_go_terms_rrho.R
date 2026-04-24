library(tidyverse)
library(openxlsx)
library(dplyr)

rm(list=ls())

project_path <- "/home/marinevernier/Documents/projets/"
setwd(project_path)

## Path 
go_results_female.path <- "female_cpid_multiregion/data/3__RRHO2/rrvgo_overlaps_rrho/GO_summary.xlsx"
go_resultats_male.path <- "cpid_multiregion/data/3__RRHO2/rrvgo_overlaps_rrho/GO_summary.xlsx"
output_data_path <- "female_cpid_multiregion/data/0_comparison_male_female/overlaps_rrho/"
output_path <- "female_cpid_multiregion/graphs_results/0_comparison_male_female/overlaps_rrho/"

## Load 
data_female <- read.xlsx(go_results_female.path)
data_male <- read.xlsx(go_resultats_male.path)

# Filtre des comparaisons communes
common_comparisons <- intersect(
  unique(data_male$comparison),
  unique(data_female$comparison)
)

data_male   <- data_male   %>% filter(comparison %in% common_comparisons)
data_female <- data_female %>% filter(comparison %in% common_comparisons)

intersections_list <- lapply(common_comparisons, function(comp) {
  
  male_terms <- data_male %>%
    filter(comparison == comp) %>%
    pull(Description) %>%
    unique()
  
  female_terms <- data_female %>%
    filter(comparison == comp) %>%
    pull(Description) %>%
    unique()
  
  intersect(male_terms, female_terms)
})

# nommer la liste
names(intersections_list) <- common_comparisons

intersections_df <- bind_rows(
  lapply(names(intersections_list), function(comp) {
    
    terms <- intersections_list[[comp]]
    
    if (length(terms) == 0) return(NULL)
    
    data.frame(
      comparison = comp,
      Description = terms
    )
  })
)

intersections_full <- intersections_df %>%
  left_join(
    data_male %>%
      select(comparison, Description, p.adjust) %>%
      rename(padj_male = p.adjust),
    by = c("comparison", "Description")
  ) %>%
  left_join(
    data_female %>%
      select(comparison, Description, p.adjust) %>%
      rename(padj_female = p.adjust),
    by = c("comparison", "Description")
  )

sapply(intersections_list, length)

nb_terms_female <- data_female %>%
  group_by(comparison) %>%
  summarise(
    n_terms_female = n(),
    .groups = "drop"
  )

nb_terms_male <- data_male %>%
  group_by(comparison) %>%
  summarise(
    n_terms_male = n(),
    .groups = "drop"
  )

  
intersections_summary <- intersections_df %>%
  group_by(comparison) %>%
  summarise(
    Descriptions = paste0('"', Description, '"', collapse = ", "),
    n_terms_intersection = n()
  ) 


write.xlsx(intersections_full, file = paste0(output_data_path, "intersection_terms_overlap_rrho_male_vs_female.xlsx"))

write.xlsx(intersections_summary, file = paste0(output_data_path, "intersection_terms_overlap_rrho_male_vs_female_summary.xlsx"))




