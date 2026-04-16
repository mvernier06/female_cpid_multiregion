library(tidyverse)
library(ggplot2)
library(ggalluvial)

rm(list=ls())

#### PATHS ####
project.path <- "/home/marinevernier/Documents/projets/"
palette_alluvial.path <- "female_cpid_multiregion/data/2__differential_expression_analysis/alluvial_palette.rds"

plot.path <- "female_cpid_multiregion/graphs_results/0_comparison_male_female/alluvial/"
dir.create(plot.path)

setwd(project.path)

pal <- readRDS(palette_alluvial.path)

#### number of DEGs per alluvial pattern for each reg ####
get_pattern_proportions <- function(regionlist, sex, repository){
  
  res_list <- list()
  
  for(reg in regionlist){
    
    alluvial_patterns.path <- paste0(repository,"alluvial_patterns_", reg, ".Rdata")
    load(alluvial_patterns.path)
    
    df <- df_new %>%
      dplyr::select(diffexpressed_alltp) %>%
      group_by(diffexpressed_alltp) %>%
      mutate(n = n())
    
    df <- get("df_new") %>% 
      dplyr::group_by(diffexpressed_tp1, diffexpressed_tp2, diffexpressed_tp3,
                      diffexpressed_alltp) %>% 
      summarise(n = n(), .groups="drop") %>%
      mutate(region = reg,
             sex = sex)
    
    res_list[[reg]] <- df
  }
  
  bind_rows(res_list)
}

df_male <- get_pattern_proportions(c("Hb", "Ins", "NAc"), "male", "cpid_multiregion/data/2__differential_expression_analysis/")
df_female <- get_pattern_proportions(c("Hb", "Ins", "Nac"), "female", "female_cpid_multiregion/data/2__differential_expression_analysis/")

df_female$region <- dplyr::recode(df_female$region,
                                  "Nac" = "NAc")

df_all <- bind_rows(df_male, df_female)
df_all <- df_all %>%
  group_by(region, sex) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()



pattern_order <- df_all %>%
  arrange(desc(prop)) %>%
  pull(diffexpressed_alltp)

df_all$diffexpressed_alltp <- factor(df_all$diffexpressed_alltp,
                                     levels = pattern_order)
df <- df_all
reg <- 'Ins'
sex
plot_bar_patterns <- function(df, reg){
  
  df_reg <- df %>% 
    dplyr::filter(region == reg)
  
  pal_reg <- pal[names(pal) %in% levels(df_reg$diffexpressed_alltp)]
  
  ggplot(df_reg,
         aes(x = sex,
             y = prop,
             fill = diffexpressed_alltp)) +
    
    geom_bar(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.2) +
    
    # annotations
    geom_text(aes(label = ifelse(prop > 0.04,
                                 paste0(diffexpressed_alltp, " : ", round(prop, 2)),
                                 "")),
              position = position_stack(vjust = 0.5),
              size = 3) +
    
    scale_fill_manual(values = pal_reg) +
    labs(y = "Proportion", x = "") +
    ggtitle(paste("Pattern proportions -", reg)) +
    theme_bw()
}

p <- plot_bar_patterns(df_all, "Ins")
print(p)
ggsave(p, file = paste0(plot.path, "pattern_prop_Ins.png"))
p <- plot_bar_patterns(df_all, "Hb")
print(p)
ggsave(p, file = paste0(plot.path, "pattern_prop_Hb.png"))
p <- plot_bar_patterns(df_all, "NAc")
print(p)
ggsave(p, file = paste0(plot.path, "pattern_prop_NAc.png"))


df_prop <- df_all %>%
  dplyr::select(region, sex, diffexpressed_alltp, prop)
df_wide <- df_prop %>%
  pivot_wider(names_from = sex, values_from = prop, values_fill = 0)
df_flow <- df_wide %>%
  dplyr::select(region, pattern_male = diffexpressed_alltp, male) %>%
  dplyr::left_join(
    df_wide %>%
      dplyr::select(region, pattern_female = diffexpressed_alltp, female),
    by = "region"
  ) %>%
  dplyr::mutate(flow = male * female)

plot_alluvial_sex <- function(df, reg){
  
  df_reg <- df %>% dplyr::filter(region == reg)
  
  ggplot(df_reg,
         aes(axis1 = pattern_male,
             axis2 = pattern_female,
             y = flow)) +
    
    # FLOWS (liens) → autre couleur
    geom_alluvium(fill = "grey70", alpha = 0.5) +
    
    # STRATA (barres) → palette principale
    geom_stratum(aes(fill = after_stat(stratum)),
                 width = 1/5,
                 color = "black") +
    
    # geom_text(stat = "stratum",
    #           aes(label = after_stat(stratum))) +
    
    scale_x_discrete(limits = c("Male", "Female")) +
    
    # 🎨 palette des patterns
    scale_fill_manual(values = pal) +
    
    ggtitle(paste("Pattern redistribution:", reg)) +
    theme_bw()
}

p <- plot_alluvial_sex(df_flow, "NAc")
print(p)
