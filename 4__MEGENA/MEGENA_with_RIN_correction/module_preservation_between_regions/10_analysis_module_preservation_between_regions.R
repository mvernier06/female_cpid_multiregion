#### Library ####
library(dplyr)
library(ggplot2)
library(tidytext)

rm(list=ls())

#### PATHS ####
project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

## INPUT ## 
ACC_vs_Hb.path <- "data/4__MEGENA/MEGENA_with_RIN_correction/preservation_between_regions/ACC_vs_Hb.Rdata"
ACC_vs_Ins.path <- "data/4__MEGENA/MEGENA_with_RIN_correction/preservation_between_regions/ACC_vs_Ins.Rdata"
ACC_vs_Nac.path <- "data/4__MEGENA/MEGENA_with_RIN_correction/preservation_between_regions/ACC_vs_Nac.Rdata"
Hb_vs_Nac.path <- "data/4__MEGENA/MEGENA_with_RIN_correction/preservation_between_regions/Hb_vs_Nac.Rdata"
Ins_vs_Hb.path <- "data/4__MEGENA/MEGENA_with_RIN_correction/preservation_between_regions/Ins_vs_Hb.Rdata"
Ins_vs_Nac.path <- "data/4__MEGENA/MEGENA_with_RIN_correction/preservation_between_regions/Ins_vs_Nac.Rdata"

## LOAD AN RENAME ## 
load(ACC_vs_Hb.path)
ACC_vs_Hb <- df
load(ACC_vs_Ins.path)
ACC_vs_Ins <- df
load(ACC_vs_Nac.path)
ACC_vs_Nac <- df
load(Hb_vs_Nac.path)
Hb_vs_Nac <- df
load(Ins_vs_Hb.path)
Ins_vs_Hb <- df
load(Ins_vs_Nac.path)
Ins_vs_Nac <- df

# Output
plot.path <- "graphs_results/4__MEGENA/MEGENA_with_RIN_correction/preservation_between_regions/"


all_res <- bind_rows(
  ACC_vs_Hb,
  ACC_vs_Ins,
  ACC_vs_Nac,
  Hb_vs_Nac,
  Ins_vs_Hb,
  Ins_vs_Nac
)

all_res <- all_res %>%
  mutate(
    pair = paste0(pmin(ref, test), "_vs_", pmax(ref, test)),
    direction = paste0(ref, "_to_", test)
  )

categorize <- function(z) {
  case_when(
    z > 10 ~ "high",
    z > 2 ~ "moderate",
    TRUE ~ "low"
  )
}

all_res <- all_res %>%
  mutate(category = categorize(Zsummary))

summary_stats <- all_res %>%
  group_by(pair, direction) %>%
  summarise(
    n_modules = n(),
    high = sum(category == "high"),
    moderate = sum(category == "moderate"),
    low = sum(category == "low"),
    .groups = "drop"
  )

print(summary_stats)

plot_data <- all_res %>%
  group_by(pair, direction, category) %>%
  summarise(n = n(), .groups = "drop")

plot_data$category <- factor(plot_data$category,
                             levels = c("low", "moderate", "high"))

ggplot(plot_data, aes(x = category, y = n, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3,
            size = 3) +
  
  facet_wrap(~pair) +
  theme_bw() +
  labs(title = "Module preservation by direction",
       x = "Category",
       y = "Number of modules",
       fill = "Direction")

plot_data_prop <- plot_data %>%
  group_by(pair, direction) %>%   
  mutate(prop = n / sum(n))

ggplot(plot_data_prop, aes(x = category, y = prop, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  
  geom_text(aes(label = round(prop, 3)),
            position = position_dodge(width = 0.9),
            vjust = -0.3,
            size = 3) +
  
  facet_wrap(~pair) +
  theme_bw() +
  labs(title = "Module preservation by direction, in proportions",
       x = "Category",
       y = "Number of modules",
       fill = "Direction")
