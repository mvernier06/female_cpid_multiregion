#### Library ####
library(dplyr)
library(ggplot2)
library(tidytext)

rm(list=ls())

#### PATHS ####
project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

## INPUT ## 
ACC_refWith.path <- 'data/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN/ACC/ACC_modulePreservation_refWith.Rdata'
ACC_refWithout.path <- 'data/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN/ACC/ACC_modulePreservation_refwithout.Rdata'
Hb_refWith.path <- 'data/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN/Hb/Hb_modulePreservation_refWith.Rdata'
Hb_refWithout.path <- 'data/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN/Hb/Hb_modulePreservation_refwithout.Rdata'
Ins_refWith.path <- 'data/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN/Ins/Ins_modulePreservation_refWith.Rdata'
Ins_refWithout.path <- 'data/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN/Ins/Ins_modulePreservation_refwithout.Rdata'
Nac_refWith.path <- 'data/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN/Nac/Nac_modulePreservation_refWith.Rdata'
Nac_refWithout.path <- 'data/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN/Nac/Nac_modulePreservation_refwithout.Rdata'


## LOAD AN RENAME ## 
load(ACC_refWith.path)
ACC_refWith <- res_refWith
load(ACC_refWithout.path)
ACC_refWithout <- res_refwithout

load(Hb_refWith.path)
Hb_refWith <- res_refWith
load(Hb_refWithout.path)
Hb_refWithout <- res_refwithout

load(Ins_refWith.path)
Ins_refWith <- res_refWith
load(Ins_refWithout.path)
Ins_refWithout <- res_refwithout

load(Nac_refWith.path)
Nac_refWith <- res_refWith
load(Nac_refWithout.path)
Nac_refWithout <- res_refwithout

# Output
plot.path <- "graphs_results/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN_correction/"

add_metadata <- function(df, reference, region) {
  df %>%
    mutate(
      reference = reference,
      region = region
    )
}

ACC_with <- add_metadata(ACC_refWith, "with", region = "ACC")
ACC_without <- add_metadata(ACC_refWithout, "without", region="ACC")

Hb_with <- add_metadata(Hb_refWith, "with", region = "Hb")
Hb_without <- add_metadata(Hb_refWithout, "without", region="Hb")

Ins_with <- add_metadata(Ins_refWith, "with", region = "Ins")
Ins_without <- add_metadata(Ins_refWithout, "without", region="Ins")

Nac_with <- add_metadata(Nac_refWith, "with", region = "Nac")
Nac_without <- add_metadata(Nac_refWithout, "without", region="Nac")

all_pres <- bind_rows(
  ACC_with, 
  ACC_without, 
  Hb_with, 
  Hb_without, 
  Ins_with, 
  Ins_without, 
  Nac_with, 
  Nac_without
)


ggplot(all_pres, aes(x = module, y = Zsummary, fill = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~reference, scales = "free_x" ) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Module preservation (Zsummary)",
       y = "Zsummary",
       x = "Module")

ggplot(all_pres,
       aes(x = reorder_within(module, Zsummary, reference),
           y = Zsummary,
           fill = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(reference~ region, scales = "free_x") +
  scale_x_reordered() +
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "Module preservation (Zsummary)",
       y = "Zsummary",
       x = "Module")


categorize <- function(z) {
  case_when(
    z > 10 ~ "high",
    z > 2 ~ "moderate",
    TRUE ~ "low"
  )
}

all_pres <- all_pres %>%
  mutate(category = categorize(Zsummary))

summary_stats <- all_pres %>%
  group_by(region, reference) %>%
  summarise(
    n_modules = n(),
    high = sum(category == "high"),
    moderate = sum(category == "moderate"),
    low = sum(category == "low")
  )
print(summary_stats)


plot_data <- all_pres %>%
  group_by(region, reference, category) %>%
  summarise(n = n(), .groups = "drop")
plot_data$category <- factor(plot_data$category,
                             levels = c("low", "moderate", "high"))
totals <- plot_data %>%
  group_by(region, reference) %>%
  summarise(total = sum(n), .groups = "drop")

ggplot(plot_data, aes(x = category, y = n, fill = reference)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3,
            size = 3) +
  
  facet_wrap(~region) +
  theme_bw() +
  labs(title = "Module preservation categories",
       x = "Category",
       y = "Number of modules",
       fill = "Reference") +
  scale_fill_manual(values = c("with" = "lightpink", "without" = "lightblue"))
file_name <- paste0(plot.path, "module_preservation_by_categories.png")
ggsave(plot=last_plot(), file_name)

plot_data_prop <- plot_data %>%
  group_by(region, reference) %>%
  mutate(prop = n / sum(n))

ggplot(plot_data_prop, aes(x = category, y = prop, fill = reference)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  
  geom_text(aes(label = round(prop,3)),
            position = position_dodge(width = 0.9),
            vjust = -0.3,
            size = 3) +
  
  facet_wrap(~region) +
  theme_bw() +
  labs(title = "Module preservation categories, in proportions",
       x = "Category",
       y = "Proportion",
       fill = "Reference") +
  scale_fill_manual(values = c("with" = "lightpink", "without" = "lightblue"))
file_name <- paste0(plot.path, "module_preservation_by_categories_proportions.png")
ggsave(plot=last_plot(), file_name)


ggplot(all_pres, aes(x = size, y = Zsummary, color = reference)) +
  geom_point(alpha = 0.5, size = 2) +
  
  facet_wrap(~region) + #  scales = "free"
  
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  
  theme_bw() +
  labs(title = "Module preservation vs module size",
       x = "Module size",
       y = "Zsummary",
       color = "Reference") +
  
  scale_color_manual(values = c("with" = "lightpink",
                                "without" = "lightblue"))
file_name <- paste0(plot.path, "module_preservation_module_size.png")
ggsave(plot=last_plot(), file_name)

ggplot(all_pres, aes(x = size, y = Zsummary, color = region)) +
  geom_point(alpha = 0.5, size = 2) +
  
  facet_wrap(~reference) + #  scales = "free"
  
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  
  theme_bw() +
  labs(title = "Module preservation vs module size",
       x = "Module size",
       y = "Zsummary",
       color = "Regions") +
  
  scale_color_manual(values = c("Ins" = "lightpink",
                                "Hb" = "lightblue",
                                "Nac" = "lightgreen", 
                                "ACC" = "yellow"))
file_name <- paste0(plot.path, "module_preservation_module_size_reg.png")
ggsave(plot=last_plot(), file_name)



ggplot(all_pres, aes(x = size, y = Zsummary, color = reference)) +
  geom_point(alpha = 0.4, size = 2) +
  facet_wrap(~region) +
  coord_cartesian(xlim = c(0, 200)) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "Module preservation vs module size",
       x = "Module size",
       y = "Zsummary",
       color = "Reference") +
  scale_color_manual(values = c("with" = "lightpink",
                                "without" = "lightblue"))
file_name <- paste0(plot.path, "module_preservation_module_size_zoom.png")
ggsave(plot=last_plot(), file_name)
