#### Library ####
library(dplyr)
library(ggplot2)
library(tidytext)

rm(list=ls())

#### PATHS ####
project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

## INPUT ## 
Nac_modulePreservation_refMale.path <- "data/4__MEGENA/preservation_vs_male/Nac/Nac_modulePreservation_refMale.Rdata"
Nac_modulePreservation_refFemale.path <- "data/4__MEGENA/preservation_vs_male/Nac/Nac_modulePreservation_refFemale.Rdata"
Hb_modulePreservation_refMale.path <- "data/4__MEGENA/preservation_vs_male/Hb/Hb_modulePreservation_refMale.Rdata"
Hb_modulePreservation_refFemale.path <- "data/4__MEGENA/preservation_vs_male/Hb/Hb_modulePreservation_refFemale.Rdata"
Ins_modulePreservation_refMale.path <- "data/4__MEGENA/preservation_vs_male/Ins/Ins_modulePreservation_refMale.Rdata"
Ins_modulePreservation_refFemale.path <- "data/4__MEGENA/preservation_vs_male/Ins/Ins_modulePreservation_refFemale.Rdata"

## LOAD AN RENAME ## 
load(Nac_modulePreservation_refFemale.path)
Nac_modulePreservation_refFemale <- res_refFemale
load(Nac_modulePreservation_refMale.path)
Nac_modulePreservation_refMale <- res_refMale

load(Hb_modulePreservation_refFemale.path)
Hb_modulePreservation_refFemale <- res_refFemale
load(Hb_modulePreservation_refMale.path)
Hb_modulePreservation_refMale <- res_refMale

load(Ins_modulePreservation_refFemale.path)
Ins_modulePreservation_refFemale <- res_refFemale
load(Ins_modulePreservation_refMale.path)
Ins_modulePreservation_refMale <- res_refMale

add_metadata <- function(df, reference, region) {
  df %>%
    mutate(
      reference = reference,
      region = region
    )
}

Nac_female <- add_metadata(Nac_modulePreservation_refFemale, "female", "Nac")
Nac_male   <- add_metadata(Nac_modulePreservation_refMale, "male", "Nac")

Hb_female  <- add_metadata(Hb_modulePreservation_refFemale, "female", "Hb")
Hb_male    <- add_metadata(Hb_modulePreservation_refMale, "male", "Hb")

Ins_female  <- add_metadata(Ins_modulePreservation_refFemale, "female", "Ins")
Ins_male    <- add_metadata(Ins_modulePreservation_refMale, "male", "Ins")

all_pres <- bind_rows(
  Nac_female,
  Nac_male,
  Hb_female,
  Hb_male,
  Ins_female,
  Ins_male
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
  facet_wrap(~reference + region, scales = "free_x") +
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
  scale_fill_manual(values = c("female" = "lightpink", "male" = "lightblue"))
  

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
  scale_fill_manual(values = c("female" = "lightpink", "male" = "lightblue"))

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
  
  scale_color_manual(values = c("female" = "lightpink",
                                "male" = "lightblue"))

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
                                "Nac" = "lightgreen"))

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
  scale_color_manual(values = c("female" = "lightpink",
                                "male" = "lightblue"))
