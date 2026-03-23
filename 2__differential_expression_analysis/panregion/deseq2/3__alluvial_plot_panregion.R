#### ALLUVIAL PLOT FOR PANREGION (~reg+group) THROUGH TIMEPOINTS ####

rm(list=ls())

#### LIBRARIES ####
library(tidyverse)
library(ggalluvial)

#### PATHS ####

project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

## INPUT ## 
df_tp1.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_tp1.rds"
df_tp2.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_tp2.rds"
df_tp3.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_annotated_genes_tp3.rds"

## PALETTE ##  
palette.path <- "data/2__differential_expression_analysis/alluvial_palette.rds" # recupere la palette utilisé en region unique 

## OUTPUT ##
alluvial_plot.path <- "graphs_results/panregion/deseq2_reg_group/alluvial_panregion_reg_group.png" 

#### FORMATTING DATA ####

df_tp1<- readRDS(df_tp1.path)
df_tp2<- readRDS(df_tp2.path)
df_tp3<- readRDS(df_tp3.path)

df_tp1 <- df_tp1 %>% select(label, TP1 = diffexpressed)
df_tp2 <- df_tp2 %>% select(label, TP2 = diffexpressed) 
df_tp3 <- df_tp3 %>% select(label, TP3 = diffexpressed) 

# Join all tp 
df_alltp <- df_tp1 %>%
  full_join(df_tp2, by = "label") %>%
  full_join(df_tp3, by = "label")


#### GET THE NUMBER OF GENES PER PATTERN ####
## SORT ALL GENES BY PATTERN ##
df_pattern <- df_alltp %>%
  group_by(label) %>%
  mutate(diffexpressed_alltp = paste(TP1, TP2, TP3, sep = "_")) %>%
  ungroup() %>% 
  group_by(diffexpressed_alltp) %>%  
  summarise(n = n())

## REMOVE NS_NS_NS PATTERN ##
df_pattern <- df_pattern %>%
  filter(diffexpressed_alltp != "ns_ns_ns")

## TP SEPARATION FOR ALLUVIAL ##
df_pattern_sep<- df_pattern %>%
  separate(diffexpressed_alltp, into = c("TP1", "TP2", "TP3"), sep = "_")

df_pattern_sep <- df_pattern_sep %>%
  mutate(
    TP1 = factor(TP1, levels = c("UP", "ns", "DOWN")),  
    TP2 = factor(TP2, levels = c("UP", "ns", "DOWN")),
    TP3 = factor(TP3, levels = c("UP", "ns", "DOWN"))
  )

#### GET PALETTE ####
palette <- readRDS(palette.path)

# get only present patters
pal <- palette[names(palette) %in% df_pattern$diffexpressed_alltp]


#### ALLUVIAL PLOT ####
alluvial_panregion <-  ggplot(df_pattern_sep,
                              aes(y = n, axis1 = TP1, axis2 = TP2, axis3 = TP3)) +
  geom_alluvium(aes(fill = df_pattern$diffexpressed_alltp), alpha = 0.6) + 
  geom_stratum(alpha = .25, width = 1/4, reverse = TRUE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  labs(y = "Number of genes", x = "Timepoints",
       title = "Alluvial plot of genes expression across timepoints (~reg+RIN+group)") +
  guides(fill = "none") +
  scale_x_discrete(limits = c("TP1", "TP2", "TP3")) +
  scale_fill_manual(values=unname(pal)) +
  theme_bw()

alluvial_panregion

#### SAVE ####
ggsave(alluvial_plot.path, plot=alluvial_panregion, width=1900, height=1200,
       units="px", scale=2)
