#### VIOLIN PLOT PANREGION - script adapté pour étudié seulement ~reg + group ####

rm(list=ls())

#### LIBRAIRIES ####
library(tidyverse)
library(ggplot2)
library(dplyr)

#### PATHS ####
project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

## INPUT ## 
panregion_tp1.path <- "data/2__differential_expression_analysis/panregion_common_regions/deseq2/design_reg_group/annotation_deg/panregion_deg_tp1.rds"
panregion_tp2.path <- "data/2__differential_expression_analysis/panregion_common_regions/deseq2/design_reg_group/annotation_deg/panregion_deg_tp2.rds"
panregion_tp3.path <- "data/2__differential_expression_analysis/panregion_common_regions/deseq2/design_reg_group/annotation_deg/panregion_deg_tp3.rds"
deglist.path <- "data/2__differential_expression_analysis/deglist.Rdata"

## OUTPUT ## 
violin.path1 <- "graphs_results/panregion_common_regions/deseq2_reg_group/violin_panregion/violin_panregion_tp1.png"
violin.path2 <- "graphs_results/panregion_common_regions/deseq2_reg_group/violin_panregion/violin_panregion_tp2.png"
violin.path3 <- "graphs_results/panregion_common_regions/deseq2_reg_group/violin_panregion/violin_panregion_tp3.png"

#### DATA PROCESSING ####
panregion_tp1 <- readRDS(panregion_tp1.path)
panregion_tp2 <- readRDS(panregion_tp2.path)
panregion_tp3 <- readRDS(panregion_tp3.path)
load(deglist.path)


# TP1
panregion_tp1$region <- "panregion"
deg_Hb_tp1$region <- "Hb"
deg_Ins_tp1$region <- "Ins"
deg_Nac_tp1$region <- "Nac"

df_tp1 <- bind_rows(panregion_tp1, deg_Hb_tp1, deg_Ins_tp1, deg_Nac_tp1)



violin_plot_tp1 <- df_tp1 %>%
  ggplot(aes(x = region, y = log2fc, fill = region)) +
  geom_violin() +
  stat_summary(aes(x = region, y = log2fc, shape = "Mean"), fun = mean,
               geom = "point", size = 3, color = "black", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  stat_summary(aes(x = region, y = log2fc, shape = "Median"), fun = median,
               geom = "point", size = 3, color = "red", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  scale_shape_manual(name = "Stats", values = c("Mean" = 20, "Median" = 17)) +
  scale_fill_viridis_d(alpha = 0.3) +
  labs(title = "Comparison between panregion and individual regions at Timepoint 1",
       x = "Regions", y = "log2FoldChange", fill = "Regions") +
  theme_minimal() +
  scale_y_continuous(limits = c(-1.5, 1.65)) +
  theme(plot.title = element_text(size = 12)) 
violin_plot_tp1
ggsave(plot=last_plot(), violin.path1)


# TP2
panregion_tp2$region <- "panregion"
deg_Hb_tp2$region <- "Hb"
deg_Ins_tp2$region <- "Ins"
deg_Nac_tp2$region <- "Nac"

df_tp2 <- bind_rows(panregion_tp2, deg_Hb_tp2, deg_Ins_tp2, deg_Nac_tp2)


violin_plot_tp2 <- df_tp2 %>%
  ggplot(aes(x = region, y = log2fc, fill = region)) +
  geom_violin() +
  stat_summary(aes(x = region, y = log2fc, shape = "Mean"), fun = mean,
               geom = "point", size = 3, color = "black", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  stat_summary(aes(x = region, y = log2fc, shape = "Median"), fun = median,
               geom = "point", size = 3, color = "red", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  scale_shape_manual(name = "Stats", values = c("Mean" = 20, "Median" = 17)) +
  scale_fill_viridis_d(alpha = 0.3) +
  labs(title = "Comparison between panregion and individual regions at Timepoint 2",
       x = "Regions", y = "log2FoldChange", fill = "Regions") +
  theme_minimal() +
  scale_y_continuous(limits = c(-1.05, 1.3)) +
  theme(plot.title = element_text(size = 12)) 
violin_plot_tp2
ggsave(plot=last_plot(), violin.path2)


# TP3
panregion_tp3$region <- "panregion"
deg_Hb_tp3$region <- "Hb"
deg_Ins_tp3$region <- "Ins"
deg_Nac_tp3$region <- "Nac"

df_tp3 <- bind_rows(panregion_tp3, deg_Hb_tp3, deg_Ins_tp3, deg_Nac_tp3)



violin_plot_tp3 <- df_tp3 %>%
  ggplot(aes(x = region, y = log2fc, fill = region)) +
  geom_violin() +
  stat_summary(aes(x = region, y = log2fc, shape = "Mean"), fun = mean,
               geom = "point", size = 3, color = "black", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  stat_summary(aes(x = region, y = log2fc, shape = "Median"), fun = median,
               geom = "point", size = 3, color = "red", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  scale_shape_manual(name = "Stats", values = c("Mean" = 20, "Median" = 17)) +
  scale_fill_viridis_d(alpha = 0.3) +
  labs(title = "Comparison between panregion and individual regions at Timepoint 3",
       x = "Regions", y = "log2FoldChange", fill = "Regions") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.8, 1.1)) +
  theme(plot.title = element_text(size = 12)) 
violin_plot_tp3
ggsave(plot=last_plot(), violin.path3)
