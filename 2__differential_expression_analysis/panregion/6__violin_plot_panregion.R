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
rg_panregion_tp1.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_tp1.rds"
rg_panregion_tp2.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_tp2.rds"
rg_panregion_tp3.path <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg/panregion_deg_tp3.rds"

## OUTPUT ## 
violin.path1 <- "graphs_results/panregion/deseq2_reg_group/violin_panregion/violin_panregion_tp1.png"
violin.path2 <- "graphs_results/panregion/deseq2_reg_group/violin_panregion/violin_panregion_tp2.png"
violin.path3 <- "graphs_results/panregion/deseq2_reg_group/violin_panregion/violin_panregion_tp3.png"

#### DATA PROCESSING ####
rg_panregion_tp1 <- readRDS(rg_panregion_tp1.path)
rg_panregion_tp2 <- readRDS(rg_panregion_tp2.path)
rg_panregion_tp3 <- readRDS(rg_panregion_tp3.path)


# TP1 
stat_tests_tp1 <- tibble(
  t_test_p = t.test(rg_panregion_tp1$log2fc, mu = 0)$p.value,
  wilcox_p = wilcox.test(rg_panregion_tp1$log2fc, mu = 0)$p.value
)

stat_tests_tp1 <- stat_tests_tp1 %>%
  mutate(
    significance_t_test = case_when(
      t_test_p < 0.001 ~ "***",
      t_test_p < 0.01 ~ "**",
      t_test_p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    significance_wilcox = case_when(
      wilcox_p < 0.001 ~ "***",
      wilcox_p < 0.01 ~ "**",
      wilcox_p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

ggplot(rg_panregion_tp1, aes(x = "", y = log2fc)) +
  geom_violin(fill = "lightblue") +
  stat_summary(aes(x = "", y = log2fc, shape = "Mean"), fun = mean,
               geom = "point", size = 3, color = "black", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  stat_summary(aes(x = "", y = log2fc, shape = "Median"), fun = median,
               geom = "point", size = 3, color = "red", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  scale_shape_manual(name = "Stats", values = c("Mean" = 20, "Median" = 17)) +
  scale_fill_viridis_d(alpha = 0.3) +
  labs(title = "Distribution log2FC (reg_group TP1)",
       y = "log2FoldChange", x = "") +
  scale_y_continuous(limits = c(-0.5, 0.6)) +
  theme(plot.title = element_text(size = 12)) + 
  geom_text(data = stat_tests_tp1, aes (x="", y = max(rg_panregion_tp1$log2fc) + 0.1,
                                        label = significance_t_test), color = "black") +
  geom_text(data = stat_tests_tp1, aes (x="", y = max(rg_panregion_tp1$log2fc) + 0.15,
                                        label = significance_wilcox), color = "red")
ggsave(plot=last_plot(), violin.path1)


# TP2
stat_tests_tp2 <- tibble(
  t_test_p = t.test(rg_panregion_tp2$log2fc, mu = 0)$p.value,
  wilcox_p = wilcox.test(rg_panregion_tp2$log2fc, mu = 0)$p.value
)

stat_tests_tp2 <- stat_tests_tp2 %>%
  mutate(
    significance_t_test = case_when(
      t_test_p < 0.001 ~ "***",
      t_test_p < 0.01 ~ "**",
      t_test_p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    significance_wilcox = case_when(
      wilcox_p < 0.001 ~ "***",
      wilcox_p < 0.01 ~ "**",
      wilcox_p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

ggplot(rg_panregion_tp2, aes(x = "", y = log2fc)) +
  geom_violin(fill = "lightblue") +
  stat_summary(aes(x = "", y = log2fc, shape = "Mean"), fun = mean,
               geom = "point", size = 3, color = "black", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  stat_summary(aes(x = "", y = log2fc, shape = "Median"), fun = median,
               geom = "point", size = 3, color = "red", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  scale_shape_manual(name = "Stats", values = c("Mean" = 20, "Median" = 17)) +
  scale_fill_viridis_d(alpha = 0.3) +
  labs(title = "Distribution log2FC (reg_group TP2)",
       y = "log2FoldChange", x = "") +
  scale_y_continuous(limits = c(-0.5, 0.6)) +
  theme(plot.title = element_text(size = 12)) + 
  geom_text(data = stat_tests_tp2, aes (x="", y = max(rg_panregion_tp2$log2fc) + 0.1,
                                        label = significance_t_test), color = "black") +
  geom_text(data = stat_tests_tp2, aes (x="", y = max(rg_panregion_tp2$log2fc) + 0.15,
                                        label = significance_wilcox), color = "red")
ggsave(plot=last_plot(), violin.path2)


# TP3
stat_tests_tp3 <- tibble(
  t_test_p = t.test(rg_panregion_tp3$log2fc, mu = 0)$p.value,
  wilcox_p = wilcox.test(rg_panregion_tp3$log2fc, mu = 0)$p.value
)

stat_tests_tp3 <- stat_tests_tp3 %>%
  mutate(
    significance_t_test = case_when(
      t_test_p < 0.001 ~ "***",
      t_test_p < 0.01 ~ "**",
      t_test_p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    significance_wilcox = case_when(
      wilcox_p < 0.001 ~ "***",
      wilcox_p < 0.01 ~ "**",
      wilcox_p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

ggplot(rg_panregion_tp3, aes(x = "", y = log2fc)) +
  geom_violin(fill = "lightblue") +
  stat_summary(aes(x = "", y = log2fc, shape = "Mean"), fun = mean,
               geom = "point", size = 3, color = "black", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  stat_summary(aes(x = "", y = log2fc, shape = "Median"), fun = median,
               geom = "point", size = 3, color = "red", position = position_dodge(0.9),
               inherit.aes = FALSE) +
  scale_shape_manual(name = "Stats", values = c("Mean" = 20, "Median" = 17)) +
  scale_fill_viridis_d(alpha = 0.3) +
  labs(title = "Distribution log2FC (reg_group TP3)",
       y = "log2FoldChange", x = "") +
  scale_y_continuous(limits = c(-0.5, 0.6)) +
  theme(plot.title = element_text(size = 12)) + 
  geom_text(data = stat_tests_tp3, aes (x="", y = max(rg_panregion_tp3$log2fc) + 0.1,
                                        label = significance_t_test), color = "black") +
  geom_text(data = stat_tests_tp3, aes (x="", y = max(rg_panregion_tp3$log2fc) + 0.15,
                                        label = significance_wilcox), color = "red")
ggsave(plot=last_plot(), violin.path3)
