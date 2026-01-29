library(tidyverse)
library(readODS)
library(ggplot2)
library(DESeq2)

rm(list=ls())

#### PATHS ####
raw_counts_path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
coldata_path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/count_data/coldata.ods"
plot.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/graphs_results/1__count_matrix_operation/"



counts <- read_csv(raw_counts_path)
coldata <- read_ods(coldata_path)


counts_long <- counts %>%
  pivot_longer(
    cols = -MGI.symbol,              
    names_to = "sample",
    values_to = "expression"
  )


counts <- counts %>% column_to_rownames("MGI.symbol")
"Ucma" %in% rownames(counts)

Camk2a <- counts_long %>%
  filter(MGI.symbol == "Camk2a")
Camk2a_annot <- Camk2a %>%
  left_join(coldata, by = "sample")

ggplot(Camk2a_annot, aes(x = reg, y = expression, label = sample)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(vjust = -0.5, size = 3) +
  theme_classic() +
  labs(
    title = "Camk2a – détection d'échantillons aberrants",
    x = "Région",
    y = "Expression"
  )

ggplot(Camk2a_annot, aes(x = reg, y = log2(expression + 1), label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_classic()

Camk2a_annot %>%
  arrange(expression) %>%
  mutate(sample_ord = factor(sample, levels = sample)) %>%
  ggplot(aes(x = sample_ord, y = log2(expression + 1), color = reg)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  labs(
    x = "Sample",
    y = "log2(Camk2a + 1)",
    color = "Région"
  )


Drd2 <- counts_long %>%
  filter(MGI.symbol == "Drd2")
Drd2_annot <- Drd2 %>%
  left_join(coldata, by = "sample")

ggplot(Drd2_annot, aes(x = reg, y = expression, label = sample)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(vjust = -0.5, size = 3) +
  theme_classic() +
  labs(
    title = "Drd2 – détection d'échantillons aberrants",
    x = "Région",
    y = "Expression"
  )
ggplot(Drd2_annot, aes(x = reg, y = log2(expression + 1), label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_classic()


genes_glut <- c("Slc17a7", "Slc17a6", "Camk2a", "Camk2b", "Ntng2")

signature <- counts_long %>%
  filter(MGI.symbol %in% genes_glut) %>%
  group_by(sample) %>%
  summarise(mean_expr = mean(log2(expression + 1)))

signature %>%
  left_join(coldata, by = "sample") %>%
  ggplot(aes(x = reg, y = mean_expr, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_classic()







raw_counts <- read_csv("/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/CPID_sham_vs_cuff_betaprior.csv")
raw_counts$Geneid <-  str_replace(raw_counts$Geneid, "\\..*", "") 
"ENSG00000104888" %in% raw_counts






