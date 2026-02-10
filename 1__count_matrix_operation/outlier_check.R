## Dans cette partie je vérifie si les "outliers" identifiés par PCA et T-sne sont vraiment bizarres 

library(tidyverse)
library(readODS)
library(ggplot2)
library(DESeq2)

rm(list=ls())

setwd("/home/marinevernier/Documents/projets/")

#### PATHS ####
raw_counts_path <- "female_cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
coldata_path <- "female_cpid_multiregion/data/counts_m39_M32/coldata.ods"
plot.path <- "female_cpid_multiregion/graphs_results/1__count_matrix_operation/"
coverage.path <- "female_cpid_multiregion/data/counts_m39_M32/coverage_per_sample.tsv"


counts <- read_csv(raw_counts_path)
coldata <- read_ods(coldata_path)
coverage <- read_tsv(coverage.path)
coverage$sample <- gsub("-", ".", coverage$sample)
coverage$sample <- sub(
  "^([A-Za-z]+)([0-9]+)$",
  "\\1.\\2",
  coverage$sample
)


counts_long <- counts %>%
  pivot_longer(
    cols = -MGI.symbol,
    names_to = "sample",
    values_to = "expression"
  )

genes_present_testes <- c("Grp", "Drd2", "Drd1", "Camk2a", "Camk2b", "Mog", "Th", "Thy1")

counts <- counts %>% column_to_rownames("MGI.symbol")
"vmat" %in% rownames(counts)


Mog <- counts_long %>%
  filter(MGI.symbol == "Mog")
Mog_annot <- Mog %>%
  left_join(coldata, by = "sample")

ggplot(Mog_annot, aes(x = reg, y = expression, label = sample)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(vjust = -0.5, size = 3) +
  theme_classic() +
  labs(
    title = "Mog – détection d'échantillons aberrants",
    x = "Région",
    y = "Expression"
  )

plots <- list()

for (gene in genes_present_testes) {
  
  counts_gene <- counts_long %>%
    filter(MGI.symbol == gene)
  
  count_cov <- counts_gene %>%
    left_join(coldata, by = "sample") %>%
    left_join(coverage, by = "sample") %>%
    mutate(
      expression_covnorm = (expression / dedup_bam_reads) * 1e6
    )
  
  plots[[gene]] <- ggplot(
    count_cov,
    aes(x = reg, y = expression_covnorm, label = sample)
  ) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(vjust = -0.5, size = 3) +
    theme_classic() +
    labs(
      title = paste0(gene, " – expression normalisée par la couverture"),
      x = "Région",
      y = "CPM (dedup reads)"
    )
}
plots[["Thy1"]]



############################################################################################################

all(colnames(counts) == coldata$sample)

mean_intra_cor <- function(counts, coldata, region,
                            exclude_sample = NULL,
                            method = "pearson") {
  
  meta_sub <- coldata %>% filter(reg == region)
  
  if (!is.null(exclude_sample)) {
    meta_sub <- meta_sub %>%
      filter(!sample %in% exclude_sample)
  }
  
  if (nrow(meta_sub) < 3) return(NA_real_)
  
  counts_sub <- counts[, meta_sub$sample, drop = FALSE]
  
  cor_mat <- cor(counts_sub, method = method)
  
  mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
}



outlier_sample <- c("Ins.1837", "Hb.1839")

regions <- unique(coldata$reg)

cor_summary <- tibble(
  reg = regions,
  cor_with_outlier = purrr::map_dbl(
    regions,
    ~ mean_intra_cor(counts, coldata, .x, method = "spearman")
  ),
  cor_without_outlier = purrr::map_dbl(
    regions,
    ~ mean_intra_cor(counts, coldata, .x,
                     exclude_sample = outlier_sample,
                     method = "spearman")
  )
) %>%
  mutate(
    delta = cor_without_outlier - cor_with_outlier
  )

cor_summary

