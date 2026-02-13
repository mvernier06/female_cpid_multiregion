library(tidyverse)
library(readODS)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(M3C)
library(ggpubr)

rm(list=ls())

setwd("/home/marinevernier/Documents/projets/")

#### PATHS ####
raw_counts_path <- "female_cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
coldata_path <- "female_cpid_multiregion/data/counts_m39_M32/coldata.ods"
plot.path <- "female_cpid_multiregion/graphs_results/1__count_matrix_operation/quality_check/"
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

setwd(plot.path)

rownames(counts) <- counts$MGI.symbol
counts$MGI.symbol <- NULL

## normalization because PCA shows distances
vst_counts <- varianceStabilizingTransformation(as.matrix(counts))

# PCA classique
pca_res <- prcomp(t(vst_counts), scale. = TRUE)

# Scores des samples
scores <- as.data.frame(pca_res$x[, 1:10])  # PC1 à PC10
scores$sample <- rownames(scores)

# Ajouter RIN
scores <- scores %>%
  left_join(coldata[, c("sample", "RIN")], by = "sample")

cor_results <- lapply(scores[, 1:10], function(pc) {
  cor.test(pc, scores$RIN, method = "pearson")
})

cor_results

cor_df <- data.frame(
  PC = paste0("PC", 1:10),
  correlation = sapply(cor_results, function(x) x$estimate),
  p_value = sapply(cor_results, function(x) x$p.value)
)
cor_df$padj <- p.adjust(cor_df$p_value, method = "BH")
cor_df

cor_df$PC <- factor(cor_df$PC,
                    levels = paste0("PC", 1:10))
cor_df$significant <- cor_df$padj < 0.05

ggplot(cor_df, aes(x = PC, y = correlation, fill = significant)) +
  geom_col() +
  geom_text(aes(label = signif(padj, 2)), vjust = -0.5) +
  scale_fill_manual(values = c("FALSE" = "grey70",
                               "TRUE"  = "red")) +
  labs(
    title = "Pearson correlation between PCA axes and RIN (FDR adjusted)",
    y = "Pearson correlation",
    x = "Principal Components"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(plot=last_plot(), "pearson_correlation_RIN_PCA.png")


cor_results <- lapply(scores[, 1:10], function(pc) {
  cor.test(pc, scores$RIN, method = "spearman")
})

cor_results

cor_df <- data.frame(
  PC = paste0("PC", 1:10),
  correlation = sapply(cor_results, function(x) x$estimate),
  p_value = sapply(cor_results, function(x) x$p.value)
)
cor_df$padj <- p.adjust(cor_df$p_value, method = "BH")
cor_df
cor_df$PC <- factor(cor_df$PC,
                    levels = paste0("PC", 1:10))
cor_df$significant <- cor_df$padj < 0.05
ggplot(cor_df, aes(x = PC, y = correlation, fill = significant)) +
  geom_col() +
  geom_text(aes(label = signif(padj, 2)), vjust = -0.5) +
  scale_fill_manual(values = c("FALSE" = "grey70",
                               "TRUE"  = "red")) +
  labs(
    title = "Spearman correlation between PCA axes and RIN (FDR adjusted)",
    y = "Spearman correlation",
    x = "Principal Components"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(plot=last_plot(), "spearman_correlation_RIN_PCA.png")

summary(lm(PC1 ~ RIN, data = scores))

#####################################################################################################
###########" same without outlier pour voir si il influence la PCA et les correlation #####################################

## suppression de l'outlier
outlier = c("Ins.1837", "Nac.1837")
coldata <- coldata %>%
  filter(!sample %in% outlier)
counts <- counts[, coldata$sample, drop = FALSE]
coverage <- coverage[coldata$sample, drop = FALSE,]

## normalization because PCA shows distances
vst_counts <- varianceStabilizingTransformation(as.matrix(counts))

# PCA classique
pca_res <- prcomp(t(vst_counts), scale. = TRUE)
# scores <- as.data.frame(pca_res$x)
# scores$sample <- rownames(scores)
# scores <- scores %>%
#   left_join(coldata, by = "sample")
# var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
# ggplot(scores, aes(x = PC1, y = PC2, color=reg)) +
#   geom_point(size = 2) +
#   labs(
#     x = paste0("PC1 (", round(var_explained[1]*100, 1), "%)"),
#     y = paste0("PC2 (", round(var_explained[2]*100, 1), "%)"),
#     title = "PCA on VST normalized counts, without ouliers"
#   ) +
#   theme_minimal()


# Scores des samples
scores <- as.data.frame(pca_res$x[, 1:10])  # PC1 à PC10
scores$sample <- rownames(scores)

# Ajouter RIN
scores <- scores %>%
  left_join(coldata[, c("sample", "RIN")], by = "sample")

cor_results <- lapply(scores[, 1:10], function(pc) {
  cor.test(pc, scores$RIN, method = "pearson")
})

cor_results

cor_df <- data.frame(
  PC = paste0("PC", 1:10),
  correlation = sapply(cor_results, function(x) x$estimate),
  p_value = sapply(cor_results, function(x) x$p.value)
)
cor_df$padj <- p.adjust(cor_df$p_value, method = "BH")
cor_df

cor_df$PC <- factor(cor_df$PC,
                    levels = paste0("PC", 1:10))
cor_df$significant <- cor_df$padj < 0.05

ggplot(cor_df, aes(x = PC, y = correlation, fill = significant)) +
  geom_col() +
  geom_text(aes(label = signif(padj, 2)), vjust = -0.5) +
  scale_fill_manual(values = c("FALSE" = "grey70",
                               "TRUE"  = "red")) +
  labs(
    title = "Pearson correlation between PCA axes and RIN (FDR adjusted)",
    y = "Pearson correlation",
    x = "Principal Components"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(plot=last_plot(), "pearson_correlation_RIN_PCA_without_outlier.png")


cor_results <- lapply(scores[, 1:10], function(pc) {
  cor.test(pc, scores$RIN, method = "spearman")
})

cor_results

cor_df <- data.frame(
  PC = paste0("PC", 1:10),
  correlation = sapply(cor_results, function(x) x$estimate),
  p_value = sapply(cor_results, function(x) x$p.value)
)
cor_df$padj <- p.adjust(cor_df$p_value, method = "BH")
cor_df
cor_df$PC <- factor(cor_df$PC,
                    levels = paste0("PC", 1:10))
cor_df$significant <- cor_df$padj < 0.05
ggplot(cor_df, aes(x = PC, y = correlation, fill = significant)) +
  geom_col() +
  geom_text(aes(label = signif(padj, 2)), vjust = -0.5) +
  scale_fill_manual(values = c("FALSE" = "grey70",
                               "TRUE"  = "red")) +
  labs(
    title = "Spearman correlation between PCA axes and RIN (FDR adjusted)",
    y = "Spearman correlation",
    x = "Principal Components"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(plot=last_plot(), "spearman_correlation_RIN_PCA_without_outlier.png")

###############################################################################################
### Etude du RIN pour region ###

coldata <- read_ods(coldata_path)

comparisons <- combn(unique(coldata$reg), 2, simplify = FALSE)
p <- ggplot(coldata, aes(x = reg, y = RIN, fill = reg)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  
  # # Test global
  # stat_compare_means(method = "kruskal.test",
  #                    label.y = max(coldata$RIN) + 0.8) +
  
  # Comparaisons 2 à 2
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test",
                     p.adjust.method = "BH",
                     label = "p.signif",
                     hide.ns = TRUE) +
  
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = "Pairwise comparison of RIN across regions",
    x = "Region",
    y = "RIN"
  )

p
ggsave(plot=last_plot(), "RIN_boxplots_region.png")



