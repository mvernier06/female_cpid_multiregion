library(tidyverse)
library(readODS)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(M3C)
library(ggpubr)

rm(list=ls())

project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion"
setwd(project_path)

#### PATHS ####
raw_counts_path <- "data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv" # attention, il n'y a pas les outliers 
coldata_path <- "data/counts_m39_M32/coldata.ods"
plot.path <- "graphs_results/1__count_matrix_operation/quality_check/"
coverage.path <- "data/counts_m39_M32/coverage_per_sample.tsv"


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
pca_res <- prcomp(t(vst_counts))
 
pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("sample") %>%
  left_join(
    coldata %>% select(sample, timepoint, group, RIN, reg),
    by = "sample"
  )

pca_df$timepoint <- as.factor(pca_df$timepoint)
ggplot(pca_df, aes(PC1, PC2, color = RIN, label = sample)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    title = "PCA ",
    x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2,1], 0), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2,2], 0), "%)")
  )
ggsave(plot=last_plot(), "PCA_colored_by_RIN.png")

# Scores des samples
scores <- as.data.frame(pca_res$x[, 1:10])  # PC1 à PC10
scores$sample <- rownames(scores)

# Ajouter RIN
scores <- scores %>%
  left_join(coldata[, c("sample", "RIN", "reg", "timepoint", "group")], by = "sample")

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

ggsave(plot=last_plot(), "spearman_correlation_RIN_PCA_without.png")


ggplot(scores, aes(x = RIN, y = PC1)) + # , color = reg
  geom_point(size = 3, alpha = 0.8) +
  # geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    title = "PC1 projection according to RIN score",
    x = "RIN",
    y = "PC1 score"
  ) +
  theme_minimal()

ggsave("PC1_vs_RIN_scatter.png", width = 6, height = 5)

ggplot(scores, aes(RIN, PC1, color=reg)) +
  geom_point() +
  labs(
    title = "PC1 projection according to RIN score",
    x = "RIN",
    y = "PC1 score"
  ) +
  geom_smooth(method="lm")
ggsave(plot=last_plot(), "PC1_vs_RIN_linear_regression.png")

###############################################################################################
### Etude du RIN pour region ###

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


ggplot(coldata, aes(x = group, y = RIN)) +
  geom_boxplot() +
  geom_jitter(width = 0.1)
ggsave(plot=last_plot(), "RIN_boxplots_group.png")


comparisons <- combn(unique(coldata$reg), 2, simplify = FALSE)

p <- ggplot(coldata, aes(x = reg, y = RIN, fill = reg)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  
  facet_grid( group ~ timepoint) + # group ~ timepoint
  
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = "Pairwise comparison of RIN across regions",
    x = "Region",
    y = "RIN"
  )
p
ggsave(plot = last_plot(), "pairwise_comparison_RIN_across_reg_timepoint_group.png")
coldata$cond <- interaction(coldata$reg, coldata$group, coldata$timepoint)

comparisons <- combn(unique(coldata$cond), 2, simplify = FALSE)

ggplot(coldata, aes(x = cond, y = RIN, fill = reg)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2) +
  
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  
  theme_classic() +
  labs(x = "Region / Group / Timepoint")

####################################################################################################################
### Test de GLM-PCA ###

outliers <- c("Ins.1837", "Nac.1837", "Hb.1839", "Hb.2049")
coldata <- coldata %>%
  filter(!sample %in% outliers)
print(coldata$sample==colnames(counts))

library("glmpca")
gpca <- glmpca(counts, L=2)
gpca.dat <- gpca$factors
gpca.dat$group <- coldata$group
gpca.dat$reg <- coldata$reg
gpca.dat$timepoint <- coldata$timepoint
gpca.dat$RIN <- coldata$RIN

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = RIN, shape = group)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
