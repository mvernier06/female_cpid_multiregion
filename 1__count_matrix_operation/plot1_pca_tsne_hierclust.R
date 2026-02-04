library(tidyverse)
library(FactoMineR)
library(factoextra)
library(M3C)
library(DESeq2)
library(corrplot)
library(readODS)
library(ggrepel)

rm(list=ls())

setwd("//wsl.localhost/Ubuntu/home/marinevernier/projets/cpid_multiregion/")

#### PATHS ####
raw_counts_path <- "female_cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
coldata_path <- "female_cpid_multiregion/data/count_data/coldata.ods"
plot.path <- "female_cpid_multiregion/graphs_results/1__count_matrix_operation/pca/"
# output.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/graphs_results/1__count_matrix_operation/"

coldata <- read_ods(coldata_path)

raw_counts <- read.csv(raw_counts_path) %>%
  column_to_rownames("MGI.symbol")

## normalization because PCA shows distances
vst_counts <- varianceStabilizingTransformation(as.matrix(raw_counts))


##### PCA ON RAW COUNTS ####
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample


pca(vst_counts, dotsize=3) +
 labs( title="PCA on VST normalized counts") # PCA classique
res <- M3C(vst_counts, removeplots = TRUE, iters=25,
          objective='PAC', fsize=8, lthick=1, dotsize=1.25) # M3C fait un clusreing consensus 
res

res$plots[[1]]
setwd(plot.path)
res$plots[[2]]
ggsave(plot=last_plot(), "PAC_raw_counts.png")
res$plots[[4]]
ggsave(plot=last_plot(), "RCSI_raw_counts.png")
res$plots[[3]]
ggsave(plot=last_plot(), "pval_raw_counts.png")
data <- res$realdataresults[[5]]$ordered_data
annon <- res$realdataresults[[5]]$ordered_annotation
ccmatrix <- res$realdataresults[[5]]$consensus_matrix
head(annon)


pca(data,labels=annon$consensuscluster,legendtextsize = 10,axistextsize = 10,dotsize=2) +
  labs( title="PCA on VST normalized counts ") # affichange de la pca avec le clustering consensus 
ggsave(plot=last_plot(), "PCA_raw_counts.png")

# On rejoute les colonnes qui nous intéressent pour colorier la PCA 
annon <- annon %>% rownames_to_column("sample")
annon <- annon %>%
  left_join(
    coldata %>%  select(sample, reg, group),
    by = "sample"
  )
# PCA sur les données ordonnées
pca <- prcomp(t(data), scale. = FALSE)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample") %>%
  left_join(annon, by = "sample")

ggplot(pca_df, aes(x = PC1, y = PC2, color = reg)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(
    title = "PCA colored by regions ",
    x = paste0("PC1 (", round(100*summary(pca)$importance[2,1],0), "%)"),
    y = paste0("PC2 (", round(100*summary(pca)$importance[2,2],0), "%)")
  )
ggsave(plot = last_plot(), "PCA_raw_counts_colored_by_region.PNG")

ggplot(pca_df, aes(x = PC1, y = PC2, color = reg)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = subset(pca_df, sample == "Ins.1837"),
    aes(label = sample),
    size = 4,
    fontface = "bold"
  ) +
  theme_classic() +
  labs(
    title = "PCA colored by regions",
    x = paste0("PC1 (", round(100*summary(pca)$importance[2,1],0), "%)"),
    y = paste0("PC2 (", round(100*summary(pca)$importance[2,2],0), "%)"))

ggsave("PCA_raw_counts_colored_by_region_labeled.PNG",
       width = 8, height = 6, dpi = 300)

###############################################################################################
#### HEAT MAP (fig1)####
#setwd("/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/graphs_results/1__count_matrix_operation/")
df <- scale(vst_counts)
df

# sort columns by region, tp and sham/cuff
new_order <- coldata %>% arrange(reg, timepoint, desc(group))
df2 <- df[, rownames(new_order)]

# Calculate the correlation matrix
cor_matrix <- cor(df2)
# Create a basic correlation heatmap using corrplot
png(file="heatmap_raw_counts_corrplot.png", width = 4000, height = 4000, units = "px", res = 100)
corrplot(cor_matrix, 
         method = "color", 
         col.lim = c(0,1), 
         col = COL2('RdBu', 200), 
         is.corr = FALSE, 
         tl.col = "grey10"
)
dev.off()

anno <- new_order
# palettes
reg_cols   <- setNames(RColorBrewer::brewer.pal(length(unique(anno$reg)), "Set2"),
                       unique(anno$reg))

tp_cols    <- setNames(RColorBrewer::brewer.pal(length(unique(anno$timepoint)), "Dark2"),
                       unique(anno$timepoint))

group_cols <- c("sham" = "#4DAF4A", "cuff" = "#E41A1C")

png(file="heatmap_regions.png", width = 4000, height = 4000, units = "px", res = 100)
label_colors <- reg_cols[anno$reg]
corrplot(
  cor_matrix,
  method = "color",
  col.lim = c(0,1),
  col = COL2('RdBu', 200),
  is.corr = FALSE,
  tl.col = label_colors,
  tl.cex = 0.6,        # taille texte
  tl.srt = 90          # rotation
)
legend("topright", legend = names(reg_cols), fill = reg_cols,
       title = "Region", cex = 1.2, bty = "n")
dev.off()

png(file="heatmap_time_points.png", width = 4000, height = 4000, units = "px", res = 100)
label_colors <- tp_cols[anno$timepoint]
corrplot(
  cor_matrix,
  method = "color",
  col.lim = c(0,1),
  col = COL2('RdBu', 200),
  is.corr = FALSE,
  tl.col = label_colors,
  tl.cex = 0.6,        # taille texte
  tl.srt = 90          # rotation
)

legend("topright", legend = names(tp_cols), fill = tp_cols,
       title = "Time points", cex = 1.2, bty = "n")
dev.off()

png(file="heatmap_group.png", width = 4000, height = 4000, units = "px", res = 100)
label_colors <- group_cols[anno$group]
corrplot(
  cor_matrix,
  method = "color",
  col.lim = c(0,1),
  col = COL2('RdBu', 200),
  is.corr = FALSE,
  tl.col = label_colors,
  tl.cex = 0.6,        # taille texte
  tl.srt = 90          # rotation
)

legend("bottomright", legend = names(group_cols), fill = group_cols,
       title = "Group", cex = 1.2, bty = "n")
dev.off()

#################################################################
#### T-SNE ON RAW COUNTS ####
vst_counts
setwd(output.path)
set.seed(123)
t <- tsne(vst_counts, labels=coldata$reg, legendtextsize = 15,axistextsize = 15,dotsize=3,
          colvec = "black") +
  labs( title="T-sne on VST normalized counts")
t
ggsave(plot=t, "tsne_raw_counts.png")
?tsne
t$data$sample <- coldata$sample
t_labeled <- t +
  geom_text_repel(
    data = subset(t$data, sample == "Ins.1837"),
    aes(label = sample),
    size = 4,
    fontface = "bold",
    color = "black"
  )

t_labeled
ggsave("tsne_raw_counts_Ins1837_labeled.png",
       plot = t_labeled,
       width = 8, height = 6, dpi = 300)

t_labeled_all <- t +
  geom_text_repel(
    data = t$data,
    aes(label = sample),
    size = 3,           # taille réduite pour 163 samples
    max.overlaps = Inf, # permet de forcer l'affichage de tous
    segment.size = 0.3  # petite ligne de liaison
  )

t_labeled_all

ggsave("tsne_raw_counts_all_samples_labeled.png",
       plot = t_labeled_all,
       width = 12, height = 8, dpi = 300)

t_labeled <- t +
  geom_text_repel(
    data = subset(t$data, sample == "Hb.1839"),
    aes(label = sample),
    size = 4,
    fontface = "bold",
    color = "black"
  )

t_labeled
ggsave("tsne_raw_counts_Hb.1839_labeled.png",
       plot = t_labeled,
       width = 8, height = 6, dpi = 300)


library(patchwork)

plots <- list()

for(i in 1:9){   
  set.seed(i)    # ensures each run is different but reproducible
  
  ti <- tsne(
    vst_counts,
    labels = coldata$reg,
    legendtextsize = 15,
    axistextsize = 15,
    dotsize = 3
  ) +
    labs(title = paste("T-SNE run", i, "on VST normalized counts"))
  ti$data$sample <- coldata$sample
  t_labeled_all <- ti +
    geom_text_repel(
      data = ti$data,
      aes(x = X1, y = X2, label = sample),  # coordonnées explicites
      size = 3,
      max.overlaps = Inf,
      segment.size = 0.3,
      box.padding = 0.4,     # espace autour du texte
      point.padding = 0.3,   # espace autour du point
      seed = 42              # positions reproductibles
    )
  
  plots[[i]] <- t_labeled_all
}

wrap_plots(plots, ncol = 3)
plots[[1]]

###################################################################################################################
#################################      PCA intra region      ######################################################

setwd("/home/marinevernier/projets/cpid_multiregion/female_cpid_multiregion/graphs_results/1__count_matrix_operation/pca/pca_intra_reg/")

sample_ins <- coldata %>%
  filter(reg == "Ins") %>%
  pull(sample)

counts_ins <- vst_counts[, sample_ins]

pca_ins <- prcomp(t(counts_ins), scale. = FALSE)

pca_df <- as.data.frame(pca_ins$x) %>%
  rownames_to_column("sample") %>%
  left_join(
    coldata %>% select(sample, timepoint, group),
    by = "sample"
  )

pca_df$timepoint <- as.factor(pca_df$timepoint)
ggplot(pca_df, aes(PC1, PC2, color = timepoint, label = sample)) +
  geom_point(size = 3) +
  # geom_text_repel(
  #   size = 3,
  #   max.overlaps = Inf,
  #   show.legend = FALSE
  # ) +
  geom_text_repel(
    data = subset(pca_df, sample == "Ins.1837"),
    aes(label = sample),
    size = 3
  )+
  theme_classic() +
  labs(
    title = "PCA intra-région Ins",
    x = paste0("PC1 (", round(100 * summary(pca_ins)$importance[2,1], 0), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_ins)$importance[2,2], 0), "%)")
  )
ggsave(plot=last_plot(), "PCA_Ins_labeled_time_point.png")

##

sample_acc <- coldata %>%
  filter(reg == "ACC") %>%
  pull(sample)

counts_acc <- vst_counts[, sample_acc]

pca_acc <- prcomp(t(counts_acc), scale. = FALSE)

pca_df <- as.data.frame(pca_acc$x) %>%
  rownames_to_column("sample") %>%
  left_join(
    coldata %>% select(sample, timepoint, group),
    by = "sample"
  )

pca_df$timepoint <- as.factor(pca_df$timepoint)
ggplot(pca_df, aes(PC1, PC2, color = timepoint, label = sample)) +
  geom_point(size = 3) +
  # geom_text_repel(
  #   size = 3,
  #   max.overlaps = Inf,
  #   show.legend = FALSE
  # ) +
  # geom_text_repel(
  #   data = subset(pca_df, sample == "Ins.1837"),
  #   aes(label = sample),
  #   size = 3
  # )+
  theme_classic() +
  labs(
    title = "PCA intra-région ACC",
    x = paste0("PC1 (", round(100 * summary(pca_acc)$importance[2,1], 0), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_acc)$importance[2,2], 0), "%)")
  )
ggsave(plot=last_plot(), "PCA_Acc_time_point.png")

###

sample_Hb <- coldata %>%
  filter(reg == "Hb") %>%
  pull(sample)

counts_Hb <- vst_counts[, sample_Hb]

pca_Hb <- prcomp(t(counts_Hb), scale. = FALSE)

pca_df <- as.data.frame(pca_Hb$x) %>%
  rownames_to_column("sample") %>%
  left_join(
    coldata %>% select(sample, timepoint, group),
    by = "sample"
  )

pca_df$timepoint <- as.factor(pca_df$timepoint)
ggplot(pca_df, aes(PC1, PC2, color = group, label = sample)) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  # geom_text_repel(
  #   data = subset(pca_df, sample == "Ins.1837"),
  #   aes(label = sample),
  #   size = 3
  # )+
  theme_classic() +
  labs(
    title = "PCA intra-région Hb",
    x = paste0("PC1 (", round(100 * summary(pca_Hb)$importance[2,1], 0), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_Hb)$importance[2,2], 0), "%)")
  )
ggsave(plot=last_plot(), "PCA_Hb_labeled_group.png")

##

sample_nac <- coldata %>%
  filter(reg == "Nac") %>%
  pull(sample)

counts_nac <- vst_counts[, sample_nac]

pca_nac <- prcomp(t(counts_nac), scale. = FALSE)

pca_df <- as.data.frame(pca_nac$x) %>%
  rownames_to_column("sample") %>%
  left_join(
    coldata %>% select(sample, timepoint, group),
    by = "sample"
  )

pca_df$timepoint <- as.factor(pca_df$timepoint)
ggplot(pca_df, aes(PC1, PC2, color = group, label = sample)) +
  geom_point(size = 3) +
  # geom_text_repel(
  #   size = 3,
  #   max.overlaps = Inf,
  #   show.legend = FALSE
  # ) +
  geom_text_repel(
    data = subset(pca_df, sample == "Nac.1837"),
    aes(label = sample),
    size = 3
  )+
  theme_classic() +
  labs(
    title = "PCA intra-région Nac",
    x = paste0("PC1 (", round(100 * summary(pca_nac)$importance[2,1], 0), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_nac)$importance[2,2], 0), "%)")
  )
ggsave(plot=last_plot(), "PCA_nac_labeled_group.png")

