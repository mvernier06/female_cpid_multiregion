library(tidyverse)
library(FactoMineR)
library(factoextra)
library(M3C)
library(DESeq2)
library(corrplot)
library(readODS)

rm(list=ls())

counts_female <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
counts_male <- "/home/marinevernier/Documents/cpid_multiregion/cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
coldata_female_path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/count_data/coldata.ods"
coldata_male_path <- "/home/marinevernier/Documents/cpid_multiregion/cpid_multiregion/data/2__differential_expression_analysis/coldata.csv"
output_path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/graphs_results/0_comparison_male_female/"
setwd(output_path)

df_female <- read_csv(counts_female)
coldata_female <- read_ods(coldata_female_path)
df_male <- read_csv(counts_male)
coldata_male <- read_csv(coldata_male_path)

coldata_female <- coldata_female %>%
  mutate(sex = "female")

coldata_male <- coldata_male %>%
  mutate(sex = "male")

df_female <- df_female %>%
  column_to_rownames("MGI.symbol")  

df_male <- df_male %>%
  column_to_rownames("MGI.symbol")

common_genes <- intersect(rownames(df_female), rownames(df_male))

counts_all <- cbind(
  df_female[common_genes, ],
  df_male[common_genes, ]
)

coldata_all <- bind_rows(coldata_female, coldata_male) %>%
  mutate(
    reg = if_else(reg == "Nac", "NAc", reg)
  ) %>%
  column_to_rownames("sample")
unique(coldata_all$reg)

all(colnames(counts_all) == rownames(coldata_all))

dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData = coldata_all,
  design = ~ sex + reg + group
)

vsd <- vst(dds, blind = TRUE)

pca <- prcomp(t(assay(vsd)), scale. = FALSE)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample") %>%
  left_join(
    coldata_all %>% rownames_to_column("sample"),
    by = "sample"
  )

ggplot(pca_df, aes(PC1, PC2, color = reg, shape = sex)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    title = "PCA colored by regions, points by sex",
    x = paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "%)"),
    y = paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "%)") # ces deux lignes servent à afficher le pourcentage de variance expliquée par les PC 
)
ggsave(plot = last_plot(), "PCA_region_sex_all.png")
summary(pca)
ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  facet_wrap(~ sex) +
  theme_classic()
ggsave(plot = last_plot(), "PCA_group_male_vs_female.png")

ggplot(pca_df, aes(PC1, PC2, color = sex)) +
  geom_point(size = 3) +
  theme_classic()
ggsave(plot = last_plot(), "PCA_sex.png")

# En affichant seulement les regions communes : 
regions_communes <- pca_df %>%
  distinct(reg, sex) %>%
  dplyr::count(reg) %>%
  filter(n == 2) %>%
  pull(reg)
pca_df_filt <- pca_df %>%
  filter(reg %in% regions_communes)
ggplot(pca_df_filt, aes(PC1, PC2, color = reg, shape = sex)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    title = "PCA displaying only commune regions",
    x = paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "%)"),
    y = paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "%)")
  )
ggsave(plot = last_plot(), "PCA_commune_regions.png")


##################################################################################################################################
#                                 Comparaison de gènes spécifique des chromosomes X Y                                           #

"Xist" %in% rownames(df_female)
"Xist" %in% rownames(df_male)

dds_female <- DESeqDataSetFromMatrix(
  countData = df_female,
  colData = coldata_female %>% column_to_rownames("sample"),
  design = ~ 1
)

vsd_female <- vst(dds_female, blind = TRUE)

dds_male <- DESeqDataSetFromMatrix(
  countData = df_male,
  colData = coldata_male %>% column_to_rownames("sample"),
  design = ~ 1
)

vsd_male <- vst(dds_male, blind = TRUE)

xist_df <- bind_rows(
  tibble(
    sample = colnames(assay(vsd_female)),
    expression = assay(vsd_female)["Xist", ],
    sex = "female"
  ),
  tibble(
    sample = colnames(assay(vsd_male)),
    expression = if ("Xist" %in% rownames(assay(vsd_male)))
      assay(vsd_male)["Xist", ] else NA,
    sex = "male"
  )
)

ggplot(xist_df, aes(x = sex, y = expression, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "Expression de Xist",
    y = "Expression (VST)",
    x = ""
  )
ggsave(plot = last_plot(), "Xlist_boxplot_male_female.png")

"Ddx3y" %in% rownames(df_female)
"Ddx3y" %in% rownames(df_male)

Ddx3y <- bind_rows(
  tibble(
    sample = colnames(assay(vsd_male)),
    expression = assay(vsd_male)["Ddx3y", ],
    sex = "male"
  ),
  tibble(
    sample = colnames(assay(vsd_female)),
    expression = if ("Ddx3y" %in% rownames(assay(vsd_female)))
      assay(vsd_female)["Ddx3y", ] else NA,
    sex = "female"
  )
)
ggplot(Ddx3y, aes(x = sex, y = expression, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "Expression de Ddx3y",
    y = "Expression (VST)",
    x = ""
  )
ggsave(plot = last_plot(), "Ddx3y_boxplot_male_female.png")



Uty <- bind_rows(
  tibble(
    sample = colnames(assay(vsd_male)),
    expression = assay(vsd_male)["Uty", ],
    sex = "male"
  ),
  tibble(
    sample = colnames(assay(vsd_female)),
    expression = if ("Uty" %in% rownames(assay(vsd_female)))
      assay(vsd_female)["Uty", ] else NA,
    sex = "female"
  )
)
ggplot(Ddx3y, aes(x = sex, y = expression, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "Expression de Uty",
    y = "Expression (VST)",
    x = ""
  )
ggsave(plot = last_plot(), "Uty_boxplot_male_female.png")

######### Même chose mais en prenant les counts avant filtrage #############

df_unfiltered_female <- read_csv("/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts.csv")
df_unfiltered_male <- read_csv("/home/marinevernier/Documents/cpid_multiregion/cpid_multiregion/data/2__differential_expression_analysis/annotated_counts.csv")

coldata_female <- coldata_female %>%
  mutate(sex = "female")
coldata_female <- coldata_female %>% column_to_rownames("sample")

coldata_male <- coldata_male %>%
  mutate(sex = "male")
coldata_male <- coldata_male %>% column_to_rownames("sample")

df_unfiltered_female <- df_unfiltered_female %>%
  column_to_rownames("MGI.symbol")  
df_unfiltered_female <- df_unfiltered_female[,-1]

df_unfiltered_male <- df_unfiltered_male %>%
  column_to_rownames("MGI.symbol")
df_unfiltered_male <- df_unfiltered_male[,-1]

# Pour les femelles
df_samples_female <- df_unfiltered_female %>%
  select(all_of(rownames(coldata_female)))

# Pour les mâles
df_samples_male <- df_unfiltered_male %>%
  select(all_of(rownames(coldata_male)))


"Xist" %in% rownames(df_unfiltered_female)
"Xist" %in% rownames(df_unfiltered_male)


dds_female <- DESeqDataSetFromMatrix(
  countData = df_samples_female,
  colData = coldata_female,
  design = ~ 1
)

vsd_female <- vst(dds_female, blind = TRUE)

dds_male <- DESeqDataSetFromMatrix(
  countData = df_samples_male,
  colData = coldata_male,
  design = ~ 1
)

vsd_male <- vst(dds_male, blind = TRUE)

xist_df <- bind_rows(
  tibble(
    sample = colnames(assay(vsd_female)),
    expression = assay(vsd_female)["Xist", ],
    sex = "female"
  ),
  tibble(
    sample = colnames(assay(vsd_male)),
    expression = if ("Xist" %in% rownames(assay(vsd_male)))
      assay(vsd_male)["Xist", ] else NA,
    sex = "male"
  )
)

ggplot(xist_df, aes(x = sex, y = expression, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "Expression de Xist",
    y = "Expression (VST)",
    x = ""
  )
ggsave(plot = last_plot(), "Xlist_boxplot_male_female_unfiltered.png")

"Ddx3y" %in% rownames(df_samples_female)
"Ddx3y" %in% rownames(df_samples_male)

Ddx3y <- bind_rows(
  tibble(
    sample = colnames(assay(vsd_male)),
    expression = assay(vsd_male)["Ddx3y", ],
    sex = "male"
  ),
  tibble(
    sample = colnames(assay(vsd_female)),
    expression = if ("Ddx3y" %in% rownames(assay(vsd_female)))
      assay(vsd_female)["Ddx3y", ] else NA,
    sex = "female"
  )
)
ggplot(Ddx3y, aes(x = sex, y = expression, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "Expression de Ddx3y",
    y = "Expression (VST)",
    x = ""
  )
ggsave(plot = last_plot(), "Ddx3y_boxplot_male_female_unfiltered.png")

Uty <- bind_rows(
  tibble(
    sample = colnames(assay(vsd_male)),
    expression = assay(vsd_male)["Uty", ],
    sex = "male"
  ),
  tibble(
    sample = colnames(assay(vsd_female)),
    expression = if ("Uty" %in% rownames(assay(vsd_female)))
      assay(vsd_female)["Uty", ] else NA,
    sex = "female"
  )
)
ggplot(Ddx3y, aes(x = sex, y = expression, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "Expression de Uty",
    y = "Expression (VST)",
    x = ""
  )
ggsave(plot = last_plot(), "Uty_boxplot_male_female_unfiltered.png")


#####################################################################################################################################

