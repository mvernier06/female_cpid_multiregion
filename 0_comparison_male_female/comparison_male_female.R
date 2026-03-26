library(tidyverse)
library(FactoMineR)
library(factoextra)
library(M3C)
library(DESeq2)
library(corrplot)
library(readODS)


rm(list=ls())

counts_female <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
counts_male <- "/home/marinevernier/Documents/projets/cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
coldata_female_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/counts_m39_M32/coldata_without_outliers.ods"
coldata_male_path <- "/home/marinevernier/Documents/projets/cpid_multiregion/data/2__differential_expression_analysis/coldata.csv"
output_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/0_comparison_male_female/"
setwd(output_path)

df_female <- read_csv(counts_female)
coldata_female <- read_ods(coldata_female_path)
df_male <- read_csv(counts_male)
coldata_male <- read_csv(coldata_male_path)

coldata_female <- coldata_female %>%
  mutate(sex = "female")

coldata_male <- coldata_male %>%
  mutate(sex = "male")

df_female <- df_female %>% column_to_rownames("MGI.symbol")
df_male <- df_male %>% column_to_rownames("MGI.symbol")

common_genes <- intersect(rownames(df_female), rownames(df_male))


counts_all <- cbind(
  df_female[common_genes, ],
  df_male[common_genes, ]
)


coldata_all <- bind_rows(coldata_female, coldata_male) %>%
  mutate(reg = if_else(reg == "Nac", "NAc", reg)) %>%
  column_to_rownames("sample")

# Vérification
stopifnot(all(colnames(counts_all) == rownames(coldata_all)))


dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData = coldata_all,
  design = ~ sex + reg + group
)

vsd <- vst(dds, blind = TRUE)


### PCA 1 : toutes les régions ###

pca_all <- prcomp(t(assay(vsd)), scale. = FALSE)
pca_all_df <- as.data.frame(pca_all$x) %>%
  rownames_to_column("sample") %>%
  left_join(coldata_all %>% rownames_to_column("sample"), by = "sample")

ggplot(pca_all_df, aes(PC1, PC2, color = reg, shape = sex)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    title = "PCA sur toutes les régions",
    x = paste0("PC1 (", round(100*summary(pca_all)$importance[2,1], 1), "%)"),
    y = paste0("PC2 (", round(100*summary(pca_all)$importance[2,2], 1), "%)")
  )

ggsave(plot = last_plot(), filename = "PCA_all_reg.png", width = 6, height = 5)


### PCA 2 : régions présentes dans les deux sexes ###


regions_common <- coldata_all %>%
  group_by(reg) %>%
  summarise(sexes = n_distinct(sex)) %>%
  filter(sexes == 2) %>%
  pull(reg)


coldata_common <- coldata_all %>% filter(reg %in% regions_common)
vsd_common <- vsd[, rownames(coldata_common)]

# PCA sur régions communes
pca_common <- prcomp(t(assay(vsd_common)), scale. = FALSE)
pca_common_df <- as.data.frame(pca_common$x) %>%
  rownames_to_column("sample") %>%
  left_join(coldata_common %>% rownames_to_column("sample"), by = "sample")

ggplot(pca_common_df, aes(PC1, PC2, color = reg, shape = sex)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    title = "PCA sur régions communes aux deux sexes",
    x = paste0("PC1 (", round(100*summary(pca_common)$importance[2,1], 1), "%)"),
    y = paste0("PC2 (", round(100*summary(pca_common)$importance[2,2], 1), "%)")
  )

ggsave(plot = last_plot(), filename = "PCA_common_reg.png", width = 6, height = 5)




ggplot(pca_common_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  facet_wrap(~ sex) +
  theme_classic()
ggsave(plot = last_plot(), "PCA_common_reg_group_male_vs_female.png")

ggplot(pca_common_df, aes(PC1, PC2, color = sex)) +
  geom_point(size = 3) +
  theme_classic()
ggsave(plot = last_plot(), "PCA_common_reg_sex.png")


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

df_unfiltered_female <- read_csv("/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts.csv")
df_unfiltered_male <- read_csv("/home/marinevernier/Documents/projets/cpid_multiregion/data/2__differential_expression_analysis/annotated_counts.csv")

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
###       Comparaison des n() et DEG des males et femelles ###


table(coldata_female$reg, coldata_female$timepoint, coldata_female$group)
table(coldata_male$reg, coldata_male$timepoint, coldata_male$group)

table(coldata_all$sex)
table(coldata_all$sex,
      coldata_all$reg,
      coldata_all$timepoint,
      coldata_all$group)

df_plot <- coldata_all %>%
  as.data.frame() %>%
  dplyr::count(sex, reg, timepoint, group)


ggplot(df_plot,
       aes(x = interaction( timepoint, reg),
           y = n,
           fill = sex)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8)) +
  facet_wrap(~ group, nrow = 1) +   
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 12, face = "bold")) +
  labs(x = "Region_Timepoint",
       y = "N samples",
       fill = "Sex")

common_reg <- c("Ins", "Hb", "NAc")
df_common <- df_plot %>%
  filter(reg %in% common_reg)
df_cuff <- df_common %>%
  filter(group == "cuff")

ggplot(df_cuff,
       aes(x = timepoint,
           y = n,
           fill = sex)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8)) +
  facet_wrap(~ reg) +
  scale_y_continuous(limits = c(0, 8),
                     breaks = 0:8) +
  theme_minimal() +
  labs(title = "Cuff (common regions only)",
       x = "Timepoint",
       y = "N samples",
       fill = "Sex")
ggsave(plot=last_plot(), "number_of_sample_for_cuff_male_and_female.png")

df_sham <- df_common %>%
  filter(group == "sham")

ggplot(df_sham,
       aes(x = timepoint,
           y = n,
           fill = sex)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8)) +
  facet_wrap(~ reg) +
  scale_y_continuous(limits = c(0, 8),
                     breaks = 0:8) +
  theme_minimal() +
  labs(title = "Sham (common regions only)",
       x = "Timepoint",
       y = "N samples",
       fill = "Sex")
ggsave(plot=last_plot(), "number_of_sample_for_sham_male_and_female.png")

deglist_male <- "/home/marinevernier/Documents/projets/cpid_multiregion/data/2__differential_expression_analysis/deglist.Rdata"
deglist_female <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/deglist.Rdata"


# -------- LOAD MALE --------
male_env <- new.env()
load(deglist_male, envir = male_env)

male_objs <- ls(male_env, pattern = "^deg_")

for(obj in male_objs){
  assign(paste0(obj, "_male"), male_env[[obj]])
}


# -------- LOAD FEMALE --------
female_env <- new.env()
load(deglist_female, envir = female_env)

female_objs <- ls(female_env, pattern = "^deg_")

for(obj in female_objs){
  assign(paste0(obj, "_female"), female_env[[obj]])
}

deg_NAc_tp1_female <- deg_Nac_tp1_female
deg_NAc_tp2_female <- deg_Nac_tp2_female
deg_NAc_tp3_female <- deg_Nac_tp3_female

all_deg_objs <- ls(pattern = "^deg_")

deg_counts <- map_df(all_deg_objs, function(obj){
  
  df <- get(obj)
  
  tibble(
    object = obj,
    n_deg = nrow(df)
  )
})

deg_counts <- deg_counts %>%
  separate(object,
           into = c("deg", "reg", "timepoint", "sex"),
           sep = "_") %>%
  select(reg, timepoint, sex, n_deg)


deg_counts <- deg_counts %>%
  mutate(reg = case_when(
    tolower(reg) == "nac" ~ "NAc",  # corrige Nac/NAC/nac
    TRUE ~ reg                     # les autres restent identiques
  ))

deg_counts <- deg_counts %>% filter(reg %in% common_reg)

ggplot(deg_counts,
       aes(x = timepoint,
           y = n_deg,
           fill = sex)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8)) +
  facet_wrap(~ reg) +
  theme_minimal() +
  labs(title = "Number of DEG per region and timepoint",
       x = "Timepoint",
       y = "Number of DEG",
       fill = "Sex")
ggsave(plot=last_plot(), "number_of_DEG_male_vs_female.png")

ggplot(deg_counts,
       aes(x = timepoint,
           y = n_deg,
           color = reg,
           group = reg)) +
  geom_point(aes(pch=reg)) +
  geom_line(aes(group=reg)) +
  facet_wrap(~ sex) +
  labs(title="DEG kinetics per TP within regions",
       y = "Number of DEG",
       color = "Region", shape = "Region") +
  theme_bw() +
  theme(axis.title.x=element_blank())

ggsave("kinetics_DEG_by_sex.png", width = 8, height = 5)
##################
## l2fc kinetic ##
##################
deg_fc <- map_df(all_deg_objs, function(obj){
  
  df <- get(obj)
  
  df %>%
    mutate(object = obj) %>%
    select(object, log2fc)
})

deg_fc <- deg_fc %>%
  separate(object,
           into = c("deg", "reg", "timepoint", "sex"),
           sep = "_") %>%
  select(reg, timepoint, sex, log2fc)
deg_fc <- deg_fc %>%
  mutate(reg = case_when(
    tolower(reg) == "nac" ~ "NAc",
    TRUE ~ reg
  )) %>%
  filter(reg %in% common_reg)
deg_fc <- deg_fc %>%
  mutate(fcabs = abs(log2fc))

ggplot(deg_fc,
       aes(x = reg,
           y = fcabs,
           color = timepoint,
           fill = timepoint)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ sex) +
  labs(y = "|Log2FC|",
       title = "Absolute Fold Change values per TP for each region",
       fill = "Timepoint",
       color = "Timepoint") +
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_color_manual(values = c("#8594B6", "#8A4F21", "#465553")) +
  scale_fill_manual(values = c("#F3F6FB", "#EDD9BA", "#D9E8E5")) +
  theme_bw() +
  theme(axis.title.x = element_blank())
ggsave(plot=last_plot(), "l2fc_kinetics_by_sex.png")

#######################
###    upset pot    ###
#######################
setwd("upset_plot")

timepoints <- c("tp1", "tp2", "tp3")

for(tp in timepoints){
  
  cat("Processing", tp, "\n")
  
  listInput <- list()
  
  for(reg in common_reg){
    
    # noms objets
    obj_male   <- paste0("deg_", reg, "_", tp, "_male")
    obj_female <- paste0("deg_", reg, "_", tp, "_female")
    
    # récupérer data
    deg_male   <- get(obj_male)
    deg_female <- get(obj_female)
    
    # extraire gènes
    listInput[[paste0(reg, "_male")]]   <- deg_male$label
    listInput[[paste0(reg, "_female")]] <- deg_female$label
  }
  
  # ---------- DEBUG ----------
  print(sapply(listInput, length))
  
  # ---------- PLOT ----------
  png(
    file = paste0("deg_upset_", tp, "_male_vs_female.png"),
    width = 1200, height = 800, res = 120
  )

  print(
    upset(
      fromList(listInput),
      nsets = length(listInput),   
      keep.order = TRUE,
      mainbar.y.label = paste0("DEGs intersections (Male vs Female, ", tp, ")"),
      sets.x.label = "Number of DEGs",
      order.by = "freq"
    )
  )

  dev.off()
}

for(reg in common_reg){
  
  listInput <- list()
  
  for(tp in timepoints){
    
    deg_male   <- get(paste0("deg_", reg, "_", tp, "_male"))
    deg_female <- get(paste0("deg_", reg, "_", tp, "_female"))
    
    listInput[[paste0(tp, "_male")]]   <- deg_male$label
    listInput[[paste0(tp, "_female")]] <- deg_female$label
  }
  
  png(
    file = paste0("deg_upset_", reg, "_timepoints_sex.png"),
    width = 1400, height = 900, res = 120
  )
  
  print(
    upset(
      fromList(listInput),
      nsets = length(listInput),
      keep.order = TRUE,
      order.by = "freq"
    )
  )
  
  dev.off()
}
