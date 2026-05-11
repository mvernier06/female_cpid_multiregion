rm(list = ls()) # rm R working space

library(MEGENA)
library(openxlsx)

project_path <- "/home/marinevernier/Documents/projets/"
setwd(project_path)

reg <- "Ins"
# PATHS ####
megena_results.path <- paste0("female_cpid_multiregion/data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/MEGENA.Results_", reg, ".Rdata")
modtable.path <- paste0("female_cpid_multiregion/data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/modtable_", reg, ".Rdata")
genes_chrX.path <- "female_cpid_multiregion/data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/genes_chrX_MGI_symbols.txt"
raw_counts_female.path <- "female_cpid_multiregion/data/counts_m39_M32/cpid_multireg_counts.txt"
raw_counts_male.path <- "cpid_multiregion/data/count_data/outliers_removed/un_normalised_counts.csv"
annot_table.path <- "female_cpid_multiregion/data/counts_m39_M32/annotation_final.csv"

output.path <- "female_cpid_multiregion/data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/infos_chrXY/"




### Formating female data ###
raw_counts_female <- read.table(raw_counts_female.path, header = TRUE,
                                sep = "\t",
                                comment.char = "#",
                                stringsAsFactors = FALSE,
                                check.names = FALSE)

#recuperer les gènes des chromosomes X et Y pour plus tard 
# Gènes chrX
genes_chrX <- raw_counts_female %>%
  filter(Chr == "chrX") %>%
  mutate(Geneid = str_replace(Geneid, "\\..*", "")) %>%
  select(Geneid) %>%
  distinct()

# Gènes chrY
genes_chrY <- raw_counts_female %>%
  filter(Chr == "chrY") %>%
  mutate(Geneid = str_replace(Geneid, "\\..*", "")) %>%
  select(Geneid) %>%
  distinct()

# Vérification
cat("Nombre de gènes chrX :", nrow(genes_chrX), "\n")
cat("Nombre de gènes chrY :", nrow(genes_chrY), "\n")



counts_female <- raw_counts_female[, !(colnames(raw_counts_female) %in% c(
  "Chr", "Start", "End", "Strand", "Length"
))]

# Formating colnames
colnames(counts_female) <- sub(
  "_R2\\.dedup\\.bam$",
  "",
  basename(colnames(counts_female))
)
colnames(counts_female) <- sub(
  "^([A-Za-z]+)([0-9]+)$",
  "\\1.\\2",
  colnames(counts_female)
)
counts_female$Geneid <-  str_replace(counts_female$Geneid, "\\..*", "") 

# annotate genes
ens2symbol <- read.csv(annot_table.path, sep=";")
ens2symbol <- ens2symbol %>% dplyr::select(!X)
counts_female <- inner_join(counts_female, ens2symbol, by=c("Geneid"="Gene.stable.ID"))
counts_female <- counts_female %>% relocate(MGI.symbol, .before=Geneid) 
counts_female <- counts_female[!duplicated(counts_female$MGI.symbol),] # remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)


genes_chrX <- genes_chrX %>%
  inner_join(ens2symbol, by = c("Geneid" = "Gene.stable.ID")) %>%
  dplyr::select(MGI.symbol) %>%
  distinct()
print(paste0("Number of genes on chrX after annotation: ", nrow(genes_chrX)))
genes_chrY <- genes_chrY %>%
  inner_join(ens2symbol, by = c("Geneid" = "Gene.stable.ID")) %>%
  dplyr::select(MGI.symbol) %>%
  distinct()
print(paste0("Number of genes on chrY after annotation: ", nrow(genes_chrY)))

# Keeping only Insula samples and excluding outlier
outlier <- "Ins.1837"
counts_female <- counts_female %>%
  distinct(MGI.symbol, .keep_all = TRUE) %>%
  column_to_rownames("MGI.symbol") %>%
  select(contains(paste0(reg, "."))) %>%
  select(!contains(outlier))
print(dim(counts_female))


### Formating male data ###
raw_counts_male <- read.csv(raw_counts_male.path, sep = ";")
raw_counts_male$X <- str_replace(raw_counts_male$X, "\\..*", "")
raw_counts_male <- inner_join(raw_counts_male, ens2symbol, by = c("X" = "Gene.stable.ID"))  
raw_counts_male <- raw_counts_male[!duplicated(raw_counts_male$MGI.symbol),] # remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)

print(dim(raw_counts_male))
filtered_counts_male <- raw_counts_male %>% filter(MGI.symbol %in% rownames(counts_female))
print(dim(filtered_counts_male))

counts_male <- filtered_counts_male %>%
  distinct(MGI.symbol, .keep_all = TRUE) %>%
  column_to_rownames("MGI.symbol") %>%
  select(contains(paste0(reg, "."))) 

print(dim(counts_male))

## tester le nombre de gènes Y qui sont à 0 
counts_genesY_male <- counts_male %>% rownames_to_column("MGI.symbol") %>% filter(MGI.symbol %in% genes_chrY$MGI.symbol)
counts_long <- counts_genesY_male %>%
  pivot_longer(
    cols = -MGI.symbol,
    names_to = "sample",
    values_to = "count"
  )
ggplot(counts_long, aes(x = count)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(
    title = "Distribution des counts - gènes chrY (male)",
    x = "Counts",
    y = "Fréquence")
mat <- as.matrix(counts_genesY_male[ , -1])

df_counts_summary <- data.frame(
  MGI.symbol = counts_genesY_male$MGI.symbol,
  n_zero = rowSums(mat == 0),
  n_non_zero = rowSums(mat != 0)
)



#keep genes with at least 1 cpm in >= 20% of samples in a single project ###
cpm_female <- apply(counts_female, 2, function(x) (x/sum(x))*1000000)
counts_female_filtered <- counts_female[rowSums(cpm_female < 1) < (dim(cpm_female)[2]/5),]
sum(genes_chrX$MGI.symbol %in% rownames(counts_female_filtered))

# CPM chrX
cpm_X <- cpm_female[rownames(cpm_female) %in% genes_chrX$MGI.symbol, ]
n_samples_X <- rowSums(cpm_X < 1)
df_X <- data.frame(
  gene = rownames(cpm_X),
  n_samples_cpm_inf_1 = n_samples_X,
  prop_samples_cpm_inf_1 = n_samples_X/ncol(cpm_X), 
  pass_filter = (n_samples_X < 8), 
  stringsAsFactors = FALSE
)
sum(df_X$n_samples_cpm_inf_1 < 8)
write.xlsx(df_X, file = file.path(output.path, "/infos_chrXY/genes_chrX_in_female.xlsx"))

# CPM chrY
cpm_Y <- cpm_female[rownames(cpm_female) %in% genes_chrY$MGI.symbol, ]
n_samples_Y <- rowSums(cpm_Y < 1)
df_Y <- data.frame(
  gene = rownames(cpm_Y),
  n_samples_cpm_inf_1 = n_samples_Y,
  prop_samples_cpm_inf_1 = n_samples_Y/ncol(cpm_Y), 
  stringsAsFactors = FALSE
)
sum(df_Y$n_samples_cpm_inf_1 < 8)

cpm_male <- apply(counts_male, 2, function(x)(x/sum(x))*1000000)
counts_male_filtered <- counts_male[rowSums(cpm_male < 1) < (dim(cpm_male)[2]/5),]

# CPM chrX male
cpm_X <- cpm_male[rownames(cpm_male) %in% genes_chrX$MGI.symbol, ]
n_samples_X <- rowSums(cpm_X < 1)
df_X <- data.frame(
  gene = rownames(cpm_X),
  n_samples_cpm_inf_1 = n_samples_X,
  prop_samples_cpm_inf_1 = n_samples_X/ncol(cpm_X), 
  pass_filter = (n_samples_X < 8), 
  stringsAsFactors = FALSE
)
sum(df_X$n_samples_cpm_inf_1 < 8)
write.xlsx(df_X, file = file.path(output.path, "/infos_chrXY/genes_chrX_in_male.xlsx"))

# CPM chrY male
cpm_Y <- cpm_male[rownames(cpm_male) %in% genes_chrY$MGI.symbol, ]
n_samples_Y <- rowSums(cpm_Y < 1)
df_Y <- data.frame(
  gene = rownames(cpm_Y),
  n_samples_cpm_inf_1 = n_samples_Y,
  prop_samples_cpm_inf_1 = n_samples_Y/ncol(cpm_Y), 
  pass_filter = (n_samples_Y < 8), 
  stringsAsFactors = FALSE
)
sum(df_Y$n_samples_cpm_inf_1 < 8)
write.xlsx(df_Y, file = file.path(output.path, "/infos_chrXY/genes_chrY_in_male.xlsx"))

counts_female_filtered <- rownames_to_column(counts_female_filtered, "MGI.symbol")
counts_male_filtered <- rownames_to_column(counts_male_filtered, "MGI.symbol")

# joining both datasets 
counts_both <- inner_join(counts_female_filtered, counts_male_filtered, by = "MGI.symbol")
counts_both <- column_to_rownames(counts_both, "MGI.symbol")
print(dim(counts_both))

load(megena_results.path)
load(modtable.path)

genes_chrX <- read.table(
  genes_chrX.path,
  header = TRUE,
  stringsAsFactors = FALSE
)
genelist_chrX <- list()
genelist_chrX[["genes"]] <- genes_chrX$MGI.symbol


module_list <- summary.output$modules

total.genes <- vcount(g) 

# vecteur des gènes chrX
genesX <- genelist_chrX[["genes"]]

# noms des modules
module_ids <- names(module_list)

# calcul pour chaque module
res_list <- lapply(module_ids, function(mod) {
  
  genes_in_module <- module_list[[mod]]
  
  module_size <- length(genes_in_module)
  
  # gènes X dans le module
  genesX_in_module <- genes_in_module[genes_in_module %in% genesX]
  
  n_genesX <- length(genesX_in_module)
  prop_genesX <- n_genesX / module_size
  
  data.frame(
    module.id = mod,
    module_size = module_size,
    n_genesX = n_genesX,
    prop_genesX = prop_genesX,
    genesX_list = paste(genesX_in_module, collapse = ","),  # <- ici
    stringsAsFactors = FALSE
  )
})

# combiner
result_df <- bind_rows(res_list)

# ajouter les infos depuis modtable
result_df <- result_df %>%
  left_join(
    modtable %>%
      select(module.id, module.parent, generation),
    by = "module.id"
  )

print(result_df)

write.xlsx(result_df,
           file =  file.path(output.path,"module_description_XY.xlsx"))

################### 
## Tester si mes genes X dans le reseaux sont DEG 
deg_male <- readRDS("/home/marinevernier/Documents/projets/cpid_multiregion/data/2__differential_expression_analysis/annotated_counts_filtered.rds")


sum(
  deg_male[deg_male$MGI.symbol == "Rpl3-ps1",
           grep("pval", names(deg_male))] < 0.05,
  na.rm = TRUE
)
deg_male[deg_male$MGI.symbol == "Rpl3-ps1",
         grep("pval", names(deg_male))] 


# Colonnes contenant des p-values
pval_cols <- grep("pval", names(deg_male), value = TRUE)

# Tableau récapitulatif
resume_sig_male <- data.frame(
  gene = genes_chrX$MGI.symbol,
  nb_colonnes_significatives_male = sapply(genes_chrX$MGI.symbol, function(g) {
    
    # Sous-table pour le gène
    tmp <- deg_male[deg_male$MGI.symbol == g, pval_cols]
    
    # Compte des pvalues < 0.05
    sum(tmp < 0.05, na.rm = TRUE)
  })
)

resume_sig_male


deg_female <- readRDS("/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts_filtered.rds")
# Colonnes contenant des p-values
pval_cols <- grep("pval", names(deg_female), value = TRUE)

# Tableau récapitulatif
resume_sig_female <- data.frame(
  gene = genes_chrX$MGI.symbol,
  nb_colonnes_significatives_female = sapply(genes_chrX$MGI.symbol, function(g) {
    
    # Sous-table pour le gène
    tmp <- deg_female[deg_female$MGI.symbol == g, pval_cols]
    
    # Compte des pvalues < 0.05
    sum(tmp < 0.05, na.rm = TRUE)
  })
)

resume_sig_female


resume_sig_all <- inner_join(resume_sig_male, resume_sig_female, by = "gene")
resume_sig_all

###############################################
## plotter l'expression des genes X 
raw_counts_female <- read.table(raw_counts_female.path, header = TRUE,
                                sep = "\t",
                                comment.char = "#",
                                stringsAsFactors = FALSE,
                                check.names = FALSE)

#recuperer les gènes des chromosomes X et Y pour plus tard 
# Gènes chrX
genes_chrX <- raw_counts_female %>%
  filter(Chr == "chrX") %>%
  mutate(Geneid = str_replace(Geneid, "\\..*", "")) %>%
  select(Geneid) %>%
  distinct()
genes_chrX <- genes_chrX %>%
  inner_join(ens2symbol, by = c("Geneid" = "Gene.stable.ID")) %>%
  dplyr::select(MGI.symbol) %>%
  distinct()

cpm_female_X <- cpm_female[rownames(cpm_female) %in% genes_chrX$MGI.symbol, ]

cpm_male_X <- cpm_male[rownames(cpm_male) %in% genes_chrX$MGI.symbol, ]

# Moyenne d'expression par gène
female_mean <- rowMeans(cpm_female_X, na.rm = TRUE)
male_mean   <- rowMeans(cpm_male_X, na.rm = TRUE)

# Construire un data.frame
df_plot <- bind_rows(
  data.frame(
    gene = names(female_mean),
    mean_expression = female_mean,
    sex = "Female"
  ),
  data.frame(
    gene = names(male_mean),
    mean_expression = male_mean,
    sex = "Male"
  )
)

# Plot
ggplot(df_plot, aes(x = mean_expression)) +
  geom_histogram(
    bins = 50,
    fill = "steelblue",
    color = "black"
  ) +
  stat_bin(
    bins = 50,
    geom = "text",
    aes(label = ifelse(after_stat(count) > 0,
                       after_stat(count), "")),
    vjust = -0.3,
    size = 3
  ) +
  facet_wrap(~ sex, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Mean gene expression (CPM)",
    y = "Number of genes",
    title = "Distribution of mean gene expression per gene"
  )

# Nombre de samples avec CPM < 1 pour chaque gène
female_detected <- rowSums(cpm_female_X < 1, na.rm = TRUE)
male_detected   <- rowSums(cpm_male_X < 1, na.rm = TRUE)

# Seuil = 20% du nombre total de samples
female_threshold <- 0.20 * ncol(cpm_female_X)
male_threshold   <- 0.20 * ncol(cpm_male_X)

# Data frame
df_detected <- bind_rows(
  data.frame(
    gene = names(female_detected),
    n_samples = female_detected,
    sex = "Female"
  ),
  data.frame(
    gene = names(male_detected),
    n_samples = male_detected,
    sex = "Male"
  )
)

# Data frame pour les lignes verticales
df_thresholds <- data.frame(
  sex = c("Female", "Male"),
  threshold = c(female_threshold, male_threshold)
)

# Plot
ggplot(df_detected, aes(x = n_samples)) +
  geom_histogram(
    bins = 30,
    fill = "steelblue",
    color = "black"
  ) +
  stat_bin(
    bins = 30,
    geom = "text",
    aes(label = ifelse(after_stat(count) > 0,
                       after_stat(count), "")),
    vjust = -0.3,
    size = 3
  ) +
  geom_vline(
    data = df_thresholds,
    aes(xintercept = threshold),
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  facet_wrap(~ sex, scales = "free_x") +
  theme_bw() +
  labs(
    x = "Number of samples with CPM < 1",
    y = "Number of genes",
    title = "Genes weakly expressed across samples"
  )

