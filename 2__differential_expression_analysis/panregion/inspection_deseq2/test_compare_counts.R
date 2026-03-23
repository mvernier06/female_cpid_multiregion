rm(list=ls())
#### LIBRAIRIES ####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(conflicted)
library(readODS)

#### PATHS #### 
project.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project.path)

go_tp3 <- "data/2__differential_expression_analysis/panregion/go/go_enrich/go_enrich_panregion_tp3.rds"
deg_tp3 <- "data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/deseq2_panregion/panregion_betaprior_rg_tp3.rds"

metadata.path <- "data/counts_m39_M32/coldata_without_outliers.ods"

plot.path <- "graphs_results/panregion/correction_counts/myelin_sheath_genes_counts/counts_myelin_sheath_cuff_tp3.png"


#### PREPROCESS ####

deg_tp3 <- readRDS(deg_tp3)

### MODIFICATIONS METADATA ###
metadata <- read_ods(metadata.path)
rownames(metadata) <- metadata$sample


## CHOSE RELATED TP ##
metadata_tp <- metadata[which(metadata$timepoint == 3),] %>% arrange(sample)   ## TO MODIFY
list <- metadata_tp$sample[which(metadata_tp$timepoint == 3 & metadata_tp$group == "cuff")]
## TO MODIFY
counts_tp3 <- deg_tp3[which(colnames(deg_tp3) %in% list,)]

deg_tp3_list <- as.data.frame(deg_tp3$MGI.symbol)
colnames(deg_tp3_list) <- "MGI.symbol"
counts_tp3 <- bind_cols(deg_tp3_list, counts_tp3)

#### GO ####
go_tp3 <- as.data.frame(readRDS(go_tp3))

# Filtrer les termes GO contenant "myelin sheath"
myelin_df <- go_tp3[grep("myelin sheath", go_tp3$Description, ignore.case = TRUE), ]

# Extraire les gènes associés
genes_list <- unique(unlist(strsplit(myelin_df$geneID, "/")))

# Afficher la liste des gènes
print(genes_list)

counts_tp3 <- counts_tp3[which(counts_tp3$MGI.symbol %in% genes_list), ]


# Extraire les noms des régions à partir des noms de colonnes
regions <- sub("\\..*", "", colnames(counts_tp3)[-1])  # On enlève 'MGI.symbol' et garde les préfixes avant le point

# Créer une liste unique des régions
unique_regions <- unique(regions)

# Initialiser un DataFrame pour stocker les résultats
result <- data.frame(MGI.symbol = counts_tp3$MGI.symbol)

# Pour chaque région, sommer les counts des échantillons correspondants
for (region in unique_regions) {
  # Identifier les colonnes de cette région
  region_cols <- grep(paste0("^", region, "\\."), colnames(counts_tp3), value = TRUE)
  
  # Faire la moyenne des counts par gène
  result[[region]] <- rowMeans(counts_tp3[, region_cols, drop = FALSE])
}

# Affichage du résultat
print(result)

# Transformation en format long
result_long <- result %>%
  pivot_longer(
    cols = -MGI.symbol,
    names_to = "Region",
    values_to = "Counts"
  )

# Line plot pour tous les gènes sur un seul graphique
plot <-ggplot(result_long, aes(x = Region, y = Counts, group = MGI.symbol, color = MGI.symbol)) +
  geom_line(linewidth = 1) +                      # Courbe pour chaque gène
  geom_point(size = 2) +                     # Points pour marquer les valeurs exactes
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Counts par région pour chaque gène (TP3 - cuff)",
       x = "Régions",
       y = "Counts",
       color = "Gènes  (myelin sheath)") + 
  scale_y_continuous(limits = c(0,70000))

plot

ggsave(plot.path, plot = plot, bg="white", width=1500, height=1000, units="px", scale=2)
