library(tidyverse)
library(edgeR)
library(readODS)

rm(list=ls())

#### PATHS ####
raw_counts_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
coldata_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/counts_m39_M32/coldata_without_outliers.ods"
plot_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/leading_logFC"
dir.create(plot_path)
setwd(plot_path)

#### LOAD DATA ####
counts <- read_csv(raw_counts_path) %>% column_to_rownames("MGI.symbol")
coldata <- read_ods(coldata_path)
coldata$reg <- factor(coldata$reg)
coldata$timepoint <- factor(coldata$timepoint)


# Check up
ngenes <- nrow(counts)
coldata$cond <- as.factor(paste(coldata$reg, coldata$group, coldata$timepoint, sep = "_"))
group <- coldata$cond

#### edgeR ####
y <- DGEList(counts, group = group)
y <- calcNormFactors(y)

#### Dimension reduction (MDS) coloured by region ####
png(filename = "leading_logFC_reg.png", height = 750, width = 1250, units = "px", res = 150)
plotMDS(y, top = 200, group = group, cex = 1, pch = 19, 
        col = as.numeric(factor(coldata$reg))) 

col_used <- as.numeric(factor(coldata$reg))
legend("topright", 
       legend = levels(coldata$reg), 
       fill = unique(col_used),  # Utiliser les couleurs uniques
       cex = 0.8)
title(main = "Multidimensional scaling plot of leading logFC colored by regions")
dev.off()

#### Dimension reduction (MDS) coloured by timepoint ####
png(filename = "leading_logFC_tp.png", height = 750, width = 1250, units = "px", res = 150)
plotMDS(y, top = 200, group = group, cex = 1, pch = 19, 
        col = as.numeric(factor(coldata$timepoint))) 
legend("bottomright", 
       legend = paste0("TP", levels(coldata$timepoint)), 
       fill = unique(col_used),
       text.col = 1:5, 
       cex = 0.8)
title(main = "Multidimensional scaling plot of leading logFC colored by TP")
dev.off()
