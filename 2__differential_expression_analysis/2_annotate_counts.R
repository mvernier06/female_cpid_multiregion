library(tidyverse)

setwd("/home/marinevernier/Documents/projets/")

# PATHS ####
counts.path <- "female_cpid_multiregion/data/2__differential_expression_analysis/CPID_sham_vs_cuff_betaprior.csv"
annot_table.path <- "female_cpid_multiregion/data/counts_m39_M32/annotation_final.csv"
output.path <- "female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts.csv" # where to save results

# Annotate counts
counts <- read.csv(counts.path)
counts$Geneid <-  str_replace(counts$Geneid, "\\..*", "") 

ens2symbol <- read.csv(annot_table.path, sep=";")
ens2symbol <- ens2symbol %>% dplyr::select(!X)

annotated_counts <- inner_join(counts, ens2symbol, by=c("Geneid"="Gene.stable.ID"))
annotated_counts <- annotated_counts %>% relocate(MGI.symbol, .before=Geneid) 

# remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)
annotated_counts <- annotated_counts[!duplicated(annotated_counts$MGI.symbol),]

write_csv(annotated_counts, output.path)
