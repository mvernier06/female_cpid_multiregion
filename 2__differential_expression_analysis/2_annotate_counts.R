library(tidyverse)

# PATHS ####
counts.path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/CPID_sham_vs_cuff_betaprior.csv"
annot_table.path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/data/count_data/annotation_final.csv"
output.path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts.csv" # where to save results

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
