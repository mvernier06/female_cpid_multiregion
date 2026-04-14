rm(list = ls()) # rm R working space

library(MEGENA)
library(tidyverse)
library(edgeR)

setwd("/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/")

# Choose a region : Ins
reg <- "Ins"
print(paste0("Running MEGENA for region: ", reg))

# PATHS
pfn.path <- paste0("data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/pfn_", reg, ".Rdata")
output.path <- paste0("data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/MEGENA.Results_", reg, ".Rdata") # file to save results

# import PFN
load(pfn.path)

# input parameters
n.cores <-10; # number of cores/threads to call for PCP
doPar <-TRUE; # do we want to parallelize?
method = "pearson" # method for correlation. either pearson or spearman.
FDR.cutoff = 0.05 # FDR threshold to define significant correlations upon shuffling samples.
module.pval = 0.05 # module significance p-value. Recommended is 0.05.
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs.
hub.perm = 100; # number of permutations for calculating connectivity significance p-value.

# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2
########### 

# perform MCA clustering
MEGENA.output <- do.MEGENA(g,
                           mod.pval = module.pval,
                           hub.pval = hub.pval,
                           remove.unsig = TRUE,
                           min.size = 50,
                           max.size = vcount(g),
                           doPar = doPar,
                           num.cores = n.cores,
                           n.perm = hub.perm,
                           save.output = FALSE)

summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = module.pval,
                                       hub.pvalue = hub.pval,
                                       min.size = 50,
                                       max.size = vcount(g),
                                       annot.table = annot.table,
                                       id.col = id.col,
                                       symbol.col = symbol.col,
                                       output.sig = TRUE)

if (!is.null(annot.table))
{
  # update annotation to map to gene symbols
  V(g)$name <- paste(annot.table[[symbol.col]][match(V(g)$name,annot.table[[id.col]])],V(g)$name,sep = "|")
  summary.output <- output[c("mapped.modules","module.table")]
  names(summary.output)[1] <- "modules"
}

cat("Objects before save:\n")
print(ls())

save(summary.output,
     MEGENA.output,
     g,
     file=output.path)