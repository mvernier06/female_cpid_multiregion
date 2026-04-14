rm(list = ls()) # rm R working space

library(MEGENA)
library(tidyverse)
library(edgeR)

setwd("/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/")

# Choose a region : Ins
reg <- "Ins"
print(paste("Doing correlation matrix for region:", reg))

# PATHS ####
CTF_normalized.path <- paste0("data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/CTF_normalized_counts_", reg, ".Rdata") # CTF normalized counts for the queried region
output.path <- paste0("data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/ijw_", reg, ".Rdata") # file to save results

# import CTF normalized counts
load(CTF_normalized.path)

# input parameters
cor.perm = 10 # number of permutations for calculating FDRs for all correlation pairs.

# generate correlation matrix
ijw <- calculate.correlation(CTF_normalized,
                             doPerm = cor.perm,
                             output.corTable = FALSE,
                             output.permFDR = FALSE)
cat("Objects before save:\n")
print(ls())
# save correlation matrix
save(ijw, file=output.path)