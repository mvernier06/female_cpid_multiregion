rm(list = ls()) # rm R working space

library(MEGENA)
library(tidyverse)
library(edgeR)

setwd("/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/")

# Choose a region 
reg <- "Hb"

# PATHS ####
CTF_normalized.path <- paste0("data/4__MEGENA/CTF_normalized_counts_", reg, ".Rdata") # CTF normalized counts for the queried region
output.path <- paste0("data/4__MEGENA/ijw_", reg, ".Rdata") # file to save results

# import CTF normalized counts
load(CTF_normalized.path)

# input parameters
cor.perm = 10 # number of permutations for calculating FDRs for all correlation pairs.

# generate correlation matrix
ijw <- calculate.correlation(CTF_normalized,
                             doPerm = cor.perm,
                             output.corTable = FALSE,
                             output.permFDR = FALSE)

# save correlation matrix
save(ijw, file=output.path)