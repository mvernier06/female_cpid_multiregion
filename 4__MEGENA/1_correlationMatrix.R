rm(list = ls()) # rm R working space

library(MEGENA)
library(tidyverse)
library(edgeR)

setwd("/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/")

# Choose a region (ACC, Hb, Ins, Nac)
reg <- "Hb"

# PATHS ####
logCPM_RINcorrected.path <- paste0("data/4__MEGENA/logCPM_RINcorrected_", reg, ".Rdata") # CTF normalized counts for the queried region
output.path <- paste0("data/4__MEGENA/ijw_", reg, ".Rdata") # file to save results

# import CTF normalized counts
load(logCPM_RINcorrected.path)

# input parameters
cor.perm = 10 # number of permutations for calculating FDRs for all correlation pairs.

# generate correlation matrix
ijw <- calculate.correlation(logCPM_corrected,
                             doPerm = cor.perm,
                             output.corTable = FALSE,
                             output.permFDR = FALSE)

# save correlation matrix
save(ijw, file=output.path)