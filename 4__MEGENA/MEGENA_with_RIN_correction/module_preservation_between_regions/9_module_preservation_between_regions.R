rm(list = ls()) # rm R working space

library(WGCNA)
library(tidyverse)

project.path <- "/home2020/home/inci/mvernier/cpid_multireg_female/"
setwd(project.path)

args <- commandArgs(trailingOnly = TRUE)

reg1 <- args[1]
reg2 <- args[2]

cat("Running:", reg1, "vs", reg2, "\n")


# -------- PATHS --------
base_modules_path <- "female_cpid_multiregion/data/4__MEGENA/MEGENA_with_RIN_correction/MEGENA.Results_"
base_counts_path  <- "female_cpid_multiregion/data/4__MEGENA/MEGENA_with_RIN_correction/logCPM_RINcorrected_"

output_data_base  <- "female_cpid_multiregion/data/4__MEGENA/MEGENA_with_RIN_correction/preservation_between_regions/"
output_graph_base <- "female_cpid_multiregion/graphs_results/4__MEGENA/MEGENA_with_RIN_correction/preservation_between_regions/"

# dir.create(output_data_base, recursive = TRUE, showWarnings = FALSE)
# dir.create(output_graph_base, recursive = TRUE, showWarnings = FALSE)

# # -------- FUNCTION --------
# MEGENA_to_moduleColors <- function(module_list) {
#   module_list <- module_list[order(sapply(module_list, length))]
#   moduleColors <- c()
  
#   for (mod in names(module_list)) {
#     genes <- module_list[[mod]]
#     genes <- setdiff(genes, names(moduleColors))
#     moduleColors[genes] <- mod
#   }
  
#   return(moduleColors)
# }

# # -------- LOAD REGIONS --------
# load_region <- function(reg) {
#   env <- new.env()
  
#   load(paste0(base_modules_path, reg, ".Rdata"), envir = env)
#   load(paste0(base_counts_path, reg, ".Rdata"))
  
#   moduleColors <- MEGENA_to_moduleColors(env$MEGENA.output$module.output$modules)
  
#   return(list(
#     moduleColors = moduleColors,
#     counts = logCPM_corrected
#   ))
# }

# data1 <- load_region(reg1)
# data2 <- load_region(reg2)

# # -------- COMMON GENES --------
# common_genes <- Reduce(intersect, list(
#   names(data1$moduleColors),
#   names(data2$moduleColors),
#   rownames(data1$counts),
#   rownames(data2$counts)
# ))
# cat("Number of common genes:", length(common_genes), "\n")

# counts1 <- data1$counts[common_genes, ]
# counts2 <- data2$counts[common_genes, ]

# moduleColors1 <- data1$moduleColors[common_genes]
# moduleColors2 <- data2$moduleColors[common_genes]

# datExpr1 <- t(counts1)
# datExpr2 <- t(counts2)

# multiExpr <- list(
#   reg1 = list(data = datExpr1),
#   reg2 = list(data = datExpr2)
# )

# # ==============================
# #  REG1 → REG2
# # ==============================
# mp1 <- modulePreservation(
#   multiExpr,
#   list(reg1 = moduleColors1),
#   referenceNetworks = 1,
#   nPermutations = 100,
#   randomSeed = 123
# )

# pres1 <- mp1$preservation$Z$ref.reg1$inColumnsAlsoPresentIn.reg2
# obs1  <- mp1$preservation$observed$ref.reg1$inColumnsAlsoPresentIn.reg2

# res1 <- data.frame(
#   module = rownames(pres1),
#   Zsummary = pres1[, "Zsummary.pres"],
#   size = obs1[, "moduleSize"],
#   medianRank = obs1[, "medianRank.pres"],
#   ref = reg1,
#   test = reg2
# )

# # ==============================
# # REG2 → REG1
# # ==============================
# mp2 <- modulePreservation(
#   multiExpr,
#   list(reg2 = moduleColors2),
#   referenceNetworks = 2,
#   nPermutations = 100,
#   randomSeed = 123
# )

# pres2 <- mp2$preservation$Z$ref.reg2$inColumnsAlsoPresentIn.reg1
# obs2  <- mp2$preservation$observed$ref.reg2$inColumnsAlsoPresentIn.reg1

# res2 <- data.frame(
#   module = rownames(pres2),
#   Zsummary = pres2[, "Zsummary.pres"],
#   size = obs2[, "moduleSize"],
#   medianRank = obs2[, "medianRank.pres"],
#   ref = reg2,
#   test = reg1
# )

# # -------- MERGE --------
# df <- rbind(res1, res2)

# # -------- SAVE --------
# save(df, file = paste0(output_data_base, reg1, "_vs_", reg2, ".Rdata"))

load(paste0(output_data_base, reg1, "_vs_", reg2, ".Rdata"))

  
# ---------- PLOTS ----------
p <- ggplot(df, aes(x = Zsummary, fill = ref)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 2, linetype = "dashed") +
    geom_vline(xintercept = 10, linetype = "dashed") +
    theme_minimal() +
    ggtitle(paste("Zsummary:", reg1, "vs", reg2))

ggsave(
filename = paste0(output_graph_base, reg1, "_vs_", reg2, "_density.png"),
plot = p
)

cat("Finished:", reg1, "vs", reg2, "\n")