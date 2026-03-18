rm(list = ls()) # rm R working space

# BiocManager::install(c("impute", "preprocessCore"))
# install.packages("WGCNA")

library(WGCNA)
library(tidyverse)

project.path <- "/home/marinevernier/Documents/projets/"
setwd(project.path)

reg <- "Hb"

male_modules.path <- paste0("cpid_multiregion/data/4__MEGENA/MEGENA_results/MEGENA.Results_", reg, ".Rdata")
male_counts.path <- paste0("cpid_multiregion/data/4__MEGENA/CTF_normalized_counts_",reg,".Rdata")

female_modules.path <- paste0("female_cpid_multiregion/data/4__MEGENA/MEGENA.Results_", reg, ".Rdata")
female_counts.path <- paste0("female_cpid_multiregion/data/4__MEGENA/logCPM_RINcorrected_",reg,".Rdata")

# -------- MALES --------
male_env <- new.env()
load(male_modules.path, envir = male_env)

# récupérer objets
MEGENA.output_male  <- male_env$MEGENA.output
summary.output_male <- male_env$summary.output
male_counts <-  load(male_counts.path)
male_counts <- CTF_normalized 

# -------- FEMELLES --------
female_env <- new.env()
load(female_modules.path, envir = female_env)

MEGENA.output_female  <- female_env$MEGENA.output
summary.output_female <- female_env$summary.output
female_counts <- load(female_counts.path)
female_counts <- logCPM_corrected

output_data <- paste0("female_cpid_multiregion/data/4__MEGENA/preservation_vs_male/",reg,"/")
dir.create(output_data)
output_graph <- paste0("female_cpid_multiregion/graphs_results/4__MEGENA/",reg,"/preservation_vs_male/")
dir.create(output_graph)


MEGENA_to_moduleColors <- function(module_list) {
  
  module_list <- module_list[order(sapply(module_list, length))]
  
  moduleColors <- c()
  
  for (mod in names(module_list)) {
    genes <- module_list[[mod]]
    genes <- setdiff(genes, names(moduleColors))
    moduleColors[genes] <- mod
  }
  
  return(moduleColors)
}

moduleColors_male   <- MEGENA_to_moduleColors(MEGENA.output_male$module.output$modules)
moduleColors_female <- MEGENA_to_moduleColors(MEGENA.output_female$module.output$modules)

common_genes <- Reduce(intersect, list(
  names(moduleColors_male),
  names(moduleColors_female),
  rownames(male_counts),
  rownames(female_counts)
))

male_counts   <- male_counts[common_genes, ]
female_counts <- female_counts[common_genes, ]

moduleColors_male   <- moduleColors_male[common_genes]
moduleColors_female <- moduleColors_female[common_genes]

datExpr_male   <- t(male_counts)
datExpr_female <- t(female_counts)

multiExpr <- list(
  male   = list(data = datExpr_male),
  female = list(data = datExpr_female)
)

colorList <- list(male = moduleColors_male)

mp_refMale <- modulePreservation(
  multiExpr,
  colorList,
  referenceNetworks = 1,  # male = référence
  nPermutations = 100,
  randomSeed = 123,
  verbose = 3
)

pres_refMale <- mp_refMale$preservation$Z$ref.male$inColumnsAlsoPresentIn.female
obs_refMale  <- mp_refMale$preservation$observed$ref.male$inColumnsAlsoPresentIn.female

res_refMale <- data.frame(
  module = rownames(pres_refMale),
  Zsummary = pres_refMale[, "Zsummary.pres"],
  medianRank = obs_refMale[, "medianRank.pres"],
  size = obs_refMale[, "moduleSize"],
  density = pres_refMale[, "Zdensity.pres"]
)

head(res_refMale)
filename <- paste0("female_cpid_multiregion/data/4__MEGENA/preservation_vs_male/",reg,"/", reg, "_modulePreservation_refMale.Rdata")
save(res_refMale, file = filename)

ggplot(res_refMale, aes(x = size, y = Zsummary)) +
  geom_point() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  labs(
    title = "Module preservation",
    x = "Module size",
    y = "Zsummary"
  ) +
  theme_minimal()
plot_name <- paste0("female_cpid_multiregion/graphs_results/4__MEGENA/",reg,"/preservation_vs_male/Zsummary_size_refMale.png")
ggsave(plot=last_plot(), filename = plot_name)

ggplot(res_refMale, aes(x = Zsummary)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  labs(
    title = "Distribution du Zsummary (préservation femelle vs mâle)",
    x = "Zsummary",
    y = "Nombre de modules"
  ) +
  theme_minimal()
plot_name <- paste0("female_cpid_multiregion/graphs_results/4__MEGENA/",reg,"/preservation_vs_male/distribution_Zsummary_refMale.png")
ggsave(plot=last_plot(), filename=plot_name)


#############################################################################################################
### tentative avec raw_counts ###
# male_modules.path <- paste0("cpid_multiregion/data/4__MEGENA/MEGENA_results/MEGENA.Results_", reg, ".Rdata")
# male_counts.path <- paste0("cpid_multiregion/data/2__differential_expression_analysis/annotated_counts.csv")
# 
# female_modules.path <- paste0("female_cpid_multiregion/data/4__MEGENA/MEGENA.Results_", reg, ".Rdata")
# female_counts.path <- paste0("female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts.csv")
# 
# # -------- MALES --------
# male_env <- new.env()
# load(male_modules.path, envir = male_env)
# 
# # récupérer objets
# MEGENA.output_male  <- male_env$MEGENA.output
# summary.output_male <- male_env$summary.output
# male_counts <-  read_csv(male_counts.path)
# male_counts <- male_counts %>% select(MGI.symbol, matches(paste0("^", reg, "\\.")))
# 
# # -------- FEMELLES --------
# female_env <- new.env()
# load(female_modules.path, envir = female_env)
# 
# MEGENA.output_female  <- female_env$MEGENA.output
# summary.output_female <- female_env$summary.output
# female_counts <- read_csv(female_counts.path)
# female_counts$Geneid <- NULL
# female_counts <- female_counts %>% select(MGI.symbol, matches(paste0("^", reg, "\\.")))
# 
# male_counts <- male_counts %>%
#   column_to_rownames("MGI.symbol")
# 
# female_counts <- female_counts %>%
#   column_to_rownames("MGI.symbol")
# 
# MEGENA_to_moduleColors <- function(module_list) {
#   
#   module_list <- module_list[order(sapply(module_list, length))]
#   
#   moduleColors <- c()
#   
#   for (mod in names(module_list)) {
#     genes <- module_list[[mod]]
#     genes <- setdiff(genes, names(moduleColors))
#     moduleColors[genes] <- mod
#   }
#   
#   return(moduleColors)
# }
# 
# moduleColors_male   <- MEGENA_to_moduleColors(MEGENA.output_male$module.output$modules)
# moduleColors_female <- MEGENA_to_moduleColors(MEGENA.output_female$module.output$modules)
# 
# common_genes <- Reduce(intersect, list(
#   names(moduleColors_male),
#   names(moduleColors_female),
#   rownames(male_counts),
#   rownames(female_counts)
# ))
# 
# male_counts   <- male_counts[common_genes, ]
# female_counts <- female_counts[common_genes, ]
# 
# moduleColors_male   <- moduleColors_male[common_genes]
# moduleColors_female <- moduleColors_female[common_genes]
# 
# datExpr_male   <- t(male_counts)
# datExpr_female <- t(female_counts)
# 
# multiExpr <- list(
#   male   = list(data = datExpr_male),
#   female = list(data = datExpr_female)
# )
# 
# colorList <- list(male = moduleColors_male)
# 
# mp <- modulePreservation(
#   multiExpr,
#   colorList,
#   referenceNetworks = 1,  # male = référence
#   nPermutations = 100,
#   randomSeed = 123,
#   verbose = 3
# )
# 
# pres <- mp$preservation$Z$ref.male$inColumnsAlsoPresentIn.female
# obs  <- mp$preservation$observed$ref.male$inColumnsAlsoPresentIn.female
# 
# res <- data.frame(
#   module = rownames(pres),
#   Zsummary = pres[, "Zsummary.pres"],
#   medianRank = obs[, "medianRank.pres"],
#   size = obs[, "moduleSize"],
#   density = pres[, "Zdensity.pres"]
# )
# 
# head(res)
# filename <- paste0("female_cpid_multiregion/data/4__MEGENA/preservation_vs_male/", reg, "_modulePreservation_raw_counts.Rdata")
# save(mp, file = filename)
# 
# ggplot(res, aes(x = size, y = Zsummary)) +
#   geom_point() +
#   geom_hline(yintercept = 2, linetype = "dashed") +
#   geom_hline(yintercept = 10, linetype = "dashed") +
#   labs(
#     title = "Module preservation",
#     x = "Module size",
#     y = "Zsummary"
#   ) +
#   theme_minimal()
