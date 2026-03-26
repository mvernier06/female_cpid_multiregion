rm(list = ls()) # rm R working space

# BiocManager::install(c("impute", "preprocessCore"))
# install.packages("WGCNA",repos = "https://cloud.r-project.org")

library(WGCNA)
library(tidyverse)

project.path <- "/home2020/home/inci/mvernier/cpid_multireg_female/"
setwd(project.path)

reg <- "Nac" # ACC - Hb - Ins - Nac

print(paste0("Région : ", reg))

with_RIN_correction_modules.path <- paste0("female_cpid_multiregion/data/4__MEGENA/MEGENA.Results_", reg, ".Rdata")
with_RIN_correction_counts.path <- paste0("female_cpid_multiregion/data/4__MEGENA/logCPM_RINcorrected_",reg,".Rdata")

without_RIN_correction_modules.path <- paste0("female_cpid_multiregion/data/4__MEGENA/MEGENA_without_RIN_correction/MEGENA.Results_", reg, ".Rdata")
without_RIN_correction_counts.path <- paste0("female_cpid_multiregion/data/4__MEGENA/MEGENA_without_RIN_correction/CTF_normalized_counts_",reg,".Rdata")

# -------- with RIN correction --------
with <- new.env()
load(with_RIN_correction_modules.path, envir = with)

# récupérer objets
MEGENA.output_with  <- with$MEGENA.output
summary.output_with <- with$summary.output
with_RIN_correction_counts <-  load(with_RIN_correction_counts.path)
with_RIN_correction_counts <- logCPM_corrected 

# -------- without RIN correction --------
without <- new.env()
load(without_RIN_correction_modules.path, envir = without)

# récupérer objets
MEGENA.output_without  <- without$MEGENA.output
summary.output_without <- without$summary.output
without_RIN_correction_counts <-  load(without_RIN_correction_counts.path)
without_RIN_correction_counts <- CTF_normalized
without_RIN_correction_counts <- load(without_RIN_correction_counts.path)
without_RIN_correction_counts <- CTF_normalized

output_data <- paste0("female_cpid_multiregion/data/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN/",reg,"/")
dir.create(output_data, recursive = TRUE)
print(paste0("Data will be saved in : ", output_data))
output_graph <- paste0("female_cpid_multiregion/graphs_results/4__MEGENA/MEGENA_without_RIN_correction/",reg,"/preservation_with_without_RIN/")
dir.create(output_graph, recursive = TRUE)
print(paste0("Graphs will be saved in : ", output_graph))


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

moduleColors_with  <- MEGENA_to_moduleColors(MEGENA.output_with$module.output$modules)
moduleColors_without <- MEGENA_to_moduleColors(MEGENA.output_without$module.output$modules)

common_genes <- Reduce(intersect, list(
  names(moduleColors_with),
  names(moduleColors_without),
  rownames(with_RIN_correction_counts),
  rownames(without_RIN_correction_counts)
))

with_RIN_correction_counts   <- with_RIN_correction_counts[common_genes, ]
without_RIN_correction_counts <- without_RIN_correction_counts[common_genes, ]

moduleColors_with  <- moduleColors_with[common_genes]
moduleColors_without <- moduleColors_without[common_genes]

datExpr_with  <- t(with_RIN_correction_counts)
datExpr_without <- t(without_RIN_correction_counts)

multiExpr <- list(
  with  = list(data = datExpr_with),
  without = list(data = datExpr_without)
)

colorList <- list(with = moduleColors_with)

print("Calcul de la préservation des modules avec with comme référence")
mp_refWith <- modulePreservation(
  multiExpr,
  colorList,
  referenceNetworks = 1,  # with = référence : on teste si les gènes des modules with sont aussi connectés dans les counts without 
  nPermutations = 100,
  randomSeed = 123,
  verbose = 3
)

pres_refWith <- mp_refWith$preservation$Z$ref.with$inColumnsAlsoPresentIn.without
obs_refWith  <- mp_refWith$preservation$observed$ref.with$inColumnsAlsoPresentIn.without

res_refWith <- data.frame(
  module = rownames(pres_refWith),
  Zsummary = pres_refWith[, "Zsummary.pres"],
  medianRank = obs_refWith[, "medianRank.pres"],
  size = obs_refWith[, "moduleSize"],
  density = pres_refWith[, "Zdensity.pres"]
)

head(res_refWith)
filename <- paste0("female_cpid_multiregion/data/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN/",reg,"/", reg, "_modulePreservation_refWith.Rdata")
save(res_refWith, file = filename)

ggplot(res_refWith, aes(x = size, y = Zsummary)) +
  geom_point() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  labs(
    title = "Module preservation",
    x = "Module size",
    y = "Zsummary"
  ) +
  theme_minimal()
plot_name <- paste0("female_cpid_multiregion/graphs_results/4__MEGENA/MEGENA_without_RIN_correction/",reg,"/preservation_with_without_RIN/Zsummary_size_refWith.png")
ggsave(plot=last_plot(), filename = plot_name)

ggplot(res_refWith, aes(x = Zsummary)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  labs(
    title = "Distribution du Zsummary (with vs without RIN correction)",
    x = "Zsummary",
    y = "Nombre de modules"
  ) +
  theme_minimal()
plot_name <- paste0("female_cpid_multiregion/graphs_results/4__MEGENA/MEGENA_without_RIN_correction/",reg,"/preservation_with_without_RIN/distribution_Zsummary_refWith.png")
ggsave(plot=last_plot(), filename=plot_name)

#########################
## Now with without rin ref ##
#########################
print("Calcul de la préservation des modules avec without comme référence")
colorList_without <- list(without = moduleColors_without)

mp_refwithout <- modulePreservation(
  multiExpr,
  colorList_without,
  referenceNetworks = 2, # without reference : donc on teste si les gènes des modules without sont aussi connecté dans les counts with 
  nPermutations = 100, 
  randomSeed = 123,
  verbose = 3
)

pres_refwithout <- mp_refwithout$preservation$Z$ref.without$inColumnsAlsoPresentIn.with
obs_refwithout  <- mp_refwithout$preservation$observed$ref.without$inColumnsAlsoPresentIn.with

res_refwithout <- data.frame(
  module = rownames(pres_refwithout),
  Zsummary = pres_refwithout[, "Zsummary.pres"],
  medianRank = obs_refwithout[, "medianRank.pres"],
  size = obs_refwithout[, "moduleSize"],
  density = pres_refwithout[, "Zdensity.pres"]
)

head(res_refwithout)
filename <- paste0("female_cpid_multiregion/data/4__MEGENA/MEGENA_without_RIN_correction/preservation_with_without_RIN/",reg,"/", reg, "_modulePreservation_refwithout.Rdata")
save(res_refwithout, file = filename)

ggplot(res_refwithout, aes(x = size, y = Zsummary)) +
  geom_point() +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  labs(
    title = "Module preservation",
    x = "Module size",
    y = "Zsummary"
  ) +
  theme_minimal()
plot_name <- paste0("female_cpid_multiregion/graphs_results/4__MEGENA/MEGENA_without_RIN_correction/",reg,"/preservation_with_without_RIN/Zsummary_size_refwithout.png")
ggsave(plot=last_plot(), filename = plot_name)

ggplot(res_refwithout, aes(x = Zsummary)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  labs(
    title = "Distribution du Zsummary (préservation sans RIN vs avec RIN)",
    x = "Zsummary",
    y = "Nombre de modules"
  ) +
  theme_minimal()
plot_name <- paste0("female_cpid_multiregion/graphs_results/4__MEGENA/MEGENA_without_RIN_correction/",reg,"/preservation_with_without_RIN/distribution_Zsummary_refwithout.png")
ggsave(plot=last_plot(), filename=plot_name)

res_refwith$sex_direction <- "With reference"
res_refwithout$sex_direction <- "Without reference"

df <- rbind(res_refwith, res_refwithout)

ggplot(df, aes(x = Zsummary, fill = sex_direction)) +
  geom_density(alpha = 0.4) +
  labs(title = "Comparaison des distributions de Zsummary") +
  theme_minimal()
