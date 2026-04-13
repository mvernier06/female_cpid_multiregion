rm(list = ls()) # rm R working space

library(tidyverse)

setwd("/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/")

# Choose a region : Ins
reg <- "Ins"
print(paste0("Annotating module table for region: ", reg))

#### PATHS ####
megena_results.path <- paste0("data/4__MEGENA/MEGENA_male_female_insula/MEGENA.Results_", reg, ".Rdata")
output.path <- paste0("data/4__MEGENA/MEGENA_male_female_insula/modtable_", reg, ".Rdata")

#### load env ####
load(megena_results.path)
modtable <- summary.output$module.table
module_list <- summary.output$modules

#### get generation number ####
module_generations <- list()
module_generations[["gen1"]] <- modtable$module.id[which(modtable$module.parent == "c1_1")]
modtable$generation <- ifelse(modtable$module.id %in% module_generations[["gen1"]], 1, NA)
gen <- 1

# loop until no more modules are found in the next generation
while(TRUE) {
  current_gen <- list()
  for(i in module_generations[[paste0("gen", gen)]]) {
    current_gen <- append(current_gen, modtable$module.id[which(modtable$module.parent == i)])
  }
  current_gen <- unlist(current_gen)

  if(length(current_gen) == 0) {
    break
  }
  
  gen <- gen + 1
  module_generations[[paste0("gen", gen)]] <- current_gen

  modtable$generation <- ifelse(modtable$module.id %in% current_gen, gen, modtable$generation)
}

cat("Objects before save:\n")
print(ls())

save(modtable, file=output.path)
print(paste0("Module table saved for region: ", reg))