rm(list = ls()) # rm R working space

library(MEGENA)

setwd("/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/")

# Choose a region : Ins
reg <- "Ins"
print(paste("Doing PFN for region:", reg))

# PATHS
ijw.path <- paste0("data/4__MEGENA/MEGENA_male_female_insula/without_chr_X_Y/ijw_", reg, ".Rdata")
output.path <- paste0("data/4__MEGENA/MEGENA_male_female_insula/without_chr_X_Y/pfn_", reg, ".Rdata") # file to save results

# import correlation matrix
load(ijw.path)

# parameters
doPar <- TRUE
n.cores <- 32

#### register multiple cores if needed: note that set.parallel.backend() is deprecated. 
run.par = doPar & (getDoParWorkers() == 1)
if (run.par)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}

# calculate PFN
el <- calculate.PFN(ijw[,1:3],
                    doPar = TRUE, # do we want to parallelize ?
                    num.cores = n.cores,
                    keep.track = FALSE)

# graph dataframe
g <- graph.data.frame(el,directed = FALSE)

cat("Objects before save:\n")
print(ls())

# save PFN
save(g, file=output.path)