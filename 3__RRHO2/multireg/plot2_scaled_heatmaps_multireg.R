library(RRHO2)

rm(list=ls())
setwd("/home/marinevernier/Documents/projets/female_cpid_multiregion/")

#### PATHS ####
rrho_obj.path <- "data/3__RRHO2/rrho_obj_multireg.Rdata"
output.path <- "graphs_results/3__RRHO2/multireg_scaled/"

load(rrho_obj.path)
dir.create(output.path, recursive=TRUE, showWarnings = FALSE)
setwd(output.path)

#### get higher pvalue per tp ####
timepoints <- c(1,2,3)
envlist <- ls()
for(tp in timepoints){
  rrholist <- envlist[grepl(tp, envlist)]
  maxtp <- 0
  for(sample in rrholist){
    rrho <- get(sample)
    hypermat <- rrho$hypermat
    maxvalue <- max(hypermat,na.rm=TRUE)
    if(maxvalue > maxtp){
      maxtp <- maxvalue
    }
  }
  cat("\nMax value of TP", tp, ":", maxtp)
  assign(paste0("max", tp), maxtp, envir = .GlobalEnv)
}

#### rrho2 heatmap for all tp ####
list1 <- c("ACC", "Hb", "Ins", "Nac")
list2 <- list1

# homemade loop that iterates each combination of regions
for(tp in c(1,2,3)){
  print(paste0("Computing RRHO2 heatmaps of the TP", tp))
  setwd("/home/marinevernier/Documents/projets/female_cpid_multiregion/")
  setwd(output.path)
  directory <- paste0("./TP", tp)
  dir.create(directory, showWarnings = FALSE)
  setwd(directory)
  max_pval <- get(paste0("max", tp))
  i=2
  for(reg in list1[1:length(list1)-1]){
    for(reg2 in list2[i:length(list2)]){
      rrho <- get(paste0("rrho", reg, "vs", reg2, tp))
      png(file=paste0(reg, "vs", reg2, tp, ".png"), width = 750, height = 500, units = "px", res = 100)
      RRHO2_heatmap(rrho, main = paste0(reg, " vs ", reg2, " tp", tp), maximum = max_pval)
      dev.off()
    }
    i=i+1
  }
}
