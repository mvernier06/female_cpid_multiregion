library(RRHO2)

rm(list=ls())
setwd("/home/marinevernier/Documents/projets/female_cpid_multiregion/")

#### PATHS ####
rrho_obj.path <- "data/0_comparison_male_female/RRHO/"
output.path <- "graphs_results/0_comparison_male_female/RRHO/"

common_regions <- c("Ins", "Nac", "Hb")
timepoints <- c("tp1", "tp2", "tp3")
reg <- "Ins"
tp <- "tp1"

for (reg in common_regions) {
  print(reg)
  maxreg <- 0
  for (tp in timepoints){
    print(tp)
    filename <- paste0(rrho_obj.path, reg,"_rrho_male_female_",tp,".Rdata")
    load(filename)
    # rrho_name <- paste0("RRHO_femalevsmale_l2fc_pval")
    hypermat <- RRHO_femalevsmale_l2fc_pval$hypermat
    maxvalue <- max(hypermat, na.rm = TRUE)
    if(maxvalue > maxreg){
      maxreg <- maxvalue
    }
    
  }
  cat("\nMax value of region", reg, ":", maxreg)
  assign(paste0("max", reg), maxreg, envir = .GlobalEnv)

}


common_regions <- c("Ins", "Nac", "Hb")
timepoints <- c("tp1", "tp2", "tp3")

for (reg in common_regions) {
  max_pval <- get(paste0("max",reg))
  print(max_pval)
  for (tp in timepoints) {
    
    load(file.path(rrho_obj.path, paste0(reg,"_rrho_male_female_",tp,".Rdata"))) 
    filename <- file.path(output.path,paste0("/",reg,"/scaled_",reg, "_male_vs_female_", tp, ".png"))
    png(file=filename, width = 750, height = 500, units = "px", res = 100)
    RRHO2_heatmap(RRHO_femalevsmale_l2fc_pval, main = paste0(reg, " male vs female ", tp), maximum = max_pval)
    dev.off()
  }
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