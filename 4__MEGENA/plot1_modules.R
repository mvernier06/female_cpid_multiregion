# This script is made to be called in the terminal

suppressPackageStartupMessages({
library(MEGENA)
library(tidyverse)
library(ggplot2)
})

setwd("/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/")

#### REGION AS ARGUMENT AFTER THE SCRIPT CALL ####
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("An argument must be supplied at the end of the Rscript call.
       Choose a region (case sensitive: ACC - Hb - Ins - Nac)", call.=FALSE)
}
if (args[1] != "ACC" & args[1] != "Hb" & args[1] != "Ins" & 
    args[1] != "Nac") {
  stop("The argument needs good formatting.
       Choose a region (case sensitive: ACC - Hb - Ins - Nac)", call.=FALSE)
}
reg <- args[1]




#### PATHS ####
MEGENA.Results.path <- paste0("data/4__MEGENA/MEGENA.Results_", reg, ".Rdata")
print(MEGENA.Results.path)
output.path <- paste0("graphs_results/4__MEGENA/", reg) # file to save results
dir.create(output.path, recursive = TRUE, showWarnings = FALSE)



# import megena results
load(MEGENA.Results.path)

setwd(output.path)

module.table <- summary.output$module.table
colnames(module.table)[1] <- "id" # first column of module table must be labelled as "id".

PlotAllModules <- function(){
  
  for (module in module.table$id){
    
    lst <- module.table$module.hub[which(module.table$id==module)] %>%
      strsplit(., ",") %>%
      unlist() %>%
      str_replace(., "\\(.*", "")
    
    pnet.obj <- plot_module(output = summary.output,
                            PFN = g,
                            subset.module = module,
                            layout = "kamada.kawai",
                            label.hubs.only = TRUE,
                            gene.set = list("hub.set" = lst),
                            color.code =  "red",
                            output.plot = TRUE,
                            out.dir = paste0("modulePlot", reg),
                            col.names = c("grey","grey","grey"),
                            hubLabel.col = "black",
                            hubLabel.sizeProp = 1,
                            show.topn.hubs = Inf,
                            show.legend = TRUE)
  }
  
  setwd(paste0("./modulePlot", reg))
  hierarchy.obj <- plot_module_hierarchy(module.table = module.table,label.scaleFactor = 0.15,
                                         arrow.size = 0.03,node.label.color = "blue")
  ggsave(filename="module_hierarchy.png")
  
}

PlotAllModules()