library(RRHO2)

rm(list=ls())
setwd("/home/marinevernier/Documents/projets/female_cpid_multiregion/")

#### PATHS ####
rrho_obj.path <- "data/3__RRHO2/rrho_obj_multitp.Rdata"
output.path <- "graphs_results/3__RRHO2/multitp_classic/"
dir.create(output.path, showWarnings = FALSE, recursive = TRUE)

#### classic rrho plots ####
load(rrho_obj.path)
setwd(output.path)


regions <- c("ACC", "Hb", "Ins", "Nac")
comparisons <- c("1vs2", "2vs3", "1vs3")

for (reg in regions) {
  setwd("/home/marinevernier/Documents/projets/female_cpid_multiregion/")
  setwd(output.path)
  directory <- paste0("./", reg )
  dir.create(directory, showWarnings = FALSE)
  setwd(directory)
  for (comp in comparisons) {
    
    # Nom de l'objet RRHO
    obj_name <- paste0("RRHO_", reg, "_", comp)
    
    # Récupérer l'objet
    rrho_obj <- get(obj_name)
    
    # Nom du fichier
    file_name <- paste0(reg, "_", comp, ".png")
    
    # Titre du plot
    plot_title <- paste(reg, comp)
    
    png(file = file_name,
        width = 750,
        height = 500,
        units = "px",
        res = 100)
    
    RRHO2_heatmap(rrho_obj, main = plot_title)
    
    dev.off()
  }
}
