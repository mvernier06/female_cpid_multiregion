library(tidyverse)
library(tibble)
library(ggplot2)
library(clusterProfiler)
## ANNOTATIONS DATABASE ##
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

rm(list=ls())

project.path <- "/home/marinevernier/Documents/projets/"
setwd(project.path)

deg_male_path <- "cpid_multiregion/data/2__differential_expression_analysis/deglist.Rdata"
deg_female_path <- "female_cpid_multiregion/data/2__differential_expression_analysis/deglist.Rdata"

output_path <- "female_cpid_multiregion/graphs_results/0_comparison_male_female/intersection_deg/"
dir.create(output_path)

# -------- LOAD MALE --------
male_env <- new.env()
load(deg_male_path, envir = male_env)

male_objs <- ls(male_env, pattern = "^deg_")

for(obj in male_objs){
  assign(paste0(obj, "_male"), male_env[[obj]])
}


# -------- LOAD FEMALE --------
female_env <- new.env()
load(deg_female_path, envir = female_env)

female_objs <- ls(female_env, pattern = "^deg_")

for(obj in female_objs){
  assign(paste0(obj, "_female"), female_env[[obj]])
}

deg_NAc_tp1_female <- deg_Nac_tp1_female
deg_NAc_tp2_female <- deg_Nac_tp2_female
deg_NAc_tp3_female <- deg_Nac_tp3_female

regions <- c("Ins", "NAc", "Hb")
timepoints <- c("tp1", "tp2", "tp3")

get_deg_genes <- function(df, pval_threshold = 0.05, lfc_thresh = 0){
  df %>%
    filter(pval < pval_threshold & abs(log2fc) > lfc_thresh) %>%
    pull(label)
}

intersections_list <- list()

for(tp in timepoints){
  for(reg in regions){
    
    male_name <- paste0("deg_", reg, "_", tp, "_male")
    female_name <- paste0("deg_", reg, "_", tp, "_female")
    
    if(exists(male_name) & exists(female_name)){
      
      male_df <- get(male_name)
      female_df <- get(female_name)
      
      male_genes <- get_deg_genes(male_df)
      female_genes <- get_deg_genes(female_df)
      
      inter <- intersect(male_genes, female_genes)
      
      intersections_list[[paste(reg, tp, sep = "_")]] <- inter
      
    }
  }
}
intersection_sizes <- sapply(intersections_list, length)
intersection_sizes



intersection_df <- tibble(
  comparison = names(intersections_list),
  n_genes = intersection_sizes
) %>%
  tidyr::separate(comparison, into = c("region", "timepoint"), sep = "_")

intersection_df
intersection_df$timepoint <- as.factor(intersection_df$timepoint)
  
ggplot(intersection_df, aes(x=timepoint, y=n_genes, color=region)) +
  geom_point(aes(pch=region)) +
  geom_line(aes(group=region)) +
  labs(title="DEG overlap between male and female",
       y = "Number of overlapping DEG",
       color = "Region", shape = "Region") +
  theme_bw() +
  theme(axis.title.x=element_blank())
ggsave(plot=last_plot(), paste0(output_path, "overlaping_deg_male_female.png"))


## Enrichissement 

for(name in names(intersections_list)){
  message("Doing GO analysis of ", name)
  
  genes <- intersections_list[[name]]
  
  go <- enrichGO(gene = genes, 
                       OrgDb = organism, 
                       keyType = "SYMBOL", 
                       ont ="ALL",
                       pAdjustMethod = "BH", 
                       qvalueCutoff = 1)
  
  assign(paste0("go_", name), go, envir = .GlobalEnv)
}


## SAVE RESULTS ##
go_obj <- ls()[grepl("go_",ls())]
save(list=go_obj, file= paste0(output_path, "go_obj.Rdata"))

## COUNT RESULT ##
for(name in names(intersections_list)){
  go_obj <- paste0("go_", name)
  print(paste0(name, ": ",nrow(get(go_obj))))
}
 
#### REDUCED TERMS FUNCTION #### 
go_rrvgo <- function(intersections_list, ontologies) {
  ## GO ## 
  for (name in names(intersections_list)) {
    for (ont in ontologies) {
      message(paste0("Ontology: ", ont))
      go_name <- paste0("go_", name)
      go_res <- get(go_name, envir = .GlobalEnv)
      go_ont <- go_res@result$ID[go_res@result$ONTOLOGY == ont]
      go_ont_qval <- go_res@result$qvalue[go_res@result$ONTOLOGY == ont]
      
      ## get simMAtrix ## 
      
      if (!is.null(go_ont) && length(go_ont) > 1) {
        simMatrix <- calculateSimMatrix(go_ont, 
                                        orgdb = "org.Mm.eg.db", 
                                        ont = ont, 
                                        method = "Rel")
        if (nrow(simMatrix) > 1) { 
          
          scores <- setNames(-log10(go_ont_qval), go_ont)
          
          reducedTerms <- reduceSimMatrix(simMatrix, 
                                          scores, 
                                          threshold = 0.7,
                                          orgdb = "org.Mm.eg.db")
          
          
          ## SAVE ##
          # PLOT # 
          output_path_treemap <- paste0(output_path, "treemap/", name,"/")
          dir.create(output_path_treemap, recursive = TRUE)
          output_file <- file.path(output_path_treemap, paste0("/treemap_", name, "_", ont, ".png"))
          
          png(output_file, width = 8, height = 8, units = "in", res = 600)
          treemapPlot(reducedTerms)  
          dev.off()
          
          # DATA #
          assign(paste("red_go",name, ont, sep = "_"), reducedTerms, envir = .GlobalEnv)
        }
        
        
      }
    }
  }
}
list_ontology <- c("BP", "CC", "MF")
go_rrvgo(intersections_list,  list_ontology)
