# This script is made to be called in the terminal

rm(list = ls()) # rm R working space

suppressPackageStartupMessages({
  library(tidyverse)
  library(GeneOverlap)
  library(igraph)
  library(ggraph)
})

install.packages("pheatmap", repos = "https://cloud.r-project.org")
install.packages("RColorBrewer", repos = "https://cloud.r-project.org")

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

project_dir <- "/home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/"
setwd(project_dir)

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
modtable.path <- paste0("data/4__MEGENA/MEGENA_without_RIN_correction/modtable_", reg, ".Rdata")
alluvial_enrichment.path <- paste0("data/4__MEGENA/MEGENA_without_RIN_correction/enrichment_alluvial/enrichment_alluvial_", reg, ".Rdata")
cell_types_enrichment <- paste0("data/4__MEGENA/MEGENA_without_RIN_correction/enrichment_cell_types/enrichment_cell_types_", reg, ".Rdata")
plots.path <- paste0("graphs_results/4__MEGENA/MEGENA_without_RIN_correction/", reg)
celltypes_plots.path <- paste0(plots.path, "/cell_types_network_plots_", reg)
output.path <- paste0("data/4__MEGENA/MEGENA_without_RIN_correction/enrichment_cell_types/")

dir.create(plots.path, recursive = TRUE, showWarnings = FALSE)
dir.create(celltypes_plots.path)



#### NETWORK PLOT ####
load(modtable.path)

# function to reset the network graph
init_network <- function(){
  edges <- modtable %>% dplyr::select(module.parent, module.id)
  rownames(edges) <- NULL
  colnames(edges) <- c("from", "to")
  edges$from <- str_replace(edges$from, "c1_", "M")
  edges$to <- str_replace(edges$to, "c1_", "M")
  assign("edges", edges, envir = .GlobalEnv)
  
  nodes <- modtable %>% dplyr::select(module.id, generation)
  rownames(nodes) <- NULL
  colnames(nodes) <- c("name", "generation")
  nodes$name <- str_replace(nodes$name, "c1_", "M")
  nodes <- rbind(c("M1",0), nodes)
  assign("nodes", nodes, envir = .GlobalEnv)
}
init_network()
#head(nodes)
#head(edges)

# Check for mismatches in edge list
mismatched_edges <- setdiff(unique(c(edges$from, edges$to)), nodes$name)
if(length(mismatched_edges) > 0) {
  cat("Mismatched vertex names in edge list:", mismatched_edges, "\n")
  # Optionally, you could filter out these edges or correct them
  edges <- edges[edges$from %in% nodes$name & edges$to %in% nodes$name, ]
}
#head(edges)

# Create the graph object
flareGraph <- graph_from_data_frame(edges, vertices = nodes)

# Plot the graph using ggraph with circular layout
title <- paste0("Network graph of ", reg, " modules")
p <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = TRUE) + 
  geom_edge_diagonal(aes(alpha = after_stat(index))) +
  coord_fixed() + 
  scale_edge_alpha('Direction', guide = 'none') +
  geom_node_point(size = 1.5) +
  scale_color_identity() +
  geom_node_text(aes(label = name), repel = TRUE, na.rm = TRUE, size = 2.2, max.overlaps=100, color="grey20") +  # Add labels for colored nodes
  ggforce::theme_no_axes() +
  ggtitle(title)
setwd(plots.path)
ggsave(p, filename=paste0("network_plot_", reg, ".png"))

title <- paste0("Network graph of ", reg, " modules")
p_dendrogram <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = FALSE) + 
  geom_edge_diagonal(aes(alpha = after_stat(index))) +
  scale_edge_alpha('Direction', guide = 'none') +
  geom_node_point(size = 1.5) +
  scale_color_identity() +
  ggforce::theme_no_axes() +
  ggtitle(title)
ggsave(p_dendrogram, filename=paste0("network_plot_dendrogram_", reg, ".png"))



setwd(project_dir)
#### network plot modules enriched ns-down-ns and ns-up-ns ####
load(alluvial_enrichment.path)

init_network()
#head(nodes)
#head(edges)

# Check for mismatches in edge list
mismatched_edges <- setdiff(unique(c(edges$from, edges$to)), nodes$name)
if(length(mismatched_edges) > 0) {
  cat("Mismatched vertex names in edge list:", mismatched_edges, "\n")
  # Optionally, you could filter out these edges or correct them
  edges <- edges[edges$from %in% nodes$name & edges$to %in% nodes$name, ]
}
#head(edges)

# patterns to color
t <- merge(alluvial_enrichment_ns_DOWN_ns, 
           alluvial_enrichment_ns_UP_ns, by=0) %>% 
  select(Row.names, overlap.pval.x, overlap.pval.y) %>% 
  rename(c(modules = Row.names, pval_ndn = overlap.pval.x, pval_nun = overlap.pval.y)) %>%
  mutate(signif = case_when(pval_ndn < 0.05 ~ "DOWN",
                            pval_nun < 0.05 ~ "UP",
                            pval_ndn < 0.05 & pval_nun < 0.05 ~ "BOTH",
                            TRUE ~ NA)) 
t[sort(rownames(t)),]
t$modules <- str_replace(t$modules, "c1_", "M")
t$modules <- str_replace(t$modules, ".genes", "")
#head(t)
mod_nun_enriched <- t$modules[which(t$signif=="UP")] 
mod_ndn_enriched <- t$modules[which(t$signif=="DOWN")] 
mod_ndn_nun_enriched <- t$modules[which(t$signif=="BOTH")] 

# Assign colors to nodes based on membership in lists
nodes$color <- NA
nodes$color[nodes$name %in% mod_ndn_enriched] <- "blue"
nodes$color[nodes$name %in% mod_nun_enriched] <- "red"

nodes$label <- ifelse(is.na(nodes$color), NA, nodes$name)
nodes$signif <- ifelse(nodes$color)

# Create the graph object
flareGraph <- graph_from_data_frame(edges, vertices = nodes)

# Plot the graph using ggraph with circular layout
setwd(plots.path)
title <- paste0("Network graph of ", reg, " modules")
p2 <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = TRUE) + 
  geom_edge_diagonal(aes(alpha = after_stat(index))) +
  coord_fixed() + 
  scale_edge_alpha('Direction', guide = 'none') +
  geom_node_point(aes(color = color), size = 2) +
  geom_node_text(aes(label = label), repel = TRUE, na.rm = TRUE, size = 2.5, max.overlaps=100, color="grey20") +  # Add labels for colored nodes
  scale_color_identity() +
  ggforce::theme_no_axes() +
  ggtitle(title)
ggsave(p2, filename=paste0("network_plot_colored_", reg, ".png"))

p2_dendrogram <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = FALSE) + 
  geom_edge_diagonal(aes(alpha = after_stat(index))) +
  scale_edge_alpha('Direction', guide = 'none') +
  geom_node_point(aes(color = color), size = 2) +
  geom_node_text(aes(label = label), repel = TRUE, na.rm = TRUE, size = 2.5, max.overlaps=100, color="grey20") +  # Add labels for colored nodes
  scale_color_identity() +
  ggforce::theme_no_axes() +
  ggtitle(title)
ggsave(p2_dendrogram, filename=paste0("network_plot_colored_dendrogram_", reg, ".png"))



setwd(project_dir)
#### network plot modules enriched with cellular types from Mithil ####
load(cell_types_enrichment)

# 1) Astrocytes:
t <- Astrocytes_more_than_5_fpkm %>% rownames_to_column("Row.names") %>%
  select(Row.names, overlap.pval) %>% 
  rename(c(modules = Row.names, pval = overlap.pval)) %>%
  mutate(signif = case_when(pval < 0.05 ~ "enriched",
                            TRUE ~ NA)) 
t$modules <- str_replace(t$modules, "c1_", "M")
t$modules <- str_replace(t$modules, ".genes", "")
a <- data.frame("modules"=t$modules[which(t$signif=="enriched")],
                                 "astrocytes_pval" = unlist(t$pval[which(t$signif=="enriched")]))

# 2) Endothelial cells
t <- Endothelial.Cells_more_than_5_fpkm %>% rownames_to_column("Row.names") %>%
  select(Row.names, overlap.pval) %>% 
  rename(c(modules = Row.names, pval = overlap.pval)) %>%
  mutate(signif = case_when(pval < 0.05 ~ "enriched",
                            TRUE ~ NA)) 
t$modules <- str_replace(t$modules, "c1_", "M")
t$modules <- str_replace(t$modules, ".genes", "")
e <- data.frame("modules"=t$modules[which(t$signif=="enriched")],
                                   "endothelial_cells_pval" = unlist(t$pval[which(t$signif=="enriched")]))

# 3) Microglia
t <- Microglia_more_than_5_fpkm %>% rownames_to_column("Row.names") %>%
  select(Row.names, overlap.pval) %>% 
  rename(c(modules = Row.names, pval = overlap.pval)) %>%
  mutate(signif = case_when(pval < 0.05 ~ "enriched",
                            TRUE ~ NA)) 
t$modules <- str_replace(t$modules, "c1_", "M")
t$modules <- str_replace(t$modules, ".genes", "")
mi <- data.frame("modules"=t$modules[which(t$signif=="enriched")],
                                   "microglia_pval" = unlist(t$pval[which(t$signif=="enriched")]))

# 4) Myelinating oligo
t <- Myelinating.Oligodendrocytes_more_than_5_fpkm %>% rownames_to_column("Row.names") %>%
  select(Row.names, overlap.pval) %>% 
  rename(c(modules = Row.names, pval = overlap.pval)) %>%
  mutate(signif = case_when(pval < 0.05 ~ "enriched",
                            TRUE ~ NA)) 
t$modules <- str_replace(t$modules, "c1_", "M")
t$modules <- str_replace(t$modules, ".genes", "")
my <- data.frame("modules"=t$modules[which(t$signif=="enriched")],
                               "myelinating_oligodendrocytes_pval" = unlist(t$pval[which(t$signif=="enriched")]))

# 5) Neurons
t <- Neuron_more_than_5_fpkm %>% rownames_to_column("Row.names") %>%
  select(Row.names, overlap.pval) %>% 
  rename(c(modules = Row.names, pval = overlap.pval)) %>%
  mutate(signif = case_when(pval < 0.05 ~ "enriched",
                            TRUE ~ NA)) 
t$modules <- str_replace(t$modules, "c1_", "M")
t$modules <- str_replace(t$modules, ".genes", "")
n <- data.frame("modules"=t$modules[which(t$signif=="enriched")],
                                 "neurons_pval" = unlist(t$pval[which(t$signif=="enriched")]))


cell_types_enrich_modules <- list(a, e, mi, my, n) %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="modules"), .)

cell_types_enrich_modules[,-1] <- cell_types_enrich_modules[,-1] %>% 
  apply(2, function(x) -log10(x))

cell_types_enrich_modules <- cell_types_enrich_modules %>% column_to_rownames("modules")

print(dir.exists(output.path))
write_rds(cell_types_enrich_modules, paste0(output.path, "modules_enriched_for_each_cell_type_", reg, ".rds"))


p <- pheatmap(cell_types_enrich_modules,
         color = rev(colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "OrRd")))(100)),
         display_numbers = TRUE,
         fontsize = 8,
         angle = 315,
         cluster_cols = FALSE,
         cluster_rows = FALSE)
setwd(project_dir)
setwd(plots.path)
ggsave(paste0("/cell_types_table_", reg, ".png"), p, height = 2500, width = 3000, units = "px")

#### dendrograms avec couleur pour type cellulaire (1 graph par type) ####
init_network()

network_plot <- function(enriched_data, cell_type){
  # patterns to color
  t <- enriched_data[,-4] %>% 
    select(overlap.pval) %>% 
    rownames_to_column("modules") %>%
    mutate(signif = case_when(overlap.pval < 0.05 ~ "enriched",
                              TRUE ~ NA)) 
  t$modules <- str_replace(t$modules, "c1_", "M")
  t$modules <- str_replace(t$modules, ".genes", "")
  #head(t)
  mod_enriched <- t$modules[which(t$signif=="enriched")] 

  # Assign colors to nodes based on membership in lists
  nodes$color <- NA
  nodes$color[nodes$name %in% mod_enriched] <- case_when(cell_type == "astrocytes" ~  "yellow3",
                                                        cell_type == "endothelial_cells" ~  "red",
                                                        cell_type == "microglia" ~  "pink",
                                                        cell_type == "myelinating_oligodendrocytes" ~ "purple",
                                                        cell_type == "neurons" ~ "blue")
  
  nodes$label <- ifelse(is.na(nodes$color), NA, nodes$name)
  nodes$signif <- ifelse(nodes$color)
  
  # Create the graph object
  flareGraph <- graph_from_data_frame(edges, vertices = nodes)
  
  # Plot the graph using ggraph with circular layout
  p <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = TRUE) + 
    geom_edge_diagonal(aes(alpha = after_stat(index))) +
    coord_fixed() + 
    scale_edge_alpha('Direction', guide = 'none') +
    geom_node_point(aes(color = color), size = 2) +
    geom_node_text(aes(label = label), repel = TRUE, na.rm = TRUE, size = 2.5, max.overlaps=100, color="grey20") +  # Add labels for colored nodes
    scale_color_identity() +
    ggforce::theme_no_axes() +
    ggtitle(paste0("Network graph of ", reg, " modules enriched in ", cell_type))
  setwd(project_dir)
  setwd(celltypes_plots.path)
  ggsave(p, filename=paste0("modules_", reg, "_", cell_type, ".png"))
  
  p_dendrogram <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = FALSE) + 
    geom_edge_diagonal(aes(alpha = after_stat(index))) +
    scale_edge_alpha('Direction', guide = 'none') +
    geom_node_point(aes(color = color), size = 2) +
    geom_node_text(aes(label = label), repel = TRUE, na.rm = TRUE, size = 2.5, max.overlaps=100, color="grey20") +  # Add labels for colored nodes
    scale_color_identity() +
    ggforce::theme_no_axes() +
    ggtitle(paste0("Network graph of ", reg, "modules enriched in ", cell_type))
  ggsave(p_dendrogram, filename=paste0("modules_", reg, "_", cell_type, "_d", ".png"))
}

network_plot(Astrocytes_more_than_5_fpkm, "astrocytes")
network_plot(Endothelial.Cells_more_than_5_fpkm, "endothelial_cells")
network_plot(Microglia_more_than_5_fpkm, "microglia")
network_plot(Myelinating.Oligodendrocytes_more_than_5_fpkm, "myelinating_oligodendrocytes")
network_plot(Neuron_more_than_5_fpkm, "neurons")


## Network enrichis en DEG 
setwd(project_dir)
deg_enrichment.path <- paste0("data/4__MEGENA/MEGENA_without_RIN_correction/enrichment_DEGs/enrichment_degs_", reg, ".Rdata")
load(deg_enrichment.path)
plots.path <-paste0("graphs_results/4__MEGENA/MEGENA_without_RIN_correction/",reg, "/DEG_network_plot_", reg, "/")
modtable.path <- paste0("data/4__MEGENA/MEGENA_without_RIN_correction/modtable_", reg, ".Rdata")
load(modtable.path)
tplist <- c("tp1","tp2","tp3")

for(tp in tplist){
  
  print(paste("Processing", tp))
  
  init_network()
  
  # récupérer table enrichissement
  df <- get(paste0("deg_enrichment_", reg, "_", tp))
  
  df <- df %>%
    tibble::rownames_to_column("module") %>%
    mutate(
      module = str_replace(module, "c1_", "M"),
      module = str_replace(module, ".DEGs", "")
    )
  
  # modules enrichis
  enriched_mod   <- df$module[df$overlap.pval < 0.05 & df$overlap.OR > 1]
  OR_neg <- df$module[df$overlap.pval < 0.05 & df$overlap.OR < 1]
  
  # classification des nodes
  nodes$status <- "NS"
  nodes$status[nodes$name %in% enriched_mod] <- "enriched in DEGs"
 
  nodes$label <- ifelse(nodes$status == "NS", NA, nodes$name)
  
  flareGraph <- graph_from_data_frame(edges, vertices = nodes)
  
  title <- paste0("DEG enriched modules - ", reg, " - ", tp)
  
  color_scale <- scale_color_manual(
    name = "Module enrichment",
    values = c(
      "enriched in DEGs" = "red", 
      "NS" = NA
    )
  )
  
  #### CIRCULAR NETWORK ####
  
  p <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = TRUE) + 
    geom_edge_diagonal(aes(alpha = after_stat(index))) +
    coord_fixed() + 
    scale_edge_alpha('Direction', guide = 'none') +
    geom_node_point(aes(color = status), size = 2) +
    geom_node_text(aes(label = label),
                   repel = TRUE,
                   na.rm = TRUE,
                   size = 2.5,
                   max.overlaps = 100,
                   color = "grey20") +
    color_scale + 
    ggforce::theme_no_axes() +
    ggtitle(title)
  
  ggsave(p, filename = paste0(plots.path,"/network_DEG_",reg,"_",tp,".png"),
    width = 8, height = 8,
    create.dir = TRUE
  )
  
  
  #### DENDROGRAM VERSION ####
  
  p_dendro <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = FALSE) + 
    geom_edge_diagonal(aes(alpha = after_stat(index))) +
    scale_edge_alpha('Direction', guide = 'none') +
    geom_node_point(aes(color = status), size = 2) +
    geom_node_text(aes(label = label),
                   repel = TRUE,
                   na.rm = TRUE,
                   size = 2.5,
                   max.overlaps = 100,
                   color = "grey20") +
    color_scale  +
    ggforce::theme_no_axes() +
    ggtitle(title)
  
  ggsave(p_dendro,filename = paste0(plots.path,"/network_DEG_dendrogram_",reg,"_", tp, ".png"),
    width = 8,height = 8, 
    create.dir = TRUE
  )
}

# Network enrichis en DEGs + cell types 
setwd(project_dir)
deg_cell_types_enrichment.path <- paste0("data/4__MEGENA/MEGENA_without_RIN_correction/enrichment_cell_types_DEGs/", reg, "_enrichment_cell_types_deg.Rdata")
load(deg_cell_types_enrichment.path)
modtable.path <- paste0("data/4__MEGENA/MEGENA_without_RIN_correction/modtable_", reg, ".Rdata")
load(modtable.path)
plot.path <- paste0("graphs_results/4__MEGENA/MEGENA_without_RIN_correction/", reg, "/cell_types_DEG_enrichment_network_plot/")
dir.create(plot.path)

tplist <- c("tp1","tp2","tp3")

for(tp in tplist){
  
  print(paste("Processing", tp))
  
  init_network()
  
  # récupérer table enrichissement
  df <- get(paste0("module_enriched_cell_deg_",reg))
  
  df <- df %>%
    mutate(
      module = str_replace(module, "c1_", "M"), 
      cell_type = str_replace(cell_type, "_more_than_5_fpkm", "")
    ) %>%
    filter(timepoint == tp)
  
  
  # classification des nodes
  nodes$status <- "NS"
  for (type in unique(df$cell_type)) {
    df_type <- df %>% filter(cell_type==type)
    
    nodes$status[nodes$name %in% df_type$module] <- paste0("DEGs ", tp, " and ", type)
  }
  
  nodes$label <- ifelse(nodes$status == "NS", NA, nodes$name)
  
  flareGraph <- graph_from_data_frame(edges, vertices = nodes)
  
  title <- paste0("DEG and cell type enriched modules - ", reg, " - ", tp)
  
  color_scale <- scale_color_manual(
    name = "Module enrichment",
    values = setNames(
      c("yellow3", "red", "pink", "purple", "blue", NA),
      c(
        paste0("DEGs ", tp, " and Astrocytes"),
        paste0("DEGs ", tp, " and Endothelial.Cells"),
        paste0("DEGs ", tp, " and Microglia"),
        paste0("DEGs ", tp, " and Myelinating.Oligodendrocytes"),
        paste0("DEGs ", tp, " and Neuron"),
        "NS"
      )
    ),
    na.value = NA
  )
  
  #### CIRCULAR NETWORK ####
  
  p <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = TRUE) + 
    geom_edge_diagonal(aes(alpha = after_stat(index))) +
    coord_fixed() + 
    scale_edge_alpha('Direction', guide = 'none') +
    geom_node_point(aes(color = status), size = 2) +
    geom_node_text(aes(label = label),
                   repel = TRUE,
                   na.rm = TRUE,
                   size = 2.5,
                   max.overlaps = 100,
                   color = "grey20") +
    color_scale + 
    ggforce::theme_no_axes() +
    ggtitle(title)
  
  ggsave(
    p, filename = paste0(plot.path, "/network_DEG_cell_type_", reg, "_", tp, ".png"),
    width = 8, height = 8, 
    create.dir = TRUE
  )
  
  
  #### DENDROGRAM VERSION ####
  
  p_dendro <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = FALSE) + 
    geom_edge_diagonal(aes(alpha = after_stat(index))) +
    scale_edge_alpha('Direction', guide = 'none') +
    geom_node_point(aes(color = status), size = 2) +
    geom_node_text(aes(label = label),
                   repel = TRUE,
                   na.rm = TRUE,
                   size = 2.5,
                   max.overlaps = 100,
                   color = "grey20") +
    color_scale  +
    ggforce::theme_no_axes() +
    ggtitle(title)
  
  ggsave(
    p_dendro,
    filename = paste0(plot.path, "/network_DEG_cell_type_dendrogram_", reg, "_", tp, ".png"),
    width = 8, height = 8,
    create.dir = TRUE
  )
}




## plot ENORA
#init_network()
## patterns to color
#mod_enriched <- c("M2", "M6", "M76", "M368", "M366", "M786")
#
## Assign colors to nodes based on membership in lists
#nodes$color <- NA
#nodes$color[nodes$name %in% mod_enriched] <-"red"
#
#nodes$label <- ifelse(is.na(nodes$color), NA, nodes$name)
#nodes$signif <- ifelse(nodes$color)
#
## Create the graph object
#flareGraph <- graph_from_data_frame(edges, vertices = nodes)
#
## Plot the graph using ggraph with circular layout
#p <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = TRUE) + 
#  geom_edge_diagonal(aes(alpha = after_stat(index))) +
#  coord_fixed() + 
#  scale_edge_alpha('Direction', guide = 'none') +
#  geom_node_point(aes(color = color), size = 2) +
#  geom_node_text(aes(label = label), repel = TRUE, na.rm = TRUE, size = 2.5, max.overlaps=100, color="grey20") +  # Add labels for colored nodes
#  scale_color_identity() +
#  ggforce::theme_no_axes() +
#  ggtitle(paste0("Network graph of DRN modules"))
#p
#ggsave(p, filename=paste0("modules_", reg, "_", "Enora.png"))
#
#p_dendrogram <- ggraph(flareGraph, 'igraph', algorithm = 'tree', circular = FALSE) + 
#  geom_edge_diagonal(aes(alpha = after_stat(index))) +
#  scale_edge_alpha('Direction', guide = 'none') +
#  geom_node_point(aes(color = color), size = 2) +
#  geom_node_text(aes(label = label), repel = TRUE, na.rm = TRUE, size = 2.5, max.overlaps=100, color="grey20") +  # Add labels for colored nodes
#  scale_color_identity() +
#  ggforce::theme_no_axes() +
#  ggtitle(paste0("Network graph of DRN modules"))
#ggsave(p_dendrogram, filename=paste0("modules_", reg, "_", "_d_Enora", ".png"))
#