## Attention, faut restart R avant de charger les librairies ##

library(tidyverse)
library(ggplot2)
library(ggalluvial)

rm(list=ls())

#### PATHS ####
deglist.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/deglist.Rdata"
tempcountlist.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/tempcountlist.Rdata"
plot.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/alluvial"
output.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/"

dir.create(plot.path)
load(deglist.path)
load(tempcountlist.path)
setwd(output.path)


#### GET ALL PATTERNS OF GENE EXPRESSION ####
get_alluvial_patterns <- function(regionlist){
  for(reg in regionlist){
    print(paste("Alluvial patterns", reg))
    
    deg_tp1 <- get(paste0("deg_", reg, "_tp1"))
    deg_tp2 <- get(paste0("deg_", reg, "_tp2"))
    deg_tp3 <- get(paste0("deg_", reg, "_tp3"))
    
    deg_genes <- unique(c(deg_tp1$label, deg_tp2$label, deg_tp3$label))
    
    df1 <- get(paste0("temp_counts_", reg, "_tp1"))
    df1 <- df1[which(df1$MGI.symbol %in% deg_genes),] %>%
      rename_at(vars(-MGI.symbol), ~ paste0(., '_tp1'))
    df2 <- get(paste0("temp_counts_", reg, "_tp2"))
    df2 <- df2[which(df2$MGI.symbol %in% deg_genes),] %>%
      rename_at(vars(-MGI.symbol), ~ paste0(., '_tp2'))
    df3 <- get(paste0("temp_counts_", reg, "_tp3"))
    df3 <- df3[which(df3$MGI.symbol %in% deg_genes),] %>%
      rename_at(vars(-MGI.symbol), ~ paste0(., '_tp3'))
    
    df_list <- list(df1, df2, df3)
    
    #merge all data frames in list
    df <- df_list %>% reduce(full_join, by="MGI.symbol")
    
    df_new <- df
    # key
    df_new <- df_new %>% 
      relocate(MGI.symbol, 
               log2fc_tp1, log2fc_tp2, log2fc_tp3,
               diffexpressed_tp1, diffexpressed_tp2, diffexpressed_tp3) %>%
      mutate(tp1=NULL, tp2=NULL, tp3=NULL, key=NULL) %>%
      add_column(diffexpressed_alltp=paste(.$diffexpressed_tp1, .$diffexpressed_tp2, .$diffexpressed_tp3, sep="_"))
    
    df_new$diffexpressed_alltp <- as.factor(df_new$diffexpressed_alltp)
    df_new <- df_new %>%
      arrange(desc(diffexpressed_alltp), desc(diffexpressed_tp1), desc(diffexpressed_tp2), desc(diffexpressed_tp3))
    assign(paste0("df_new_", reg), df_new, envir = .GlobalEnv)
    save(df_new, file=paste0("alluvial_patterns_", reg, ".Rdata"))
  }
}

regionlist <- c("ACC", "Hb", "Ins", "Nac")
get_alluvial_patterns(regionlist)



#### COLOR PALETTE ####
# get a list of all patterns
all_patterns_lst <- unique(c(df_new_ACC %>% .$diffexpressed_alltp %>% unique,
                             df_new_Hb %>% .$diffexpressed_alltp %>% unique,
                             df_new_Ins %>% .$diffexpressed_alltp %>% unique,
                             df_new_Nac %>% .$diffexpressed_alltp %>% unique))
all_patterns_lst <- all_patterns_lst[order(as.character(all_patterns_lst))]
write_rds(all_patterns_lst, "all_patterns_lst.rds")

# same color palette for all graphs
get_color_pal <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- get_color_pal(length(all_patterns_lst))
names(pal) <- all_patterns_lst
pal[14] <- "palegreen3" # change ns_ns_UP color because to close from ns_ns_DOWN for the BLA
write_rds(pal, "alluvial_palette.rds")




#### ALLUVIAL PLOTS ####
alluvial_plots <- function(regionlist){
  for(reg in regionlist){
    print(paste("Printing alluvial ", reg))
    
    df <- get(paste0("df_new_", reg)) %>% 
      dplyr::group_by(diffexpressed_tp1, diffexpressed_tp2, diffexpressed_tp3,
                      diffexpressed_alltp) %>% 
      summarise(n = n()) 
    
    # Define the factors
    df$diffexpressed_tp1 <- factor(df$diffexpressed_tp1, levels = c("UP", "ns", "DOWN"))
    df$diffexpressed_tp2 <- factor(df$diffexpressed_tp2, levels = c("UP", "ns", "DOWN"))
    df$diffexpressed_tp3 <- factor(df$diffexpressed_tp3, levels = c("UP", "ns", "DOWN"))
    df$diffexpressed_alltp <- factor(df$diffexpressed_alltp)
    
    print(df$diffexpressed_alltp)
    pal_reg <- pal[names(pal) %in% df$diffexpressed_alltp]
    print(length(pal_reg))
    
    df <- df %>% arrange(desc(diffexpressed_alltp))
    
    top_patterns <- df %>% arrange(desc(n))
    setwd(output.path)
    write_rds(top_patterns, paste0("top_patterns_", reg, ".rds"))
    
    p_alluvial <- ggplot(df,
                         aes(y = n,
                             axis1 = diffexpressed_tp1, 
                             axis2 = diffexpressed_tp2, 
                             axis3 = diffexpressed_tp3)) +
      geom_alluvium(aes(fill = diffexpressed_alltp)) +
      guides(fill = "none") +
      geom_stratum(alpha = .25, width = 1/7, reverse = TRUE) +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)),
                reverse = TRUE) +
      scale_x_continuous(breaks = 1:3, labels = c("TP1", "TP2", "TP3")) +
      labs(y="Genes") +
      ggtitle(paste("Pattern of gene expression within the", reg)) +
      scale_fill_manual(values=unname(pal_reg)) +
      theme_bw()
    setwd(plot.path)
    ggsave(paste0("alluvial_", reg, ".png"), plot=p_alluvial, width=1900, height=1200, units="px", scale=2)  
  }
}

alluvial_plots(regionlist)



# # plot ENORA
# reg <- "ACC"
# 
# df <- get(paste0("df_new_", reg)) %>% 
#   dplyr::group_by(diffexpressed_tp1, diffexpressed_tp2, diffexpressed_tp3,
#                   diffexpressed_alltp) %>% 
#   summarise(n = n()) 
# 
# # Define the factors
# df$diffexpressed_tp1 <- factor(df$diffexpressed_tp1, levels = c("UP", "ns", "DOWN"))
# df$diffexpressed_tp2 <- factor(df$diffexpressed_tp2, levels = c("UP", "ns", "DOWN"))
# df$diffexpressed_tp3 <- factor(df$diffexpressed_tp3, levels = c("UP", "ns", "DOWN"))
# df$diffexpressed_alltp <- factor(df$diffexpressed_alltp)
# 
# print(df$diffexpressed_alltp)
# pal_reg <- pal[names(pal) %in% df$diffexpressed_alltp]
# print(length(pal_reg))
# 
# df <- df %>% arrange(desc(diffexpressed_alltp))
# 
# top_patterns <- df %>% arrange(desc(n))
# setwd(output.path)
# 
# p_alluvial <- ggplot(df,
#                      aes(y = n,
#                          axis1 = diffexpressed_tp1, 
#                          axis2 = diffexpressed_tp2, 
#                          axis3 = diffexpressed_tp3)) +
#   geom_alluvium(aes(fill = diffexpressed_alltp)) +
#   guides(fill = "none") +
#   geom_stratum(alpha = .25, width = 1/7, reverse = TRUE) +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)),
#             reverse = TRUE) +
#   scale_x_continuous(breaks = 1:3, labels = c("TP1", "TP2", "TP3")) +
#   labs(y="Number of genes") +
#   ggtitle(paste("Pattern of gene expression within the", reg)) +
#   scale_fill_manual(values=unname(pal_reg)) +
#   theme_bw()
# p_alluvial
# setwd(plot.path)
# ggsave(paste0("alluvial_", reg, "_Enora.png"), plot=p_alluvial, width=1900, height=1200, units="px", scale=1)  