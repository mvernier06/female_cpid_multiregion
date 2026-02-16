library(tidyverse)
library(ggpubr)
library(rstatix)

rm(list=ls())


#### PATHS ####
annotated_counts.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts_filtered.rds"
plot.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/lfc_kinetics_multireg"
palette.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/alluvial_palette.rds"

dir.create(plot.path)
setwd(plot.path)
annotated_counts <- read_rds(annotated_counts.path)
pal <- read_rds(palette.path)

#### plot kinetics of lfc for each reg ####
regionlist <- c("ACC", "Hb", "Ins", "Nac")

plotLFCmultireg <- function(reg){
  cat("\n#### LFC Boxplots of", reg, "####", fill=TRUE)
  
  alluvial_patterns.path <- paste0("/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/alluvial_patterns_", reg, ".Rdata")
  top_patterns.path <- paste0("/home/marinevernier/Documents/projets/female_cpid_multiregion/data/2__differential_expression_analysis/top_patterns_", reg, ".rds")
  load(alluvial_patterns.path) # df_new
  top_patterns <- read_rds(top_patterns.path)
  
  pattern1 <- as.character(top_patterns$diffexpressed_alltp[1])
  pattern1color <- pal[pattern1]
  pattern2 <- as.character(top_patterns$diffexpressed_alltp[2])
  pattern2color <- pal[pattern2]
  pattern3 <- as.character(top_patterns$diffexpressed_alltp[3])
  pattern3color <- pal[pattern3]
  pattern4 <- as.character(top_patterns$diffexpressed_alltp[4])
  pattern4color <- pal[pattern4]
  
  patterns_ordered <- df_new %>%
    dplyr::select(diffexpressed_alltp) %>%
    group_by(diffexpressed_alltp) %>%
    mutate(n = n()) %>%
    unique %>%
    arrange(-n)
  
  IDpattern <- 1
  
  for(pattern in head(patterns_ordered$diffexpressed_alltp, 4)){
    cat("\n", pattern, fill=TRUE)
    
    pattern_color <- get(paste("pattern", IDpattern, "color", sep = ""))
    
    df_test <- annotated_counts %>% 
      select(MGI.symbol,contains(paste0(reg, "_log2fc"))) %>%
      rename(c("tp1" = paste0(reg, "_log2fc_tp1"),
               "tp2" = paste0(reg, "_log2fc_tp2"), 
               "tp3" = paste0(reg, "_log2fc_tp3"))) %>%
      pivot_longer(cols=!c(MGI.symbol), names_to = "timepoint", values_to = "lfc") %>%
      na.omit()
    df_test$diffexpressed_alltp <- "all_genes"
    df_test
    df <- df_new %>% 
      filter(diffexpressed_alltp == pattern) %>%
      select(MGI.symbol, contains("log2fc"), diffexpressed_alltp) %>%
      rename(c("tp1" = "log2fc_tp1", "tp2" =  "log2fc_tp2" , "tp3" = "log2fc_tp3" )) %>%
      pivot_longer(cols=!c(MGI.symbol, diffexpressed_alltp), names_to = "timepoint", values_to = "lfc")
    
    df2 <- merge(df_test, df, by=c("MGI.symbol", "timepoint", "lfc", "diffexpressed_alltp"), all=TRUE)
    
    # print statistics and ANOVA results
    df2 %>%
      group_by(timepoint, diffexpressed_alltp) %>%
      get_summary_stats(lfc, type = "mean_sd") %>%
      print
    df2 %>% 
      anova_test(lfc ~ timepoint * diffexpressed_alltp) %>%
      print
    
    n <- patterns_ordered$n[which(patterns_ordered$diffexpressed_alltp == pattern)]
    print(paste("Max size of pattern", pattern, ":", n))
    
    # compute upper whiskers of lfc at tp1 for all genes and pattern
    lfc_tp1_allgenes <- df2 %>% filter(timepoint == "tp1")
    lfc_tp1_pattern <- df %>% filter(timepoint == "tp1")
    ylim1 = max(boxplot.stats(lfc_tp1_allgenes$lfc)$stats[5], boxplot.stats(lfc_tp1_pattern$lfc)$stats[5])
    # compute upper whiskers of lfc at tp2 for all genes and pattern
    lfc_tp2_allgenes <- df2 %>% filter(timepoint == "tp2")
    lfc_tp2_pattern <- df %>% filter(timepoint == "tp2")
    ylim2 = max(boxplot.stats(lfc_tp2_allgenes$lfc)$stats[5], boxplot.stats(lfc_tp2_pattern$lfc)$stats[5])
    # compute upper whiskers of lfc at tp3 for all genes and pattern
    lfc_tp3_allgenes <- df2 %>% filter(timepoint == "tp3")
    lfc_tp3_pattern <- df %>% filter(timepoint == "tp3")
    ylim3 = max(boxplot.stats(lfc_tp3_allgenes$lfc)$stats[5], boxplot.stats(lfc_tp3_pattern$lfc)$stats[5])
    
    p <- ggplot(df2, aes(x=timepoint, y=lfc, fill=diffexpressed_alltp)) +
      geom_boxplot(outliers = FALSE) +
      geom_hline(yintercept=0, linetype="dotted") +
      scale_fill_manual(values=c("grey90", alpha(colour = unname(pattern_color), 0.5))) +
      labs(title=paste0(reg, " genes ", pattern, " (", n, ")"),
           fill="pattern") +
      stat_compare_means(method = "anova", label="p.signif",
                         label.y = c(ylim1, ylim2, ylim3), show.legend = F) +
      theme_bw()
    
    assign(paste0("p_lfc_", IDpattern), p, envir = .GlobalEnv)
    IDpattern <- IDpattern + 1
    
    file.name <- paste0("lfc_kinetics_", reg, "_", pattern, ".png")
    ggsave(file.name, plot=p, width=2796, height=1796, units="px")
  }
}
# plotLFCmultireg("ACC")

for(reg in regionlist){
  print(paste0("Printing lfc plots of ", reg))
  plotLFCmultireg(reg)
}
