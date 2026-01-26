library(tidyverse)
library(ggsci)
library(ggpubr)

rm(list=ls())

#### PATHS ####
annotated_counts.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts_filtered.rds"
plot.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/"
setwd(plot.path)

#### apeglm no outliers ####
annotated_counts <- read_rds(annotated_counts.path)

df <- annotated_counts %>% select(MGI.symbol,contains("log2fc"))
df <- df %>% pivot_longer(cols=!MGI.symbol,
                          values_to="fc",
                          names_to="info")
df$stat <- df$info
df <- df %>% separate_wider_delim(info,
                                  "_log2fc_", 
                                  names = c("region", "tp"))
df$fcabs <- abs(df$fc) 
head(df)
colnames(df)

p1 <- ggplot(df, aes(x = region, y = fcabs, color = tp, fill = tp)) +
  geom_boxplot(outliers = FALSE) +
  labs(y="|Log2FC|",
       title="Absolute Fold Change values per TP for each region",
       fill = "Timepoint", color = "Timepoint") +
  coord_cartesian(ylim=c(0,0.8)) + 
  scale_color_manual(values = c("#8594B6", "#8A4F21", "#465553")) +
  scale_fill_manual(values = c("#F3F6FB", "#EDD9BA", "#D9E8E5")) +
  theme_bw() +
  theme(axis.title.x = element_blank())
p1
ggsave(filename = "lfc_multireg.png", plot=p1, width = 2000, height = 1500, dpi = 300, units = "px")
