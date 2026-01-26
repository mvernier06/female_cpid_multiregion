library(tidyverse)
library(ggplot2)

rm(list=ls())

#### PATHS ####
annotated_counts.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts_filtered.rds"
plot.path <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/graphs_results/2__differential_expression_analysis/"
dir.create(plot.path)
setwd(plot.path)

#### apeglm no outliers ####
de <- readRDS(annotated_counts.path)

degs_acc <- list("tp1" = de$MGI.symbol[which(de$ACC_pval_tp1<0.05)],
                 "tp2" = de$MGI.symbol[which(de$ACC_pval_tp2<0.05)],
                 "tp3" = de$MGI.symbol[which(de$ACC_pval_tp3<0.05)])
length(degs_acc$tp1)
length(degs_acc$tp2)
length(degs_acc$tp3)


degs_hb <- list("tp1" = de$MGI.symbol[which(de$Hb_pval_tp1<0.05)],
                "tp2" = de$MGI.symbol[which(de$Hb_pval_tp2<0.05)],
                "tp3" = de$MGI.symbol[which(de$Hb_pval_tp3<0.05)])
length(degs_hb$tp1)
length(degs_hb$tp2)
length(degs_hb$tp3)

degs_ins <- list("tp1" = de$MGI.symbol[which(de$Ins_pval_tp1<0.05)],
                 "tp2" = de$MGI.symbol[which(de$Ins_pval_tp2<0.05)],
                 "tp3" = de$MGI.symbol[which(de$Ins_pval_tp3<0.05)])
length(degs_ins$tp1)
length(degs_ins$tp2)
length(degs_ins$tp3)

degs_nac <- list("tp1" = de$MGI.symbol[which(de$Nac_pval_tp1<0.05)],
                 "tp2" = de$MGI.symbol[which(de$Nac_pval_tp2<0.05)],
                 "tp3" = de$MGI.symbol[which(de$Nac_pval_tp3<0.05)])
length(degs_nac$tp1)
length(degs_nac$tp2)
length(degs_nac$tp3)


#### deg counts plot ####
df <- data.frame("TP1" = c(length(degs_acc$tp1),
                           length(degs_hb$tp1),
                           length(degs_ins$tp1),
                           length(degs_nac$tp1)
),
"TP2" = c(length(degs_acc$tp2),
          length(degs_hb$tp2),
          length(degs_ins$tp2),
          length(degs_nac$tp2)
),
"TP3" = c(length(degs_acc$tp3),
          length(degs_hb$tp3),
          length(degs_ins$tp3),
          length(degs_nac$tp3)
)
)
df$reg <- c("ACC", "Hb", "Ins", "Nac")
df

## DEGs number plot
df2 <- df %>% pivot_longer(cols=!reg, names_to = "tp", values_to = "number")
df2
plt <- ggplot(df2, aes(x=tp, y=number, color=reg)) +
  geom_point(aes(pch=reg)) +
  geom_line(aes(group=reg)) +
  labs(title="DEG kinetics per TP within regions",
       y = "Number of DEG",
       color = "Region", shape = "Region") +
  theme_bw() +
  theme(axis.title.x=element_blank())
plt
ggsave("deg_counts.png", plot=plt, width=2000, height=1500, units="px", scale=1)



## plot ENORA
#cbPalette <- c("grey30", "red", "grey30", "grey30", "grey30", "grey30")
#plt2 <- ggplot(df2, aes(x=tp, y=number, color=reg)) +
#  geom_point(aes(pch=reg)) +
#  geom_line(aes(group=reg)) +
#  scale_colour_manual(values=cbPalette) +
#  labs(title="DEGs kinetics per TP within regions",
#       y="Number of DEGs", x="")
#plt2
#ggsave("deg_counts_Enora.png", plot=plt2, width=1900, height=1200, units="px", scale=1)