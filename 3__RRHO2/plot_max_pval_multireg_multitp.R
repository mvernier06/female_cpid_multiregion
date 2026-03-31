library(tidyverse)
rm(list=ls())

project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)

load("data/3__RRHO2/rrho_obj_multireg.Rdata")
load("data/3__RRHO2/rrho_obj_multitp.Rdata")


# max(RRHO_ACC_1vs2$hypermat, na.rm = TRUE)
envlist <- ls()
rrho_multireg <- envlist[grepl(paste0("rrho"), envlist)]
rrho_multitp <- envlist[grepl(paste0("RRHO_"), envlist)]
# inverse rrho_multitp et multireg


#### Function to get all 4 quadrant values per rrho comparison ####

# functions to get limits of matrices
nrow0 <- function(x) dim(x)[1]
ncol0 <- function(x) dim(x)[2]

# multireg #
regionList <- c("ACC", "Hb", "Ins", "Nac")

for(reg in regionList){
  print(reg)
  
  # 1v2 #
  rrho <- get(paste0("RRHO_", reg, "_1vs2"))
  mat <- rrho$hypermat
  dim(mat)
  max_row <- mat %>% nrow0
  max_col <- mat %>% ncol0
  
  mat_row <- mat[0:1,]
  mat_col <- mat[,0:1]
  
  na_rows <- which(is.na(mat[0:1,]), arr.ind=TRUE)
  na_cols <- which(is.na(mat[,0:1]), arr.ind=TRUE)
  
  pval <-  c(max(mat[0:min(na_rows), 0:min(na_cols)], na.rm=TRUE),
             max(mat[0:min(na_rows), max(na_cols):max_col], na.rm=TRUE),
             max(mat[max(na_rows):max_row, 0:min(na_cols)], na.rm=TRUE),
             max(mat[max(na_rows):max_row, max(na_cols):max_col], na.rm=TRUE))
  
  assign(paste0("pval_", reg, "_1v2"), pval, envir = .GlobalEnv)
  
  # 1v3 #
  rrho <- get(paste0("RRHO_", reg, "_1vs3"))
  mat <- rrho$hypermat
  dim(mat)
  max_row <- mat %>% nrow0
  max_col <- mat %>% ncol0
  
  mat_row <- mat[0:1,]
  mat_col <- mat[,0:1]
  
  na_rows <- which(is.na(mat[0:1,]), arr.ind=TRUE)
  na_cols <- which(is.na(mat[,0:1]), arr.ind=TRUE)
  
  pval <-  c(max(mat[0:min(na_rows), 0:min(na_cols)], na.rm=TRUE),
             max(mat[0:min(na_rows), max(na_cols):max_col], na.rm=TRUE),
             max(mat[max(na_rows):max_row, 0:min(na_cols)], na.rm=TRUE),
             max(mat[max(na_rows):max_row, max(na_cols):max_col], na.rm=TRUE))
  
  assign(paste0("pval_", reg, "_1v3"), pval, envir = .GlobalEnv)
  
  # 2v3 #
  rrho <- get(paste0("RRHO_", reg, "_2vs3"))
  mat <- rrho$hypermat
  dim(mat)
  max_row <- mat %>% nrow0
  max_col <- mat %>% ncol0
  
  mat_row <- mat[0:1,]
  mat_col <- mat[,0:1]
  
  na_rows <- which(is.na(mat[0:1,]), arr.ind=TRUE)
  na_cols <- which(is.na(mat[,0:1]), arr.ind=TRUE)
  
  pval <-  c(max(mat[0:min(na_rows), 0:min(na_cols)], na.rm=TRUE),
             max(mat[0:min(na_rows), max(na_cols):max_col], na.rm=TRUE),
             max(mat[max(na_rows):max_row, 0:min(na_cols)], na.rm=TRUE),
             max(mat[max(na_rows):max_row, max(na_cols):max_col], na.rm=TRUE))
  
  assign(paste0("pval_", reg, "_2v3"), pval, envir = .GlobalEnv)
}

pval_multitp <- data.frame("ACC" = c(pval_ACC_1v2,
                                     pval_ACC_1v3,
                                     pval_ACC_2v3),
                           "Hb" = c(pval_Hb_1v2,
                                    pval_Hb_1v3,
                                    pval_Hb_2v3),
                           "Ins" = c(pval_Ins_1v2,
                                     pval_Ins_1v3,
                                     pval_Ins_2v3),
                           "Nac" = c(pval_Nac_1v2,
                                     pval_Nac_1v3,
                                     pval_Nac_2v3))
pval_multitp

# multireg #
tpList <- c(1, 2, 3)

list1 <- list2 <- regionList

# homemade loop that iterates each combination of regions
for(tp in tpList){
  print(paste0("Searching RRHO2 pval of the TP", tp))
  i=2
  for(reg in list1[1:length(list1)-1]){
    for(reg2 in list2[i:length(list2)]){
      rrho <- get(paste0("rrho", reg, "vs", reg2, tp))
      mat <- rrho$hypermat
      dim(mat)
      max_row <- mat %>% nrow0
      max_col <- mat %>% ncol0
      
      mat_row <- mat[0:1,]
      mat_col <- mat[,0:1]
      
      na_rows <- which(is.na(mat[0:1,]), arr.ind=TRUE)
      na_cols <- which(is.na(mat[,0:1]), arr.ind=TRUE)
      
      pval <-  c(max(mat[0:min(na_rows), 0:min(na_cols)], na.rm=TRUE),
                 max(mat[0:min(na_rows), max(na_cols):max_col], na.rm=TRUE),
                 max(mat[max(na_rows):max_row, 0:min(na_cols)], na.rm=TRUE),
                 max(mat[max(na_rows):max_row, max(na_cols):max_col], na.rm=TRUE))
      
      assign(paste0("intratp_pval_", reg, "vs", reg2, tp), pval, envir = .GlobalEnv)
    }
    i=i+1
  }
}


envlist <- ls()

test <- mget(c(envlist[str_detect(envlist, "intratp_pval_.*1$")],
               envlist[str_detect(envlist, "intratp_pval_.*2$")],
               envlist[str_detect(envlist, "intratp_pval_.*3$")]))
test

pval_multireg <- data.frame(unlist(test)) %>% rename("pval" = unlist.test.)
pval_multireg

pval_multitp <- pval_multitp %>% pivot_longer(cols = everything(), names_to = "reg", values_to = "pval")

pval_multireg$type <- "multireg"
pval_multitp$type <- "multitp"

df <- merge(pval_multireg, pval_multitp, by = c("pval", "type"), all=TRUE)
df
df %>% filter(type == "multireg") %>% nrow
df %>% filter(type == "multitp") %>% nrow

library(ggpubr)
ggplot(df, aes(y = pval, x = type)) +
  geom_boxplot(outliers = FALSE)


#### prendre que les pval max ####

# intrareg #
for(reg in regionList){
  print(reg)
  
  # 1v2 #
  rrho <- get(paste0("RRHO_", reg, "_1vs2"))
  mat <- rrho$hypermat
  pval <-  max(mat, na.rm=TRUE)
  assign(paste0("pval_", reg, "_1v2"), pval, envir = .GlobalEnv)
  
  # 1v3 #
  rrho <- get(paste0("RRHO_", reg, "_1vs3"))
  mat <- rrho$hypermat
  pval <-  max(mat, na.rm=TRUE)
  assign(paste0("pval_", reg, "_1v3"), pval, envir = .GlobalEnv)
  
  # 2v3 #
  rrho <- get(paste0("RRHO_", reg, "_2vs3"))
  mat <- rrho$hypermat
  pval <-  max(mat, na.rm=TRUE)
  assign(paste0("pval_", reg, "_2v3"), pval, envir = .GlobalEnv)
}

pval_multitp <- data.frame("ACC" = c(pval_ACC_1v2,
                                     pval_ACC_1v3,
                                     pval_ACC_2v3),
                           "Hb" = c(pval_Hb_1v2,
                                    pval_Hb_1v3,
                                    pval_Hb_2v3),
                           "Ins" = c(pval_Ins_1v2,
                                     pval_Ins_1v3,
                                     pval_Ins_2v3),
                           "Nac" = c(pval_Nac_1v2,
                                     pval_Nac_1v3,
                                     pval_Nac_2v3))
pval_multitp

# multireg #
for(tp in tpList){
  print(paste0("Searching RRHO2 pval of the TP", tp))
  i=2
  for(reg in list1[1:length(list1)-1]){
    for(reg2 in list2[i:length(list2)]){
      rrho <- get(paste0("rrho", reg, "vs", reg2, tp))
      mat <- rrho$hypermat
      pval <-  max(mat, na.rm=TRUE)
      
      assign(paste0("multireg_pval_", reg, "vs", reg2, tp), pval, envir = .GlobalEnv)
    }
    i=i+1
  }
}

envlist <- ls()
test <- mget(c(envlist[str_detect(envlist,  "multireg_pval_.*1$")],
               envlist[str_detect(envlist, "multireg_pval_.*2$")],
               envlist[str_detect(envlist, "multireg_pval_.*3$")]))
test

pval_multireg <- data.frame(unlist(test)) %>% rename("pval" = unlist.test.)
pval_multireg

pval_multitp <- pval_multitp %>% pivot_longer(cols = everything(), names_to = "reg", values_to = "pval")
pval_multireg$type <- "multireg" # RENAMING INTRATIP AS MULTIREG
pval_multitp$type <- "multitp" # RENAMING INTRAREG AS MULTITP

df2 <- merge(pval_multireg, pval_multitp, by = c("pval", "type"), all=TRUE)
df2
df2 %>% filter(type == "multireg") %>% nrow
df2 %>% filter(type == "multitp") %>% nrow
df2
df2 <- df2 %>%
  dplyr::select(pval, type) %>%
  group_by(type) %>%
  mutate(n = n())
df2

plot.path <- "graphs_results/3__RRHO2/pval_multireg_vs_multitp"
dir.create(plot.path)
setwd(plot.path)
library(ggpubr)


my_comparisons <- list(c("multireg", "multitp"))

# graph rrho MAX pval (Fig. 2B)
ggplot(df2, aes(y = pval, x = type, fill = type)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  geom_label(aes(y = -15, label = n), show.legend = FALSE) +
  labs(title = "RRHO comparison MAX pval",
       y = "-Log10( P-value )",
       x = element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label.y = 250) +
  scale_fill_manual(values = c("pink", "lightgreen")) +
  theme_classic()+
  guides(fill = guide_legend(title = NULL)) 
ggsave(plot = last_plot(), "rrho_pval_max_points.png", width=1250, height=1250, units="px", scale=1)
getwd()

df <- df %>%
  dplyr::select(pval, type) %>%
  group_by(type) %>%
  mutate(n = n())
df

# graph rrho ALL pval
ggplot(df, aes(y = pval, x = type)) +
  geom_boxplot(outliers = FALSE) +
  geom_label(aes(y = -15, label = n)) +
  labs(title = "RRHO comparison ALL pval",
       y = "-Log10( P-value )",
       x = element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label.y = 100) +
  theme_bw()
# ggsave(plot = last_plot(), "rrho_pval_all.png")


# #### COMPARAISON DE 25% DES PVALS MAX ####
# 
# # intrareg #
# for(reg in regionList){
#   print(reg)
#   
#   # 1v2 #
#   rrho <- get(paste0("RRHO_", reg, "_1vs2"))
#   mat <- rrho$hypermat
#   matlist <- unlist(na.omit(lapply(list(mat),sort,decreasing=TRUE)) )
#   length(matlist)
#   # on choisit de prendre le top 25% des pval
#   top_pval <- head(matlist, length(matlist)*0.25)
#   assign(paste0("top_pval_", reg, "_1v2"), top_pval, envir = .GlobalEnv)
#   
#   # 1v3 #
#   rrho <- get(paste0("RRHO_", reg, "_1vs3"))
#   mat <- rrho$hypermat
#   matlist <- unlist(na.omit(lapply(list(mat),sort,decreasing=TRUE)) )
#   top_pval <- head(matlist, length(matlist)*0.25)
#   assign(paste0("top_pval_", reg, "_1v3"), top_pval, envir = .GlobalEnv)
#   
#   # 2v3 #
#   rrho <- get(paste0("RRHO_", reg, "_2vs3"))
#   mat <- rrho$hypermat
#   matlist <- unlist(na.omit(lapply(list(mat),sort,decreasing=TRUE)) )
#   top_pval <- head(matlist, length(matlist)*0.25)
#   assign(paste0("top_pval_", reg, "_2v3"), top_pval, envir = .GlobalEnv)
# }
# 
# # multireg #
# for(tp in tpList){
#   print(paste0("Searching RRHO2 pval of the TP", tp))
#   i=2
#   for(reg in list1[1:length(list1)-1]){
#     for(reg2 in list2[i:length(list2)]){
#       rrho <- get(paste0("rrho", reg, "vs", reg2, tp))
#       mat <- rrho$hypermat
#       matlist <- unlist(na.omit(lapply(list(mat),sort,decreasing=TRUE)) )
#       top_pval <- head(matlist, length(matlist)*0.25)
#       
#       assign(paste0("multireg_top_pval_", reg, "vs", reg2, tp), top_pval, envir = .GlobalEnv)
#     }
#     i=i+1
#   }
# }
# 
# envlist <- ls()
# test <- mget(c(envlist[str_detect(envlist,  "multireg_top_pval_.*1$")],
#                envlist[str_detect(envlist, "multireg_top_pval_.*2$")],
#                envlist[str_detect(envlist, "multireg_top_pval_.*3$")]))
# test
# pval_multireg <- data.frame(unlist(test)) %>% rename("pval" = unlist.test.)
# pval_multireg
# 
# pval_multitp <- list("BLA" = c(top_pval_BLA_1v2,
#                                top_pval_BLA_1v3,
#                                top_pval_BLA_2v3),
#                      "DRN" = c(top_pval_DRN_1v2,
#                                top_pval_DRN_1v3,
#                                top_pval_DRN_2v3),
#                      "Hb" = c(top_pval_Hb_1v2,
#                               top_pval_Hb_1v3,
#                               top_pval_Hb_2v3),
#                      "Ins" = c(top_pval_Ins_1v2,
#                                top_pval_Ins_1v3,
#                                top_pval_Ins_2v3),
#                      "NAc" = c(top_pval_NAc_1v2,
#                                top_pval_NAc_1v3,
#                                top_pval_NAc_2v3),
#                      "VTA" = c(top_pval_VTA_1v2,
#                                top_pval_VTA_1v3,
#                                top_pval_VTA_2v3))
# pval_multitp <- as.data.frame(unlist(pval_multitp), row.names = c(1:length(unlist(pval_multitp))))
# pval_multitp
# colnames(pval_multitp) <- "pval"
# 
# pval_multireg$type <- "multireg"
# pval_multitp$type <- "multitp"
# 
# df2 <- merge(pval_multireg, pval_multitp, by = c("pval", "type"), all=TRUE)
# df2
# df2 %>% filter(type == "multireg") %>% nrow
# df2 %>% filter(type == "multitp") %>% nrow
# df2 %>% 
#   group_by(type) %>%
#   get_summary_stats(pval, type = "mean_sd")
# df2 %>% 
#   group_by(type) %>% 
#   compare_means(pval ~ type, .)
# 
# plot.path <- "~/CPID_multiregion/graphs_results/4__RRHO2/pval_intrareg_vs_intratp"
# # dir.create(plot.path)
# setwd(plot.path)
# library(ggpubr)
# df2
# plot.path <- "~/CPID_multiregion/graphs_results/4__RRHO2/pval_intrareg_vs_intratp"
# setwd(plot.path)
# ggplot(df2, aes(y = pval, x = type)) +
#   geom_boxplot(outliers = FALSE) +
#   labs(title = "RRHO comparison of 25% of MAX pval",
#        y = "-log10( pval )",
#        x = element_blank()) +
#   stat_compare_means(comparisons = my_comparisons, label.y = 210) +
#   theme_bw()
# # ggsave(plot = last_plot(), "rrho_pval_25%_max.png")
# 
# df2 %>% filter(type == "multireg") %>% nrow # 141960
# df2 %>% filter(type == "multitp") %>% nrow # 62661