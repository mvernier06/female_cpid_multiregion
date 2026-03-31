library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

rm(list=ls())

project_path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/"
setwd(project_path)
load("data/3__RRHO2/rrho_obj_multireg.Rdata")
load("data/3__RRHO2/rrho_obj_multitp.Rdata")

plot.path <- "graphs_results/3__RRHO2/pval_multireg_vs_multitp/"


find_rrho_split <- function(mat) {
  row_split <- which(is.na(mat[, 1]))[1]
  col_split <- which(is.na(mat[1, ]))[1]
  list(row = row_split, col = col_split)
}

# max pval per quadrant 
extract_rrho_quadrant_max <- function(rrho_obj) {
  
  mat <- rrho_obj$hypermat
  split <- find_rrho_split(mat)
  r <- split$row
  c <- split$col
  
  c(
    up_up   = max(mat[1:(r-1), 1:(c-1)], na.rm = TRUE),
    up_down = max(mat[1:(r-1), c:ncol(mat)], na.rm = TRUE),
    down_up     = max(mat[r:nrow(mat), 1:(c-1)], na.rm = TRUE),
    down_down   = max(mat[r:nrow(mat), c:ncol(mat)], na.rm = TRUE)
  )
}


regionList <- c("ACC", "Hb", "Ins", "Nac")
tp_comparisons <- c("1vs2", "1vs3", "2vs3")

multitp_list <- list()

for (reg in regionList) {
  for (tp in tp_comparisons) {
    
    obj_name <- paste0("RRHO_", reg, "_", tp)
    rrho <- get(obj_name)
    
    pval <- extract_rrho_quadrant_max(rrho)
    
    multitp_list[[obj_name]] <- data.frame(
      pval = pval,
      quadrant = names(pval),
      type = "multitp",
      comparison = obj_name
    )
  }
}

df_multitp <- bind_rows(multitp_list)

tpList <- c(1, 2, 3)


multireg_list <- list()

for (tp in tpList) {
  for (i in seq_len(length(regionList) - 1)) {
    for (j in (i + 1):length(regionList)) {
      
      reg1 <- regionList[i]
      reg2 <- regionList[j]
      
      obj_name <- paste0("rrho", reg1, "vs", reg2, tp)
      rrho <- get(obj_name)
      
      pval <- extract_rrho_quadrant_max(rrho)
      
      multireg_list[[obj_name]] <- data.frame(
        pval = pval,
        quadrant = names(pval),
        type = "multireg",
        comparison = obj_name
      )
    }
  }
}

df_multireg <- bind_rows(multireg_list)

df <- bind_rows(df_multireg, df_multitp)

df$type <- factor(df$type, levels = c("multireg", "multitp"))
df$quadrant <- factor(
  df$quadrant,
  levels = c("up_down", "down_down", "up_up", "down_up")
)




# add n per group
df <- df %>%
  group_by(type, quadrant) %>%
  mutate(n = n()) %>%
  ungroup()
save(df, file="max_pval_per_quadrant.Rdata")

wilcox_res <- df %>%
  group_by(quadrant) %>%
  summarise(
    p_value = wilcox.test(
      pval ~ type,
      exact = FALSE
    )$p.value
  )

wilcox_res
df %>%
  group_by(quadrant, type) %>%
  summarise(
    min    = min(pval, na.rm = TRUE),
    q1     = quantile(pval, 0.25, na.rm = TRUE),
    median = median(pval, na.rm = TRUE),
    mean   = mean(pval, na.rm = TRUE),
    q3     = quantile(pval, 0.75, na.rm = TRUE),
    max    = max(pval, na.rm = TRUE),
    n      = n(),
    .groups = "drop"
  )


wilcox_res <- wilcox_res %>%
  mutate(
    label = paste0("Wilcoxon p = ", signif(p_value, 3))
  )

my_comparisons <- list(c("multireg", "multitp"))


n_labels <- df %>%
  group_by(quadrant, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    x = ifelse(type == "multireg", 1, 2)
  )

ggplot(df, aes(y = pval, x = type)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  geom_label(aes(y = -15, label = n), show.legend = FALSE) +
  labs(title = "RRHO comparison MAX pval",
       y = "-Log10( P-value )",
       x = element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label.y = 90) +
  facet_wrap(~ quadrant, ncol = 2) +
  theme_classic()+
  guides(fill = guide_legend(title = NULL)) 

file_name <- paste0 (plot.path, "plot_max_pval_per_RRHO_quadrant.png")
ggsave(plot = last_plot(), file_name)

# df$quadrant <- factor(
#   df$quadrant,
#   levels = c("up_down", "down_down", "up_up", "down_up")
# )

ggplot(df, aes(x = pval, fill = type)) +
  geom_density(alpha = 0.4, adjust = 1) +
  facet_wrap(~ quadrant, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Max -log10(p-value)",
    y = "Density",
    title = "Distribution des -log10(p-values) max par quadrant"
  )
file_name <- paste0(plot.path,  "density_max_pval_per_quadrant.png")
ggsave(plot = last_plot(), file_name)

ggplot(df, aes(x = pval, fill = type)) +
  geom_histogram(bins = 70,
                 alpha = 0.6,
                 position = "identity") +
  facet_wrap(~ quadrant) +
  theme_classic() +
  labs(
    x = "Max -log10(p-value)",
    y = "Count",
    title = "Distribution des -log10(p-values) max par quadrant"
  )
file_name <- paste0(plot.path, "histogram_max_pval_per_quadrant.png")
ggsave(plot = last_plot(), file_name)


###############################################################################
## Meme chose pour concordant vs discordant ##

df2 <- df %>%
  mutate(
    class = case_when(
      quadrant %in% c("up_up", "down_down") ~ "Concordant",
      quadrant %in% c("up_down", "down_up") ~ "Discordant"
    )
  )
df2
n_labels <- df2 %>%
  group_by(class, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(x = as.numeric(factor(type)))

wilcox_res <- df2 %>%
  group_by(class) %>%
  summarise(
    p_value = wilcox.test(pval ~ type, exact = FALSE)$p.value,
    .groups = "drop"
  )
wilcox_res

plot_labels <- wilcox_res %>%
  mutate(
    label = paste0("Wilcoxon p = ", signif(p_value, 3))
  )
ggplot(df2, aes(y = pval, x = type)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  geom_label(aes(y = -15, label = n), show.legend = FALSE) +
  labs(
    title = "RRHO – MAX -log10(p-value)",
    subtitle = "Concordant vs Discordant ",
    y = "-log10(P-value)",
    x = NULL
  ) +
  stat_compare_means(comparisons = my_comparisons, label.y = 110) +
  facet_wrap(~ class, ncol = 2) +
  theme_classic()+
  guides(fill = guide_legend(title = NULL)) 

file_name <- paste0(plot.path,"max_pval_RRHO_discordants_concordants.png" )
ggsave(last_plot(), file=file_name)

df2 %>%
  group_by(class, type) %>%
  summarise(
    min    = min(pval, na.rm = TRUE),
    q1     = quantile(pval, 0.25, na.rm = TRUE),
    median = median(pval, na.rm = TRUE),
    mean   = mean(pval, na.rm = TRUE),
    q3     = quantile(pval, 0.75, na.rm = TRUE),
    max    = max(pval, na.rm = TRUE),
    n      = n(),
    .groups = "drop"
  )

ggplot(df2, aes(x = pval, fill = type)) +
  geom_density(alpha = 0.4, adjust = 1) +
  facet_wrap(~ class, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Max -log10(p-value)",
    y = "Density",
    title = "Distribution des -log10(p-values) max discordant vs concordant"
  )
file_name <- paste0(plot.path, "density_max_pval_discordant_vs_concordant.png")
ggsave(plot = last_plot(),file_name )

ggplot(df2, aes(x = pval, fill = type)) +
  geom_histogram(bins = 70,
                 alpha = 0.6,
                 position = "identity") +
  facet_wrap(~ class) +
  theme_classic() +
  labs(
    x = "Max -log10(p-value)",
    y = "Count",
    title = "Distribution des -log10(p-values) max discordant vs concordant"
  )
file_name <- paste0(plot.path, "histogram_max_pval_discordant_vs_concordant.png")
ggsave(plot = last_plot(), file_name)
