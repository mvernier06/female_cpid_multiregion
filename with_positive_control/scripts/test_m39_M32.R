library(tidyverse)
library(readODS)
library(tidyr)
library(patchwork)
library(purrr)

rm(list = ls())
getwd()
origine.path <- "/home/marinevernier/Documents/projets/cpid_multiregion/data/2__differential_expression_analysis/annotated_counts.csv"
test.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/with_positive_control/data/m39_M32_final/cpid_multireg_control_counts.txt"

origine <- read_csv(origine.path)
test <- read.table(test.path, header = TRUE,
                   sep = "\t",
                   comment.char = "#",
                   stringsAsFactors = FALSE,
                   check.names = FALSE)

test <- test[, !(colnames(test) %in% c(
  "Chr", "Start", "End", "Strand", "Length"
))]

rownames(test) <- test$Geneid
test <- test[, -1]

colnames(test) <- sub(
  "_R2\\.dedup\\.bam$",
  "",
  basename(colnames(test))
)
colnames(test) <- sub(
  "^([A-Za-z]+)([0-9]+)$",
  "\\1.\\2",
  colnames(test)
)
colnames(test) <- gsub("-", ".", colnames(test))

## Annotation 
annot_table.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/data/counts_m39_M32/annotation_final.csv"

# Annotate counts
test_counts <- test
test_counts <-  rownames_to_column(test_counts)
test_counts$rowname <-  str_replace(test_counts$rowname, "\\..*", "") 

ens2symbol <- read.csv(annot_table.path, sep=";")
ens2symbol <- ens2symbol %>% dplyr::select(!X)

test_counts <- inner_join(test_counts, ens2symbol, by=c("rowname"="Gene.stable.ID"))
test_counts <- test_counts %>% relocate(MGI.symbol, .before=rowname) 
# remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)
test_counts <- test_counts[!duplicated(test_counts$MGI.symbol),]

rownames(test_counts) <- test_counts$MGI.symbol
test_counts$MGI.symbol <- NULL
test_counts$rowname <- NULL

origine <- origine %>%
  column_to_rownames("MGI.symbol")
origine$X <- NULL

test_sample <- colnames(test_counts)
origine <- origine %>%
  select(all_of(test_sample))


diff_counts <- test_counts - origine 


correlations <- sapply(colnames(origine), function(s) {
  cor(origine[, s], test_counts[, s], method = "spearman")
})

correlations



compare_one_sample <- function(sample_name,
                               counts_origine,
                               counts_test) {
  
  df_compare <- full_join(
    counts_origine %>%
      select(MGI.symbol, origine = all_of(sample_name)),
    counts_test %>%
      select(MGI.symbol, test = all_of(sample_name)),
    by = "MGI.symbol"
  ) %>%
    mutate(diff = test - origine)
  
  df_scatter <- df_compare %>%
    filter(!is.na(origine) & !is.na(test))
  
  cor_spearman <- cor(
    df_scatter$origine,
    df_scatter$test,
    method = "spearman"
  )
  
  n_total <- nrow(df_compare)
  
  table_recap <- tibble(
    sample = sample_name,
    categorie = c(
      "diff = 0",
      "diff ≠ 0",
      "diff = NA",
      "test = NA & origine = 0",
      "origine = NA & test = 0"
    ),
    n_genes = c(
      sum(df_compare$diff == 0, na.rm = TRUE),
      sum(df_compare$diff != 0, na.rm = TRUE),
      sum(is.na(df_compare$diff)),
      sum(is.na(df_compare$test) & df_compare$origine == 0, na.rm = TRUE),
      sum(is.na(df_compare$origine) & df_compare$test == 0, na.rm = TRUE)
    )
  ) %>%
    mutate(percent = 100 * n_genes / n_total)
  
  list(
    df_compare = df_compare,
    df_scatter = df_scatter,
    cor_spearman = cor_spearman,
    table_recap = table_recap
  )
}

origine <- rownames_to_column(origine, "MGI.symbol")
test_counts <- rownames_to_column(test_counts, "MGI.symbol")

results <- lapply(test_sample, function(s) {
  compare_one_sample(
    sample_name = s,
    counts_origine = origine,
    counts_test = test_counts
  )
})

names(results) <- test_sample

s <- "Ins.1788"

results[[s]]$df_compare %>%
  filter(diff != 0 | is.na(diff)) %>%
  mutate(diff_bin = cut(diff, breaks = 20)) %>%
  mutate(diff_bin = addNA(diff_bin)) %>%
  ggplot(aes(x = diff_bin)) +
  geom_bar(fill = "steelblue", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(results[[s]]$df_scatter,
       aes(x = origine, y = test)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = paste("Comparaison origine vs test –", s),
    x = "Counts origine",
    y = "Counts test"
  )

table_recap_all <- bind_rows(
  lapply(results, function(x) x$table_recap)
)

table_recap_all

###################################################################################################################
################                   meme chose en ne gardant que les gènes après filtrage ###################


filtered_genes <- read.csv(
  "/home/marinevernier/Documents/projets/cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv",
  stringsAsFactors = FALSE
)

genes_keep <- filtered_genes$MGI.symbol
length(genes_keep)

origine_filt <- origine %>%
  filter(MGI.symbol %in% genes_keep)

test_counts_filt <- test_counts %>%
  filter(MGI.symbol %in% genes_keep)

correlations <- sapply(test_sample, function(s) {
  cor(origine_filt[, s], test_counts_filt[, s], method = "spearman")
})

correlations

results_filt <- lapply(test_sample, function(s) {
  compare_one_sample(
    sample_name = s,
    counts_origine = origine_filt,
    counts_test = test_counts_filt
  )
})

names(results_filt) <- test_sample

table_recap_filtered_all <- bind_rows(
  lapply(results_filt, function(x) x$table_recap)
)

table_recap_filtered_all

results_filt[[s]]$df_compare %>%
  filter(diff != 0 | is.na(diff)) %>%
  mutate(diff_bin = cut(diff, breaks = 20)) %>%
  mutate(diff_bin = addNA(diff_bin)) %>%
  ggplot(aes(x = diff_bin)) +
  geom_bar(fill = "steelblue", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(results_filt[[s]]$df_scatter,
       aes(x = origine, y = test)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  labs(
    title = paste("Comparaison origine vs test –", s, "only ~15000 genes"),
    x = "Counts origine",
    y = "Counts test"
  )


####################################################################################################################
################### Ploter les delta (différence d'expression entre origine et test) pour les 2 échantillons (voir si l'erreur est cohérente entre les samples)

df_origine <- origine_filt
df_test <- test_counts_filt


#df_test <- df_test %>% select(-rowname) 

df_compare_Ins_1787 <- full_join(
  df_origine %>% select(MGI.symbol, origine = 5),
  df_test %>% select(MGI.symbol, test = 5),
  by = "MGI.symbol"
) %>%
  mutate(diff = test - origine)

df_compare_Ins_1788 <- full_join(
  df_origine %>% select(MGI.symbol, origine = 6),
  df_test %>% select(MGI.symbol, test = 6),
  by = "MGI.symbol"
) %>%
  mutate(diff = test - origine)

df_delta <- inner_join(
  df_compare_Ins_1787 %>% select(MGI.symbol, delta_1787 = diff),
  df_compare_Ins_1788 %>% select(MGI.symbol, delta_1788 = diff),
  by = "MGI.symbol"
)

ggplot(df_delta, aes(x = delta_1787, y = delta_1788)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(
    title = "Comparaison des delta entre origine et test pour Ins1787 et Ins1788",
    x = "delta 1787",
    y = "delta 1788"
  ) +
  theme_minimal()

# code pour tout comparer d'un coup 
compare_two_samples <- function(df_origine, df_test, col1, col2) {
  name1 <- colnames(df_origine)[col1]
  name2 <- colnames(df_origine)[col2]
  
  df_compare_1 <- full_join(
    df_origine %>% select(MGI.symbol, origine = col1),
    df_test    %>% select(MGI.symbol, test    = col1),
    by = "MGI.symbol"
  ) %>%
    mutate(delta_1 = test - origine) 
  
  df_compare_2 <- full_join(
    df_origine %>% select(MGI.symbol, origine = col2),
    df_test    %>% select(MGI.symbol, test    = col2),
    by = "MGI.symbol"
  ) %>%
    mutate(delta_2 = test - origine) 
  
  df_delta <- inner_join(df_compare_1, df_compare_2, by = "MGI.symbol")
  
  ggplot(df_delta, aes(x = delta_1, y = delta_2)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1,
                color = "red", linetype = "dashed") +
    labs(
      title = paste("Comparaison des delta entre origine et test"), #pour", name1, "et", name2
      x = paste("delta", name1),
      y = paste("delta", name2)
    ) +
    theme_minimal()
}

compare_two_samples(df_origine, df_test, 2, 3)

sample_cols <- 2:9
pairs <- combn(sample_cols, 2, simplify = FALSE)

plots <- purrr::map(
  pairs,
  ~ compare_two_samples(df_origine, df_test, .x[1], .x[2])
)

names(plots) <- purrr::map_chr(
  pairs,
  ~ paste0(
    colnames(df_origine)[.x[1]],
    "_vs_",
    colnames(df_origine)[.x[2]]
  )
)

plots[[1]]
show_plots_page <- function(plots, page = 1, n_per_page = 6, ncol = 3) {
  start <- (page - 1) * n_per_page + 1
  end   <- min(page * n_per_page, length(plots))
  
  wrap_plots(plots[start:end], ncol = ncol)
}
show_plots_page(plots, page = 1)
show_plots_page(plots, page = 2)
show_plots_page(plots, page = 3)
show_plots_page(plots, page = 4)
show_plots_page(plots, page = 5)


####################################################################################
### test avec la couverture 
coverage.path <- "/home/marinevernier/Documents/projets/female_cpid_multiregion/with_positive_control/data/m39_M32_final/coverage_per_sample_control.tsv"

coverage <- read_tsv(coverage.path)
coverage$sample <- gsub("-", ".", coverage$sample)
  
compute_slope_pair <- function(df_origine, df_test, col1, col2) {
  name1 <- colnames(df_origine)[col1]
  name2 <- colnames(df_origine)[col2]
  
  df1 <- full_join(
    df_origine %>% select(MGI.symbol, origine = col1),
    df_test %>% select(MGI.symbol, test = col1),
    by = "MGI.symbol"
  ) %>% mutate(delta_1 = test - origine) %>% select(MGI.symbol, delta_1)
  
  df2 <- full_join(
    df_origine %>% select(MGI.symbol, origine = col2),
    df_test %>% select(MGI.symbol, test = col2),
    by = "MGI.symbol"
  ) %>% mutate(delta_2 = test - origine) %>% select(MGI.symbol, delta_2)
  
  df_delta <- inner_join(df1, df2, by = "MGI.symbol") %>% drop_na()
  
  slope <- coef(lm(delta_2 ~ delta_1, data = df_delta))[2]
  
  tibble(
    sample1 = name1,
    sample2 = name2,
    slope = slope
  )
}

sample_cols <- 2:9
pairs <- combn(sample_cols, 2, simplify = FALSE)

slopes <- purrr::map_dfr(
  pairs,
  ~ compute_slope_pair(df_origine, df_test, .x[1], .x[2])
)

slopes_cov <- slopes %>%
  left_join(coverage %>% select(sample, dedup_bam_reads),
            by = c("sample1" = "sample")) %>%
  rename(coverage1 = dedup_bam_reads) %>%
  left_join(coverage %>% select(sample, dedup_bam_reads),
            by = c("sample2" = "sample")) %>%
  rename(coverage2 = dedup_bam_reads) %>%
  mutate(
    coverage_diff = coverage2 - coverage1,
    coverage_ratio = coverage2 / coverage1  
  )


lm_slope <- lm(slope ~ coverage_ratio, data = slopes_cov)
summary(lm_slope)

ggplot(slopes_cov, aes(x = log2(coverage_ratio), y = slope)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Effet de la couverture relative sur la pente des deltas",
    x = "log2(couverture_sample2 / couverture_sample1)",
    y = "pente delta2 ~ delta1"
  ) +
  theme_minimal()

lm_summary <- summary(lm_slope)

coef_beta1 <- round(lm_summary$coefficients["coverage_ratio", "Estimate"], 2)
pval <- signif(lm_summary$coefficients["coverage_ratio", "Pr(>|t|)"], 3)
r2 <- round(lm_summary$r.squared, 2)

# Ajouter sur le plot
ggplot(slopes_cov, aes(x = log2(coverage_ratio), y = slope)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Effet de la couverture relative sur la pente des deltas",
    x = "log2(couverture_sample2 / couverture_sample1)",
    y = "pente delta2 ~ delta1"
  ) +
  theme_minimal() +
  annotate(
    "text", 
    x = min(log2(slopes_cov$coverage_ratio), na.rm = TRUE), 
    y = max(slopes_cov$slope, na.rm = TRUE), 
    label = paste0("Slope = ", coef_beta1,
                   "\nR² = ", r2,
                   "\nP-value = ", pval),
    hjust = 0, vjust = 1,
    size = 4, color = "darkred"
  )

