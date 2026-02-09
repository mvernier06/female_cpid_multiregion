library(tidyverse)
library(readODS)
library(tidyr)

rm(list = ls())

setwd("//wsl.localhost/Ubuntu/home/marinevernier/projets/cpid_multiregion/")


origine.path <- "/home/marinevernier/projets/cpid_multiregion/cpid_multiregion/data/2__differential_expression_analysis/annotated_counts.csv"
test.path <- "/home/marinevernier/projets/cpid_multiregion/female_cpid_multiregion/with_positive_control/m39_M32_magis/data_2/cpid_multireg_control_counts.txt"

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
annot_table.path <- "female_cpid_multiregion/data/count_data/annotation_final.csv"

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
  "cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv",
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

