library(tidyverse)
library(DESeq2)
library(readODS)


rm(list=ls())

setwd("//wsl.localhost/Ubuntu/home/marinevernier/projets/cpid_multiregion/")
#### PATHS ####
raw_counts.path <- "female_cpid_multiregion/with_positive_control/data/count_data_positive_control/cpid_multireg_females_counts.txt"
coldata.path <- "female_cpid_multiregion/with_positive_control/data/count_data_positive_control/coldata_control.ods"
output.path <- "female_cpid_multiregion/with_positive_control/data/2__differential_expression_analysis/CPID_sham_vs_cuff_betaprior.csv" # file to save results

#### Formatting raw counts ####
raw_counts <- read.table(raw_counts.path, header = TRUE,
                         sep = "\t",
                         comment.char = "#",
                         stringsAsFactors = FALSE,
                         check.names = FALSE)

raw_counts <- raw_counts[, !(colnames(raw_counts) %in% c(
  "Chr", "Start", "End", "Strand", "Length"
))]

rownames(raw_counts) <- raw_counts$Geneid
#raw_counts <- raw_counts[, -1]

colnames(raw_counts) <- sub(
  "_R2\\.dedup\\.bam$",
  "",
  basename(colnames(raw_counts))
)
colnames(raw_counts) <- sub(
  "^([A-Za-z]+)([0-9]+)$",
  "\\1.\\2",
  colnames(raw_counts)
)
colnames(raw_counts) <- gsub("-", ".", colnames(raw_counts))

#### Formatting metadata ####
coldata <- read_ods(coldata.path)
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample
coldata$timepoint <- as.factor(coldata$timepoint) 
coldata$reg <- as.factor(coldata$reg)
coldata$sample <- as.factor(coldata$sample)
coldata <- coldata %>% arrange(by = sample)
coldata %>% filter(timepoint == 1) %>% nrow

#### Differential expression analysis of sham vs cuff ####
Deseq2MultiReg <- function(regionList, tpList){
  
  for (tp in tpList){
    
    # Select rows of the timepoint in coldata
    metadata <- coldata[which(coldata$timepoint == tp),]
    
    # Select columns of the timepoint in raw_counts
    list <- metadata$sample[which(metadata$timepoint == tp)]
    df <- raw_counts[which(colnames(raw_counts) %in% list,)]
    
    for (region in regionList){
      
      print(paste0("Computing stastistics about TP", tp, " of ", region))
      
      # Compute pval, qval and lfc related to the difference of sham vs cuff for each timepoint per region
      dds <- DESeqDataSetFromMatrix(countData=df[, grepl(region, names(df))], 
                                    colData=metadata[str_detect(metadata$reg, region),], 
                                    design=~group)
      
      dds <- DESeq(object = dds, betaPrior = TRUE)
      
      res_tp <- results(dds, contrast=c("group", "cuff", "sham"), 
                        cooksCutoff=TRUE, independentFiltering=TRUE, alpha=0.05)
      
      # Rearrange data and save results at the end of raw_counts
      res_tp <-  as.data.frame(res_tp[, c('log2FoldChange', 'pvalue', 'padj'), drop = FALSE])
      colnames(res_tp) <- c(
        paste0(region, "_log2fc_tp", tp),
        paste0(region, "_pval_tp", tp),
        paste0(region, "_padj_tp", tp)
      )
      
      raw_counts <- cbind(raw_counts, res_tp)
    }
  }
  
  write.csv(raw_counts, output.path, row.names=FALSE)
}

#### Call to the function ####
Deseq2MultiReg(unique(coldata$reg), unique(coldata$timepoint))

#####################################################################################################################################

rm(list=ls())

counts.path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/with_positive_control/data/2__differential_expression_analysis/CPID_sham_vs_cuff_betaprior.csv"
annot_table.path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/data/count_data/annotation_final.csv"
output.path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/with_positive_control/data/2__differential_expression_analysis/annotated_counts.rds" # where to save results

# Annotate counts
counts <- read.csv(counts.path)
counts$Geneid <-  str_replace(counts$Geneid, "\\..*", "") 

ens2symbol <- read.csv(annot_table.path, sep=";")
ens2symbol <- ens2symbol %>% dplyr::select(!X)

annotated_counts <- inner_join(counts, ens2symbol, by=c("Geneid"="Gene.stable.ID"))
annotated_counts <- annotated_counts %>% relocate(MGI.symbol, .before=Geneid) 
# remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)
annotated_counts <- annotated_counts[!duplicated(annotated_counts$MGI.symbol),]

write_rds(annotated_counts, output.path)
output.path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/with_positive_control/data/2__differential_expression_analysis/annotated_counts.csv"
write_csv(annotated_counts, output.path)


#######################################################################################################################################
# attention, j'ai creer à la main deux .ods ne gardant que les samples que je veux comparer, en partant de annotated_counts.csv

origine <- "/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/with_positive_control/data/samples_interet_originel.ods"
test <-"/home/marinevernier/Documents/cpid_multiregion/female_cpid_multiregion/with_positive_control/data/samples_interet_test.ods"

df_origine <- read_ods(origine)
df_test <- read_ods(test)



df_compare_Ins_1787 <- full_join(
  df_origine %>% select(MGI.symbol, origine = 2),
  df_test %>% select(MGI.symbol, test = 2),
  by = "MGI.symbol"
) %>%
  mutate(diff = test - origine)

df_compare_plot <- df_compare_Ins_1787 %>%
  # Garder seulement les gènes avec diff != 0 ou NA
  filter(diff != 0 | is.na(diff)) %>%
  mutate(diff_bin = cut(diff, breaks = 20))  # 30 bins pour les diff numériques

# Ajouter les NA comme catégorie
df_compare_plot$diff_bin <- addNA(df_compare_plot$diff_bin)

ggplot(df_compare_plot, aes(x = diff_bin)) +
  geom_bar(fill = "steelblue", color = "black") +
  labs(
    title = "Histogramme des différences entre les résultats d'origines et le control positif",
    x = "diff (counts test - counts origine)",
    y = "Nombre de gènes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


df_scatter <- df_origine %>%
  select(MGI.symbol, origine = 2) %>%
  full_join(df_test %>% select(MGI.symbol, test = 2),
            by = "MGI.symbol") %>%
  # Supprimer les lignes avec NA
  filter(!is.na(origine) & !is.na(test))

# Scatter plot log-log
ggplot(df_scatter, aes(x = origine, y = test)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Comparaison des counts entre origine et test (contrôle positif)",
    x = "Counts origine",
    y = "Counts test"
  ) +
  theme_minimal()

summary(df_compare_Ins_1787$origine)
summary(df_compare_Ins_1787$test)

cor(df_scatter$origine, df_scatter$test, method = "pearson") 

n_total <- nrow(df_compare_Ins_1787)

table_recap <- tibble(
  categorie = c(
    "diff = 0",
    "diff ≠ 0",
    "diff = NA",
    "test = NA & origine = 0",
    "origine = NA & test = 0"
  ),
  n_genes = c(
    sum(df_compare_Ins_1787$diff == 0, na.rm = TRUE),
    sum(df_compare_Ins_1787$diff != 0, na.rm = TRUE),
    sum(is.na(df_compare_Ins_1787$diff)),
    sum(is.na(df_compare_Ins_1787$test) & df_compare_Ins_1787$origine == 0, na.rm = TRUE),
    sum(is.na(df_compare_Ins_1787$origine) & df_compare_Ins_1787$test == 0, na.rm = TRUE)
  )
) %>%
  mutate(
    percent = 100 * n_genes / n_total
  )

table_recap



df_compare_Ins_1788 <- full_join(
  df_origine %>% select(MGI.symbol, origine = 3),
  df_test %>% select(MGI.symbol, test = 3),
  by = "MGI.symbol"
) %>%
  mutate(diff = test - origine)

df_compare_plot <- df_compare_Ins_1788 %>%
  # Garder seulement les gènes avec diff != 0 ou NA
  filter(diff != 0 | is.na(diff)) %>%
  mutate(diff_bin = cut(diff, breaks = 20))  # 30 bins pour les diff numériques

# Ajouter les NA comme catégorie
df_compare_plot$diff_bin <- addNA(df_compare_plot$diff_bin)

ggplot(df_compare_plot, aes(x = diff_bin)) +
  geom_bar(fill = "steelblue", color = "black") +
  labs(
    title = "Histogramme des différences entre les résultats d'origines et le control positif",
    x = "diff (counts test - counts origine)",
    y = "Nombre de gènes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


df_scatter <- df_origine %>%
  select(MGI.symbol, origine = 3) %>%
  full_join(df_test %>% select(MGI.symbol, test = 3),
            by = "MGI.symbol") %>%
  # Supprimer les lignes avec NA
  filter(!is.na(origine) & !is.na(test))

# Scatter plot log-log
ggplot(df_scatter, aes(x = origine, y = test)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Comparaison des counts entre origine et test (contrôle positif Ins 1788)",
    x = "Counts origine",
    y = "Counts test"
  ) +
  theme_minimal()

summary(df_compare_Ins_1788$origine)
summary(df_compare_Ins_1788$test)

cor(df_scatter$origine, df_scatter$test, method = "pearson") 

n_total <- nrow(df_compare_Ins_1788)

table_recap <- tibble(
  categorie = c(
    "diff = 0",
    "diff ≠ 0",
    "diff = NA",
    "test = NA & origine = 0",
    "origine = NA & test = 0"
  ),
  n_genes = c(
    sum(df_compare_Ins_1788$diff == 0, na.rm = TRUE),
    sum(df_compare_Ins_1788$diff != 0, na.rm = TRUE),
    sum(is.na(df_compare_Ins_1788$diff)),
    sum(is.na(df_compare_Ins_1788$test) & df_compare_Ins_1788$origine == 0, na.rm = TRUE),
    sum(is.na(df_compare_Ins_1788$origine) & df_compare_Ins_1788$test == 0, na.rm = TRUE)
  )
) %>%
  mutate(
    percent = 100 * n_genes / n_total
  )

table_recap


##########################################################################################################################################
#################### En ne gardant que les gènes utilisés dans l'analyse (après filtrage des gènes faiblement exprimés)##################

filtered_genes <- "cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
origine <- "female_cpid_multiregion/with_positive_control/data/samples_interet_originel.ods"
test <-"female_cpid_multiregion/with_positive_control/data/samples_interet_test.ods"

filtered <- read_csv(filtered_genes)
counts_origine <- read_ods(origine)
counts_test <- read_ods(test)

genes_kept <- filtered %>% pull(MGI.symbol)

# filtrer origine et test
counts_origine_filt <- counts_origine %>%
  filter(MGI.symbol %in% genes_kept)

counts_test_filt <- counts_test %>%
  filter(MGI.symbol %in% genes_kept)
all(filtered$MGI.symbol %in% counts_origine_filt$MGI.symbol)
all(filtered$MGI.symbol %in% counts_test_filt$MGI.symbol)

genes_missing_in_test <- setdiff(
  filtered$MGI.symbol,
  counts_test$MGI.symbol
)

length(genes_missing_in_test)
head(genes_missing_in_test)

df_compare_Ins_1787 <- full_join(
  counts_origine_filt %>% select(MGI.symbol, origine = 2),
  counts_test_filt %>% select(MGI.symbol, test = 2),
  by = "MGI.symbol"
) %>%
  mutate(diff = test - origine)

df_compare_plot <- df_compare_Ins_1787 %>%
  # Garder seulement les gènes avec diff != 0 ou NA
  filter(diff != 0 | is.na(diff)) %>%
  mutate(diff_bin = cut(diff, breaks = 20))  # 30 bins pour les diff numériques

# Ajouter les NA comme catégorie
df_compare_plot$diff_bin <- addNA(df_compare_plot$diff_bin)

ggplot(df_compare_plot, aes(x = diff_bin)) +
  geom_bar(fill = "steelblue", color = "black") +
  labs(
    title = "Histogramme des différences entre les résultats d'origines et le control positif",
    x = "diff (counts test - counts origine)",
    y = "Nombre de gènes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


df_scatter <- counts_origine_filt %>%
  select(MGI.symbol, origine = 2) %>%
  full_join(counts_test_filt %>% select(MGI.symbol, test = 2),
            by = "MGI.symbol") %>%
  # Supprimer les lignes avec NA
  filter(!is.na(origine) & !is.na(test))

# Scatter plot log-log
ggplot(df_scatter, aes(x = origine, y = test)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Comparaison des counts entre origine et test (contrôle positif)",
    x = "Counts origine",
    y = "Counts test"
  ) +
  theme_minimal()

summary(df_compare_Ins_1787$origine)
summary(df_compare_Ins_1787$test)

cor(df_scatter$origine, df_scatter$test, method = "pearson") 

n_total <- nrow(df_compare_Ins_1787)

table_recap <- tibble(
  categorie = c(
    "diff = 0",
    "diff ≠ 0",
    "diff = NA",
    "test = NA & origine = 0",
    "origine = NA & test = 0"
  ),
  n_genes = c(
    sum(df_compare_Ins_1787$diff == 0, na.rm = TRUE),
    sum(df_compare_Ins_1787$diff != 0, na.rm = TRUE),
    sum(is.na(df_compare_Ins_1787$diff)),
    sum(is.na(df_compare_Ins_1787$test) & df_compare_Ins_1787$origine == 0, na.rm = TRUE),
    sum(is.na(df_compare_Ins_1787$origine) & df_compare_Ins_1787$test == 0, na.rm = TRUE)
  )
) %>%
  mutate(
    percent = 100 * n_genes / n_total
  )

table_recap



df_compare_Ins_1788 <- full_join(
  counts_origine_filt %>% select(MGI.symbol, origine = 3),
  counts_test_filt %>% select(MGI.symbol, test = 3),
  by = "MGI.symbol"
) %>%
  mutate(diff = test - origine)

df_compare_plot <- df_compare_Ins_1788 %>%
  # Garder seulement les gènes avec diff != 0 ou NA
  filter(diff != 0 | is.na(diff)) %>%
  mutate(diff_bin = cut(diff, breaks = 20))  # 30 bins pour les diff numériques

# Ajouter les NA comme catégorie
df_compare_plot$diff_bin <- addNA(df_compare_plot$diff_bin)

ggplot(df_compare_plot, aes(x = diff_bin)) +
  geom_bar(fill = "steelblue", color = "black") +
  labs(
    title = "Histogramme des différences entre les résultats d'origines et le control positif",
    x = "diff (counts test - counts origine)",
    y = "Nombre de gènes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


df_scatter <- counts_origine_filt %>%
  select(MGI.symbol, origine = 3) %>%
  full_join(counts_test_filt %>% select(MGI.symbol, test = 3),
            by = "MGI.symbol") %>%
  # Supprimer les lignes avec NA
  filter(!is.na(origine) & !is.na(test))

# Scatter plot log-log
ggplot(df_scatter, aes(x = origine, y = test)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Comparaison des counts entre origine et test (contrôle positif Ins 1788)",
    x = "Counts origine",
    y = "Counts test"
  ) +
  theme_minimal()

summary(df_compare_Ins_1788$origine)
summary(df_compare_Ins_1788$test)

cor(df_scatter$origine, df_scatter$test, method = "pearson") 

n_total <- nrow(df_compare_Ins_1788)

table_recap <- tibble(
  categorie = c(
    "diff = 0",
    "diff ≠ 0",
    "diff = NA",
    "test = NA & origine = 0",
    "origine = NA & test = 0"
  ),
  n_genes = c(
    sum(df_compare_Ins_1788$diff == 0, na.rm = TRUE),
    sum(df_compare_Ins_1788$diff != 0, na.rm = TRUE),
    sum(is.na(df_compare_Ins_1788$diff)),
    sum(is.na(df_compare_Ins_1788$test) & df_compare_Ins_1788$origine == 0, na.rm = TRUE),
    sum(is.na(df_compare_Ins_1788$origine) & df_compare_Ins_1788$test == 0, na.rm = TRUE)
  )
) %>%
  mutate(
    percent = 100 * n_genes / n_total
  )

table_recap

df_scatter %>%
  mutate(
    log2FC = log2((test + 1) / (origine + 1))
  ) %>%
  summarise(
    pct_big_change = mean(abs(log2FC) > 0.1) * 100,
    pct_small_change = mean(abs(log2FC) < 0.1) * 100,
    pct_very_small = mean(abs(log2FC) < 0.01) * 100
  )
