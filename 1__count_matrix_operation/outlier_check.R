## Dans cette partie je vérifie si les "outliers" identifiés par PCA et T-sne sont vraiment bizarres 

library(tidyverse)
library(readODS)
library(ggplot2)
library(DESeq2)
library(ggrepel)

rm(list=ls())

setwd("/home/marinevernier/Documents/projets/")

#### PATHS ####
raw_counts_path <- "female_cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
coldata_path <- "female_cpid_multiregion/data/counts_m39_M32/coldata.ods"
plot.path <- "female_cpid_multiregion/graphs_results/1__count_matrix_operation/outlier_check/"
coverage.path <- "female_cpid_multiregion/data/counts_m39_M32/coverage_per_sample.tsv"


counts <- read_csv(raw_counts_path)
coldata <- read_ods(coldata_path)
coverage <- read_tsv(coverage.path)
coverage$sample <- gsub("-", ".", coverage$sample)
coverage$sample <- sub(
  "^([A-Za-z]+)([0-9]+)$",
  "\\1.\\2",
  coverage$sample
)

setwd(plot.path)

counts_long <- counts %>%
  pivot_longer(
    cols = -MGI.symbol,
    names_to = "sample",
    values_to = "expression"
  )

genes_present_testes <- c("Grp", "Drd2", "Drd1", "Camk2a", "Camk2b", "Mog", "Th", "Thy1")

counts <- counts %>% column_to_rownames("MGI.symbol")
"vmat" %in% rownames(counts)


Mog <- counts_long %>%
  filter(MGI.symbol == "Mog")
Mog_annot <- Mog %>%
  left_join(coldata, by = "sample")

ggplot(Mog_annot, aes(x = reg, y = expression, label = sample)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(vjust = -0.5, size = 3) +
  theme_classic() +
  labs(
    title = "Mog – détection d'échantillons aberrants",
    x = "Région",
    y = "Expression"
  )

plots <- list()

for (gene in genes_present_testes) {
  
  counts_gene <- counts_long %>%
    filter(MGI.symbol == gene)
  
  count_cov <- counts_gene %>%
    left_join(coldata, by = "sample") %>%
    left_join(coverage, by = "sample") %>%
    mutate(
      expression_covnorm = (expression / dedup_bam_reads) * 1e6
    )
  
  
    
    p<- ggplot(
      count_cov,
      aes(x = reg, y = expression_covnorm, label = sample)
      ) +
        geom_point(size = 3, alpha = 0.8) +
        geom_text(vjust = -0.5, size = 3) +
        theme_classic() +
        labs(
          title = paste0(gene, " – expression normalisée par la couverture"),
          x = "Région",
          y = "CPM (dedup reads)"
        )
    plots[[gene]] <- p
    filename <- paste0(gene, "_expression.png")
    ggsave(plot=last_plot(), filename)
}
plots[["Thy1"]]



############################################################################################################

all(colnames(counts) == coldata$sample)

mean_intra_cor <- function(counts, coldata, region,
                            exclude_sample = NULL,
                            method = "pearson") {
  
  meta_sub <- coldata %>% filter(reg == region)
  
  if (!is.null(exclude_sample)) {
    meta_sub <- meta_sub %>%
      filter(!sample %in% exclude_sample)
  }
  
  if (nrow(meta_sub) < 3) return(NA_real_)
  
  counts_sub <- counts[, meta_sub$sample, drop = FALSE]
  
  cor_mat <- cor(counts_sub, method = method)
  
  mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
}



outlier_sample <- c("Ins.1837", "Hb.1839")

regions <- unique(coldata$reg)

cor_summary <- tibble(
  reg = regions,
  cor_with_outlier = purrr::map_dbl(
    regions,
    ~ mean_intra_cor(counts, coldata, .x, method = "spearman")
  ),
  cor_without_outlier = purrr::map_dbl(
    regions,
    ~ mean_intra_cor(counts, coldata, .x,
                     exclude_sample = outlier_sample,
                     method = "spearman")
  )
) %>%
  mutate(
    delta = cor_without_outlier - cor_with_outlier
  )

cor_summary


#################################################################################################################

outliers <- c("Ins.1837", "Hb.1839", "Hb.2049")

lib_size <- colSums(counts)

lib_df <- tibble(
  sample = names(lib_size),
  total_counts = lib_size
) %>%
  left_join(coldata, by = "sample") %>%
  mutate(
    is_outlier = sample %in% outliers
  )

ggplot(lib_df,
       aes(x = reorder(sample, total_counts),
           y = total_counts,
           fill = is_outlier)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = c("FALSE" = "grey70",
               "TRUE"  = "red3"),
    guide = "none"
  ) +
  labs(
    title = "Total library size per sample",
    x = NULL,
    y = "Total counts"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 4),
    axis.text.x = element_text(size = 8),
    plot.title  = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank()
  )
ggsave(plot=last_plot(), "library_size.png")



lib_df_corrected <- tibble(
  sample = names(lib_size),
  total_counts = lib_size
) %>%
  left_join(coldata, by = "sample") %>%
  left_join(coverage, by = "sample") %>%
  mutate(
    is_outlier = sample %in% outliers,
    counts_corrected = (total_counts/dedup_bam_reads)*1e6
  )

ggplot(lib_df_corrected,
       aes(x = reorder(sample, counts_corrected),
           y = counts_corrected,
           fill = is_outlier)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = c("FALSE" = "grey70",
               "TRUE"  = "red3"),
    guide = "none"
  ) +
  labs(
    title = "Total library size per sample, correted by coverage",
    x = NULL,
    y = "Total counts corrected by coverage"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 8),
    plot.title  = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank()
  )
ggsave(plot=last_plot(), "corrected_library_size.png")  

df <- tibble(
  sample = names(lib_size),
  total_counts = lib_size
) %>%
  left_join(coverage, by = "sample") 

ggplot(df, 
       aes(x = fastq_reads ,
       y = total_counts))+
  geom_point() +
  geom_text_repel(
    data = subset(df, sample %in% c("Ins.1837", "Hb.1839", "Hb.2049", "Ins.2073", "Ins.2088", "Ins.1850", "Ins2079", "Nac.1837", "Ins.2072", "Ins.2079", "Ins.1834", "Hb.2074", "Ins.1840")),
    aes(label = sample),
    size = 4,
    fontface = "bold"
  ) +
  labs(
       title = "Total library size and total fatsq reads",
       x = "fastq reads",
       y = "total counts"
     ) 
ggsave(plot=last_plot(), "total_counts_vs_number_reads.png")

#########################################################################################################
###################  test avec le RIN ##################







