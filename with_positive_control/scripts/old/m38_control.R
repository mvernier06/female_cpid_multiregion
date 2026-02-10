library(tidyverse)
library(readODS)

rm(list = ls())

setwd("//wsl.localhost/Ubuntu/home/marinevernier/projets/cpid_multiregion/")


v1.path <- "/home/marinevernier/projets/cpid_multiregion/female_cpid_multiregion/with_positive_control/m38/counts_control/control_v1/cpid_multireg_control_counts.txt"
v2.path <- "/home/marinevernier/projets/cpid_multiregion/female_cpid_multiregion/with_positive_control/m38/counts_control/control_v2/cpid_multireg_control_counts.txt"


v1 <- read.table(v1.path, header = TRUE,
                 sep = "\t",
                 comment.char = "#",
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

v1 <- v1[, !(colnames(v1) %in% c(
  "Chr", "Start", "End", "Strand", "Length"
))]

rownames(v1) <- v1$Geneid
v1$Geneid <- NULL

colnames(v1) <- sub(
  "_R2\\.dedup\\.bam$",
  "",
  basename(colnames(v1))
)
colnames(v1) <- sub(
  "^([A-Za-z]+)([0-9]+)$",
  "\\1.\\2",
  colnames(v1)
)
colnames(v1) <- gsub("-", ".", colnames(v1))

annot_table.path <- "female_cpid_multiregion/data/count_data/annotation_final.csv"

v1 <-  rownames_to_column(v1)
v1$rowname <-  str_replace(v1$rowname, "\\..*", "") 

ens2symbol <- read.csv(annot_table.path, sep=";")
ens2symbol <- ens2symbol %>% dplyr::select(!X)

counts_v1 <- inner_join(v1, ens2symbol, by=c("rowname"="Gene.stable.ID"))
counts_v1 <- counts_v1 %>% relocate(MGI.symbol, .before=rowname) 
# remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)
counts_v1 <- counts_v1[!duplicated(counts_v1$MGI.symbol),]


v2 <- read.table(v2.path, header = TRUE,
                 sep = "\t",
                 comment.char = "#",
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

v2 <- v2[, !(colnames(v2) %in% c(
  "Chr", "Start", "End", "Strand", "Length"
))]

rownames(v2) <- v2$Geneid
v2$Geneid <- NULL

colnames(v2) <- sub(
  "_R2\\.dedup\\.bam$",
  "",
  basename(colnames(v2))
)
colnames(v2) <- sub(
  "^([A-Za-z]+)([0-9]+)$",
  "\\1.\\2",
  colnames(v2)
)
colnames(v2) <- gsub("-", ".", colnames(v2))


v2 <-  rownames_to_column(v2)
v2$rowname <-  str_replace(v2$rowname, "\\..*", "") 



counts_v2 <- inner_join(v2, ens2symbol, by=c("rowname"="Gene.stable.ID"))
counts_v2 <- counts_v2 %>% relocate(MGI.symbol, .before=rowname) 
# remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)
counts_v2 <- counts_v2[!duplicated(counts_v2$MGI.symbol),]


counts_v1$rowname <- NULL
rownames(counts_v1) <- counts_v1$MGI.symbol
counts_v1$MGI.symbol <- NULL
counts_v2$rowname <- NULL
rownames(counts_v2) <- counts_v2$MGI.symbol
counts_v2$MGI.symbol <- NULL

diff_counts <- counts_v2 - counts_v1
summary(as.vector(diff_counts))

sample <- colnames(counts_v1)[1]

df <- data.frame(
  v1 = counts_v1[, sample],
  v2 = counts_v2[, sample]
)

ggplot(df, aes(v1, v2)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle(sample)

correlations <- sapply(colnames(counts_v1), function(s) {
  cor(counts_v1[, s], counts_v2[, s], method = "spearman")
})

correlations

# Conclusion : l'alignement et comptage sont bien déterministes 




