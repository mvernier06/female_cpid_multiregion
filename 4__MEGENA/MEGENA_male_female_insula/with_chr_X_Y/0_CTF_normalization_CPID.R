# install.packages("readODS", repos="https://cloud.r-project.org")
library(tidyverse)
library(edgeR)
library(limma)
library(readODS)

project_path <- "/home2020/home/inci/mvernier/cpid_multireg_female/"
setwd(project_path)

print("CTF normalization")
reg <- "Ins"

# PATHS ####
raw_counts_female.path <- "female_cpid_multiregion/data/counts_m39_M32/cpid_multireg_counts.txt"
raw_counts_male.path <- "male_cpid_multiregion/data/count_data/outliers_removed/un_normalised_counts.csv"
annot_table.path <- "female_cpid_multiregion/data/counts_m39_M32/annotation_final.csv"
output.path <- "female_cpid_multiregion/data/4__MEGENA/MEGENA_male_female_insula/with_chr_X_Y/" # folder to save results
dir.create(output.path)
print(paste0("output path: ", output.path))

### Formating female data ###
raw_counts_female <- read.table(raw_counts_female.path, header = TRUE,
                         sep = "\t",
                         comment.char = "#",
                         stringsAsFactors = FALSE,
                         check.names = FALSE)

#recuperer les gènes des chromosomes X et Y pour plus tard 
genes_chrXY_raw <- raw_counts_female %>%
  dplyr::filter(Chr %in% c("chrX", "chrY")) %>%
  dplyr::select(Geneid) %>%
  distinct()
print(paste0("Number of genes on chrX and chrY: ", nrow(genes_chrXY_raw)))
# enlever les versions
genes_chrXY_raw$Geneid <- str_replace(genes_chrXY_raw$Geneid, "\\..*", "")


counts_female <- raw_counts_female[, !(colnames(raw_counts_female) %in% c(
  "Chr", "Start", "End", "Strand", "Length"
))]

# Formating colnames
colnames(counts_female) <- sub(
  "_R2\\.dedup\\.bam$",
  "",
  basename(colnames(counts_female))
)
colnames(counts_female) <- sub(
  "^([A-Za-z]+)([0-9]+)$",
  "\\1.\\2",
  colnames(counts_female)
)
counts_female$Geneid <-  str_replace(counts_female$Geneid, "\\..*", "") 

# annotate genes
ens2symbol <- read.csv(annot_table.path, sep=";")
ens2symbol <- ens2symbol %>% dplyr::select(!X)
counts_female <- inner_join(counts_female, ens2symbol, by=c("Geneid"="Gene.stable.ID"))
counts_female <- counts_female %>% relocate(MGI.symbol, .before=Geneid) 
counts_female <- counts_female[!duplicated(counts_female$MGI.symbol),] # remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)

genes_chrXY <- genes_chrXY_raw %>%
  inner_join(ens2symbol, by = c("Geneid" = "Gene.stable.ID")) %>%
  dplyr::select(MGI.symbol) %>%
  distinct()
print(paste0("Number of genes on chrX and chrY after annotation: ", nrow(genes_chrXY)))

# Keeping only Insula samples and excluding outlier
outlier <- "Ins.1837"
counts_female <- counts_female %>%
  distinct(MGI.symbol, .keep_all = TRUE) %>%
  column_to_rownames("MGI.symbol") %>%
  select(contains(paste0(reg, "."))) %>%
  select(!contains(outlier))
print(dim(counts_female))


### Formating male data ###
raw_counts_male <- read.csv(raw_counts_male.path, sep = ";")
raw_counts_male$X <- str_replace(raw_counts_male$X, "\\..*", "")
raw_counts_male <- inner_join(raw_counts_male, ens2symbol, by = c("X" = "Gene.stable.ID"))  
raw_counts_male <- raw_counts_male[!duplicated(raw_counts_male$MGI.symbol),] # remove duplicates to avoid alluvial errors (ex: pattern ns_ns_ns)

print(dim(raw_counts_male))
filtered_counts_male <- raw_counts_male %>% filter(MGI.symbol %in% rownames(counts_female))
print(dim(filtered_counts_male))

counts_male <- filtered_counts_male %>%
  distinct(MGI.symbol, .keep_all = TRUE) %>%
  column_to_rownames("MGI.symbol") %>%
  select(contains(paste0(reg, "."))) 

print(dim(counts_male))

#keep genes with at least 1 cpm in >= 20% of samples in a single project ###
cpm_female <- apply(counts_female, 2, function(x) (x/sum(x))*1000000)
counts_female_filtered <- counts_female[rowSums(cpm_female < 1) < (dim(cpm_female)[2]/5),]

cpm_male <- apply(counts_male, 2, function(x)(x/sum(x))*1000000)
counts_male_filtered <- counts_male[rowSums(cpm_male < 1) < (dim(cpm_male)[2]/5),]

counts_female_filtered <- rownames_to_column(counts_female_filtered, "MGI.symbol")
counts_male_filtered <- rownames_to_column(counts_male_filtered, "MGI.symbol")

# joining both datasets 
counts_both <- inner_join(counts_female_filtered, counts_male_filtered, by = "MGI.symbol")
counts_both <- column_to_rownames(counts_both, "MGI.symbol")
print(dim(counts_both))

# CTF normalization
lib_size <- base::colSums(counts_both)
norm_factors <- calcNormFactors(object = counts_both, lib.size = lib_size, method = "TMM")
CTF_normalized <- sweep(counts_both, 2, norm_factors, "/")
print(dim(CTF_normalized))

genes_chrXY_in_data <- genes_chrXY %>%
  dplyr::filter(MGI.symbol %in% rownames(CTF_normalized))
print(paste0("Number of genes on chrX and chrY in the normalized data: ", nrow(genes_chrXY_in_data)))
write.table(
  genes_chrXY_in_data,
  file = paste0(output.path, "genes_chrXY_MGI_symbols.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

print(head(CTF_normalized))
# save R object of normalized counts
filename <- paste0("CTF_normalized_counts_", reg, ".Rdata")
save(CTF_normalized, file=paste0(output.path, filename))

print(paste0("CTF normalization done for ", reg, " region. Normalized counts saved in: ", output.path, filename))
