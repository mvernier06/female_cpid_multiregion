library(tidyverse)
library(FactoMineR)
library(factoextra)
library(M3C)
library(DESeq2)
library(corrplot)
library(readODS)


rm(list=ls())

#### PATHS ####
raw_counts_path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/raw_counts_filtered_allreg_union.csv"
coldata_path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/data/count_data/coldata.ods"


coldata <- read_ods(coldata_path)
#coldata
raw_counts <- read.csv(raw_counts_path) %>%
  column_to_rownames("MGI.symbol")
#
## normalization because PCA shows distances
vst_counts <- varianceStabilizingTransformation(as.matrix(raw_counts))


# raw_counts_path_control <- "~/Documents/cpid_multiregion/female_cpid_multiregion/with_positive_control/data/count_data_positive_control/cpid_multireg_females_counts.txt"
# 
# raw_counts <- read.table(
#   raw_counts_path_control,
#   header = TRUE,
#   sep = "\t",
#   comment.char = "#",
#   stringsAsFactors = FALSE,
#   check.names = FALSE
# )
# 
# counts <- raw_counts[, !(colnames(raw_counts) %in% c(
#   "Chr", "Start", "End", "Strand", "Length"
# ))]
# 
# rownames(counts) <- counts$Geneid
# count_matrix <- counts[, -1]
# 
# colnames(count_matrix) <- sub(
#   "_R2\\.dedup\\.bam$",
#   "",
#   basename(colnames(count_matrix))
# )
# colnames(count_matrix) <- sub(
#   "^([A-Za-z]+)([0-9]+)$",
#   "\\1.\\2",
#   colnames(count_matrix)
# )
# colnames(count_matrix) <- gsub("-", ".", colnames(count_matrix))



plot.path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/graphs_results/1__count_matrix_operation/pca/"
output.path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/data/1__count_matrix_operation/"



## NOT USED
##### PCA ON RAW COUNTS ####
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample


pca(vst_counts, dotsize=3) +
 labs( title="PCA on VST normalized counts")
res <- M3C(vst_counts, removeplots = TRUE, iters=25,
          objective='PAC', fsize=8, lthick=1, dotsize=1.25)
res



setwd(plot.path)
res$plots[[2]]
ggsave(plot=last_plot(), "PAC_raw_counts.png")
res$plots[[4]]
ggsave(plot=last_plot(), "RCSI_raw_counts.png")
res$plots[[3]]
ggsave(plot=last_plot(), "pval_raw_counts.png")
data <- res$realdataresults[[5]]$ordered_data
annon <- res$realdataresults[[5]]$ordered_annotation
ccmatrix <- res$realdataresults[[5]]$consensus_matrix
head(annon)
pca(data,labels=annon$consensuscluster,legendtextsize = 10,axistextsize = 10,dotsize=2) +
 labs( title="PCA on VST normalized counts (Group: 1-Ins, 2-BLA, 3-NAc, 4-Hb, 5-DRN+VTA)")
ggsave(plot=last_plot(), "PCA_raw_counts.png")


#### HEAT MAP (fig1)####
df <- scale(vst_counts)
df

# sort columns by region, tp and sham/cuff
new_order <- coldata %>% arrange(reg, timepoint, desc(group))
df2 <- df[, rownames(new_order)]

# Calculate the correlation matrix
cor_matrix <- cor(df2)
# Create a basic correlation heatmap using corrplot
png(file="heatmap_raw_counts_corrplot.png", width = 4000, height = 4000, units = "px", res = 100)
corrplot(cor_matrix, method = "color")
dev.off()

#### T-SNE ON RAW COUNTS ####
vst_counts
setwd(plot.path)
setwd("../")
set.seed(123)
t <- tsne(vst_counts, labels=coldata$reg, legendtextsize = 15,axistextsize = 15,dotsize=3,
          colvec = "black") +
  labs( title="T-sne on VST normalized counts")
t
ggsave(plot=t, "tsne_raw_counts.png")
?tsne


## test: run 10X to see differences
#for(i in 1:10){
#  ti <- tsne(vst_counts, labels=coldata$reg, legendtextsize = 15,axistextsize = 15,dotsize=3) +
#    labs( title="T-sne on VST normalized counts")
#  ggsave(plot = ti, paste0("tsne_raw_counts", i,".png"))
#}
