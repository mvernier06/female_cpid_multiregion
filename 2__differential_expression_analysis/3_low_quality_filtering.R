library(tidyverse)

regionlist <- c("ACC", "Hb", "Ins", "Nac")

#### PATHS ####
annotated_counts.path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/annotated_counts.rds"
output.path <- "~/Documents/cpid_multiregion/female_cpid_multiregion/data/2__differential_expression_analysis/"

annotated_counts <- read_rds(annotated_counts.path)



#### low quality filtering ####
for(reg in regionlist){
  # load unnormalised counts as dataframe 
  counts <- annotated_counts %>%
    distinct(MGI.symbol, .keep_all = TRUE) %>%
    column_to_rownames("MGI.symbol") %>%
    select(contains(paste0(reg, ".")))
  cpm <- apply(counts, 2, function(x) (x/sum(x))*1000000)
  raw_counts_filtered <- counts[rowSums(cpm < 1) < (dim(cpm)[2]/5),]
  counts_filtered <- annotated_counts[which(annotated_counts$MGI.symbol %in% rownames(raw_counts_filtered)),]
  saveRDS(counts_filtered, file=paste0(output.path, "counts_filtered_", reg, ".rds")) # save a count file for each reg
  assign(paste0("counts_filtered_", reg), counts_filtered, envir = .GlobalEnv)
}

df <- list(counts_filtered_ACC %>% select(MGI.symbol,contains("ACC")), 
           counts_filtered_Hb %>% select(MGI.symbol,contains("Hb")),
           counts_filtered_Ins %>% select(MGI.symbol,contains("Ins")),
           counts_filtered_Nac %>% select(MGI.symbol,contains("Nac"))
) %>% 
  purrr::reduce(full_join, by = "MGI.symbol")
saveRDS(df, file=paste0(output.path, "annotated_counts_filtered.rds")) # save a count file for all reg

# save filtered raw counts as csv for Magis dimension reduction #
# union des gènes de chaque région
filtered_genes <- unique(c(counts_filtered_ACC$MGI.symbol,
                           counts_filtered_Hb$MGI.symbol,
                           counts_filtered_Ins$MGI.symbol,
                           counts_filtered_Nac$MGI.symbol))

raw_counts_filtered <- annotated_counts %>%
  distinct(MGI.symbol, .keep_all = TRUE) %>%
  select(MGI.symbol, contains(".")) %>%
  filter(MGI.symbol %in% filtered_genes)

write.csv(raw_counts_filtered, file=paste0(output.path, "raw_counts_filtered_allreg_union.csv"), row.names = FALSE)


#####################################################################################################################


# print n of each file:
paste('Acc:', nrow(counts_filtered_ACC))
paste('Hb:', nrow(counts_filtered_Hb))
paste('Ins:', nrow(counts_filtered_Ins))
paste('Nac:', nrow(counts_filtered_Nac))
paste("intersect_Acc_Hb:",
      length(intersect(counts_filtered_Hb$MGI.symbol, counts_filtered_ACC$MGI.symbol)))
paste("intersect_Acc_Ins:",
      length(intersect(rownames(counts_filtered_Ins), rownames(counts_filtered_ACC))))
paste("intersect_Acc_Nac:",
      length(intersect(rownames(counts_filtered_Nac), rownames(counts_filtered_ACC))))
paste("intersect_Hb_Ins:",
      length(intersect(rownames(counts_filtered_Hb), rownames(counts_filtered_Ins))))
paste("intersect_Hb_Nac:",
      length(intersect(rownames(counts_filtered_Hb), rownames(counts_filtered_Nac))))
paste("intersect_Ins_Nac:",
      length(intersect(rownames(counts_filtered_Ins), rownames(counts_filtered_Nac))))
paste("intersect_all:",
      length(Reduce(intersect, list(rownames(counts_filtered_ACC),
                                    rownames(counts_filtered_Hb),
                                    rownames(counts_filtered_Ins),
                                    rownames(counts_filtered_Nac)))))


