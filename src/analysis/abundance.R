##Differential abundance:
#Aim: Are celltypes differential abundant between batches.
# Does the batch effect has an effect on the celltype abundance?

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

set.seed(1000)

## Show arguments
print(data)
print(outputfile)
print(params)


#libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(magrittr)
  library(tibble)
  library(dplyr)
  library(tidyr)
})


# Read in data
sce <- readRDS(file = data)
params <- readRDS(file = params)

# vars
batch <- params[["batch"]]
celltype <- params[["celltype"]]
sample <- params[["Sample"]]

### ------------ Differential abundance-------------------------###
meta_tib <- as_tibble(colData(sce)) %>% group_by_at(c(batch, celltype)) %>%
  summarize(n = n()) %>% spread(eval(batch),"n")
meta_df <- as.data.frame(eval(meta_tib))[,-1]
meta_comb <- combn(meta_df, 2, simplify=FALSE)
res <- lapply(meta_comb, function(x){
  cond1 <- names(x)[1]
  cond2 <- names(x)[2]
  rel_abund_diff <- mapply(function(cond1, cond2){
      abs(cond1 - cond2)/(cond1 + cond2)},
      x[,cond1], x[,cond2])
  rel_abund_diff
})

### -------------- save sce object ----------------------###
saveRDS(res, file = outputfile)


