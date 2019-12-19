#Cellspecific mixing score:
#Aim: Batch characterization

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
  library(CellMixS)
})


# Read in data
sce <- readRDS(file = data)
params <- readRDS(file = params)

# vars
batch <- params[["batch"]]
celltype <- params[["celltype"]]
sample <- params[["Sample"]]
k <- params[["k"]]
k_min <- params[["k_min"]]

### ------------ CellMixS-------------------------###
if( ncol(reducedDims(sce)[["PCA"]]) < 5 ){
  n_dim = ncol(reducedDims(sce)[["PCA"]])
}else{
  n_dim = 5
}

sce <- cms(sce, group = batch, k = k, k_min = 50, n_dim = n_dim)

### -------------- save sce object ----------------------###
saveRDS(sce, file = outputfile)


