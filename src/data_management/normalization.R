##Normalization and preprocessing:
#Aim: Make sure all dataset have counts based on the same way normalization,
#all defined variables and all dimension reductions.

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
  library(scran)
  library(scater)
})


# Read in data
sce <- readRDS(file = data)
params <- readRDS(file = params)

# vars
batch <- params[["batch"]]
celltype <- params[["celltype"]]
sample <- params[["Sample"]]

### ------------ Check parameter-------------------------###
if( !all(c(batch, celltype, sample) %in% names(colData(sce))) ){
  stop("Unspecified parameters - Please check your meta-file!")
}

### ------------ Standardize normalized counts-----------###
clusters <- quickCluster(sce, use.ranks=FALSE)
table(clusters)
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
sce <-  logNormCounts(sce)


### ------------ Check reduced dimensions----------------###
if( !"PCA" %in% reducedDimNames(sce) ){
  sce <- runPCA(sce, ncomponents = 10, ntop = 1000)
}

if( !"UMAP" %in% reducedDimNames(sce) ){
  sce <- runUMAP(sce)
}

if( !"TSNE" %in% reducedDimNames(sce) ){
  sce <- runTSNE(sce)
}


### -------------- save sce object ----------------------###
saveRDS(sce, file = outputfile)


