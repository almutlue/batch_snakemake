##Characterize Simulation - Reduced Dim representation :
# Aim: Add reduced dim representations to simulations

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

set.seed(1000)

## Show arguments
print(sim_sce)
print(outputfile)
print(params)
print(out_path)


#libraries
suppressPackageStartupMessages({
    library(scater)
    library(CellMixS)
    library(scran)
})

# Read in data
sim_sce <- readRDS(file = sim_sce)
params <- readRDS(file = params)
k <- params[["k"]]

## Create directory for outputfiles
dir.create(out_path, showWarnings = FALSE)


#Remove group label
sim_sce$sample_id <- gsub('\\..', '', sim_sce$sample_id)

### ------------ Standardize normalized counts-----------###
clusters <- quickCluster(sim_sce, use.ranks=FALSE)
table(clusters)
sim_sce <- computeSumFactors(sim_sce, min.mean=0.1, cluster=clusters)
sim_sce <-  logNormCounts(sim_sce)

### --------------Reduced dim representation------------###
sim_sce <- runUMAP(sim_sce, ntop = 1000)
sim_sce <- runPCA(sim_sce, ntop = 1000, ncomponents = 10)
sim_sce <- runTSNE(sim_sce, ntop = 1000)


### --------------Run cms-----------------------------###
sim_sce <- cms(sim_sce, group = "sample_id", k = k, k_min = 10, n_dim = 10)

#Save outputfile
saveRDS(sim_sce, outputfile)