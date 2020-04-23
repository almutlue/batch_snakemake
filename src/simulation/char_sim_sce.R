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
params <- readRDS(file = params)

name <- params[["dataset_name"]]
flag <- ifelse(sim_sce %in% c(paste0('out/sim/', name, '/sim_', name, '_4__2_1_sce.rds'),
                              paste0('out/sim/', name, '/sim_', name, '_4__1_0_sce.rds'),
                              paste0('out/sim/', name, '/sim_', name, '_4__1_sce.rds'),
                              "out/sim/hca/sim_hca_2__2_1_sce.rds"),
               TRUE, FALSE)

if( name %in% "pancreas"){
  flag <- TRUE
}

sim_sce <- readRDS(file = sim_sce)
k <- params[["k"]]
k

## Create directory for outputfiles
dir.create(out_path, showWarnings = FALSE)


#Remove group label
sim_sce$sample_id <- gsub('\\..', '', sim_sce$sample_id)

### ------------ Standardize normalized counts-----------###
#remove empty genes and cells
sim_sce <- sim_sce[,!colSums(assays(sim_sce)[["counts"]]) <= 100]
sim_sce <- sim_sce[!rowSums(assays(sim_sce)[["counts"]]) <= 50,]
clusters <- quickCluster(sim_sce, use.ranks=FALSE)
table(clusters)

if( !flag ){
   sim_sce <- computeSumFactors(sim_sce, min.mean=0.1, cluster=clusters)
}

sim_sce <-  logNormCounts(sim_sce)

### --------------Reduced dim representation------------###
sim_sce <- runPCA(sim_sce, ntop = 1000, ncomponents = 10, exprs_values = "counts")
sim_sce <- runTSNE(sim_sce, ntop = 1000, exprs_values = "counts")
sim_sce <- runUMAP(sim_sce, ntop = 1000, exprs_values = "logcounts")



### --------------Run cms-----------------------------###
if( !flag ){
    sim_sce <- cms(sim_sce, group = "batch_id", k = k, k_min = 50, n_dim = 10)
}else{
    sim_sce$cms <- rep(NA, ncol(sim_sce))
    sim_sce$cms_smooth <- rep(NA, ncol(sim_sce))
}


#Save outputfile
saveRDS(sim_sce, outputfile)
