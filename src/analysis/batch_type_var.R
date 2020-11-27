##Batch type variance partition:
#Aim: The idea here is to investigate how batch effects can manifest in single-cell
# data by testing what method is necessary to remove them.
# Basically, if batch effects are “simply” mean shifts of expression levels for some
# genes for all the cells in a given celltype and batch, then we should be able to remove
# them by looking at residuals from some linear model involving batch and/or cell type.
# Of course, in a real data set we don’t know the cell type in advance, so this is not
# really intended to be a way of removing batch effects in practice,
# but rather the aim is to figure out how they manifest
# (and thus how we should include them in a simulation),
# based on some data sets where we “know” the batch as well as the cell type.

# We will model the variance explained in different ways and compare the fraction 
# of unexplained variance across samples

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
    library(limma)
    library(dplyr)
    library(scater)
    library(CellMixS)
    library(BiocParallel)
    library(variancePartition)
})


# Read in data
sce <- readRDS(file = data)
params <- readRDS(file = params)

# vars
batch <- params[["batch"]]
celltype <- params[["celltype"]]
sample <- params[["Sample"]]



### ------------ Variance partitioning-------------------------###
run_var_par <- function(sce, x_nam, form){
    #remove zero counts
    if( length(which(rowSums(assays(sce)$logcounts) == 0)) > 0 ){
        sce <- sce[-which(rowSums(assays(sce)$logcounts) == 0),]
    }
    expr <- as.matrix(assays(sce)$logcounts)
    meta_sub <- as.data.frame(colData(sce)[, c(celltype, batch)])
    varPart <- fitExtractVarPartModel(expr, form, meta_sub, BPPARAM=MulticoreParam(4))
    # Add to sce
    rowData(sce)[,paste0("vp_celltype_", x_nam)] <- varPart[[celltype]]
    rowData(sce)[,paste0("vp_batch_", x_nam)] <- varPart[[batch]]
    rowData(sce)[,paste0("vp_residuals_", x_nam)] <- varPart[["Residuals"]]
    sce
}

form <- as.formula(paste0("~ (1|", celltype, ") + (1|", batch, ")"))

no_removal <- as.formula(paste0("~ (1|", celltype, ")"))
adj1 <- as.formula(paste0("~ (1|", celltype, ") + (1|", batch, ")"))
adj2 <- as.formula(paste0("~ (1|(", celltype, ":", batch, "))"))
adj3 <- as.formula(paste0("~ (1|", celltype, ") + (1|", batch, ")+",  "(1|(", celltype, ":", batch, "))"))

sce <- run_var_par(sce, "no_adjust", form = no_removal)
sce <- run_var_par(sce, "Xadj1", form = adj1)


if(length(grep("pancreas", data)) == 0){
    sce <- run_var_par(sce, "Xadj3", form = adj3)
    sce <- run_var_par(sce, "Xadj2", form = adj2)
}




### -------------- save sce object ----------------------###
saveRDS(sce, file = outputfile)


