##Batch type:
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

# We first define three functions that we will use to remove the batch effect -
#   the difference between them is which linear model we fit, and thus which residuals we use.
# Throughout, we will assume that X is a matrix with log-expression values,
# ct is a cell type vector, and bt is a batch vector.

# We will see that if the batch effect is a simple (possibly batch- and cell type-specific)
# shift in mean, one of these functions will remove it
# (adj3 should remove it regardless of the complexity/cell type specificity,
#   whereas the other functions will only remove batch effects with a simpler structure).
# If none of these functions remove the batch effect, it is more complex than a mean shift?

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

### ------------ Batch types -------------------------###
## This function will remove the batch effect if it is constant for all cell
## types, and the cell type composition is the same across batches. It will not
## work if the cell type compositions differ, since it is basically just
## equalizing the overall means across the batches
adj1 <- function(X, bt, ct) {
  mm <- model.matrix(~ bt)
  fit <- lmFit(X, mm)
  beta <- fit$coefficients
  beta[, 1] <- 0
  X - beta %*% t(mm)
}

## This function will remove the batch effect if it is constant for all cell
## types, even if the cell type composition is not the same across batches. It
## will not work for cell type-specific batch effects, since we don't include an
## interaction effect between batch and cell type
adj2 <- function(X, bt, ct) {
  mm <- model.matrix(~ ct + bt)
  fit <- lmFit(X, mm)
  beta <- fit$coefficients
  beta[, grep("^bt", colnames(beta), invert = TRUE)] <- 0
  X - beta %*% t(mm)
}

## This function will remove the batch effect even if it is cell type specific
adj3 <- function(X, bt, ct) {
  mm <- model.matrix(~ ct / bt)
  fit <- lmFit(X, mm)
  beta <- fit$coefficients
  beta[, union(grep("^ct", colnames(beta), invert = TRUE),
               grep(":bt", colnames(beta), invert = TRUE))] <- 0
  X - beta %*% t(mm)
}

# Run models
Xadj1 <- adj1(logcounts(sce), colData(sce)[, batch], colData(sce)[, celltype])
Xadj2 <- adj2(logcounts(sce), colData(sce)[, batch], colData(sce)[, celltype])
Xadj3 <- adj3(logcounts(sce), colData(sce)[, batch], colData(sce)[, celltype])

#calculate red dimensions
add_dim_adj <- function(X, name_X){
  pca <- prcomp(t(X), rank. = 10)
  reducedDims(sce)[[paste0("PCA_", name_X)]] <- pca$x
  assays(sce)$removed_batch <- X
  sce <- runUMAP(sce, ntop = 2000, exprs_values = "removed_batch",
                 name = paste0("UMAP_", name_X))
  sce
}

sce <- add_dim_adj(Xadj1, "Xadj1")
sce <- add_dim_adj(Xadj2, "Xadj2")
sce <- add_dim_adj(Xadj3, "Xadj3")

# calculate cms
sce <- cms(sce, group = batch, k = k, k_min = k_min, dim_red = "pca_Xadj1", res_name = "Xadj1")
sce <- cms(sce, group = batch, k = k, k_min = k_min, dim_red = "pca_Xadj2", res_name = "Xadj2")
sce <- cms(sce, group = batch, k = k, k_min = k_min, dim_red = "pca_Xadj3", res_name = "Xadj3")


### -------------- save sce object ----------------------###
saveRDS(sce, file = outputfile)


