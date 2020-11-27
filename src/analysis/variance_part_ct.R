##Variance Partitioning celltype specificity:
#Aim: Find out how much variance/gene can be associated to the batch effect

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
    library(variancePartition)
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

### ------------ Variance partitioning-------------------------###
#remove genes with zrero counts only
if( length(which(rowSums(assays(sce)$logcounts) == 0)) > 0 ){
    sce <- sce[-which(rowSums(assays(sce)$logcounts) == 0),]
}

#use hvg only 
dec <- modelGeneVar(sce)
chosen <- getTopHVGs(dec, prop=0.2)
sce <- sce[chosen,]

expr <- as.matrix(assays(sce)$logcounts)
meta_sub <- as.data.frame(colData(sce)[, c(celltype, batch)])
form <- as.formula(paste0("~ (1|", celltype, ") + (1|", batch, ") + (1|(", celltype, ":", batch, "))"))

varPart <- fitExtractVarPartModel(expr, form, meta_sub, BPPARAM=MulticoreParam(4))

# Add to sce
rowData(sce)$vp_celltype_ct <- varPart[[celltype]]
rowData(sce)$vp_batch_ct <- varPart[[batch]]
rowData(sce)$vp_cellt_by_batch_ct <- varPart[[paste0("(",celltype, ":", batch, ")")]]
rowData(sce)$vp_residuals_ct <- varPart[["Residuals"]]


### -------------- save sce object ----------------------###
saveRDS(sce, file = outputfile)


