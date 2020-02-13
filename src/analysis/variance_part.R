##Variance Partitioning:
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

expr <- as.matrix(assays(sce)$logcounts)
meta_sub <- as.data.frame(colData(sce)[, c(celltype, batch)])
form <- as.formula(paste0("~ (1|", celltype, ") + (1|", batch, ")"))
#meta_sub <- as.data.frame(colData(sce)[, c(batch)])
#colnames(meta_sub) <- batch
#form <- as.formula(paste0("~ (1|", batch, ")"))
varPart <- fitExtractVarPartModel(expr, form, meta_sub)
# Add to sce
rowData(sce)$vp_celltype <- varPart[[celltype]]
rowData(sce)$vp_batch <- varPart[[batch]]
rowData(sce)$vp_residuals <- varPart[["Residuals"]]


### -------------- save sce object ----------------------###
saveRDS(sce, file = outputfile)


