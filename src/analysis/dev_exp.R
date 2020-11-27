## Deviance explained cell-type-specificity:
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
    library(edgeR)
})


# Read in data
sce <- readRDS(file = data)
params <- readRDS(file = params)

# vars
batch <- params[["batch"]]
celltype <- params[["celltype"]]
sample <- params[["Sample"]]

### ------------ Deviance explained ------------------------###
#remove genes with zrero counts only
if( length(which(rowSums(assays(sce)$logcounts) == 0)) > 0 ){
    sce <- sce[-which(rowSums(assays(sce)$logcounts) == 0),]
}

getDevianceExplained <- function(sce, form.full, form.null, tagwise=TRUE){
    if(is.null(sizeFactors(sce))) sce <- scran::computeSumFactors(sce)
    dds <- DGEList(as.matrix(counts(sce)))
    dds$samples$lib.size <- 1
    CD <- as.data.frame(colData(sce))
    CD$lsize <- log(sizeFactors(sce))
    mm <- model.matrix(form.full, data=CD)
    mm0 <- model.matrix(form.null, data=CD)
    dds <- estimateDisp(dds, mm, tagwise=tagwise)
    fit <- glmFit(dds, mm)
    fit0 <- glmFit(dds, mm0)
    de <- (deviance(fit0)-deviance(fit))/deviance(fit0)
    de[which(de<0)] <- 0
    return( de )
}

form_full <- as.formula(paste0("~ lsize + ", batch, ":", celltype, " + ", celltype, " + ", batch))
form_null <- as.formula(paste0("~ lsize +  ", celltype, " + ", batch))

de <- getDevianceExplained(sce, form.full = form_full, form.null = form_null)
rowData(sce)[,"de_exp_ct_be_int"] <- de

### -------------- save sce object ----------------------###
saveRDS(sce, file = outputfile)


