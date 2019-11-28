##Simulation - Run edgeR:
# Aim: This is a preparation of the actual simulations input. 
# We want to simulate from real data and thus take the edgeR beta - coefficients 
# and logFoldchanges from real data to add sample effect/de and get back to count 
# data from there. 

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
    library(edgeR)
    library(scater)
    library(here)
    library(magrittr)
    library(dplyr)
    library(purrr)
})


# Read in data
sce <- readRDS(file = data)
params <- readRDS(file = params)

# vars
sample <- params[["sample"]]
celltype <- params[["celltype"]]

colData(sce)[,sample] <- as.factor(colData(sce)[,sample])
colData(sce)[,celltype] <- as.factor(colData(sce)[,celltype])

### ------------ EdgeR -------------------------###
#function to fit edgeR
fit_edgeR <- function(sample, celltype){
    be <- colData(sce)[,sample]
    ct <- colData(sce)[,celltype]
    levels(ct) <- paste("cell", 1:length(levels(ct)), sep = "_")
    # fit NB
    y <- DGEList(counts(sce))
    y <- calcNormFactors(y)
    mm <- model.matrix(~ 0 + be + ct)
    y <- estimateDisp(y, mm)
    y <- glmFit(y, prior.count = 0)

    # update row- & colData
    colData(sce)[,"offset"] <- c(y$offset)
    colData(sce)[,"lib.size"] <- y$samples$lib.size
    rowData(sce)[,"dispersion"] <- y$dispersion
    rowData(sce) <- DataFrame(rowData(sce), y$coefficients) %>% 
                                  set_colnames(c(colnames(rowData(sce)),
                                      paste("beta", 
                                            c(levels(colData(sce)[,sample]),
                                              levels(colData(sce)[,celltype])[-1]), 
                                            sep = ".")))
    sce
}

sce <- fit_edgeR(sample, celltype)

### -------------- save sce object ----------------------###
saveRDS(sce, file = outputfile)


