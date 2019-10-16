args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

## Show arguments
print(data)
print(outputfile)
print(meta)

#libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    })


# Read in data
sce <- readRDS(file = data)
meta <- readRDS(file = meta)

# Define parameter
## general
n_genes <- nrow(sce)
n_cells <- ncol(sce)
dataset_name <- meta[["dataset_name"]]
batch <- meta[["batch"]]
celltype <- meta[["celltype"]]
patient <- meta[["patient"]]
sample <- meta[["sample"]]
## DE
cont <- meta[["cont"]]
## cms
k <- meta[["k"]]
k_min <- meta[["k_min"]]

params <- list("n_genes" = n_genes,
               "n_cells" = n_cells,
               "dataset_name" = dataset_name,
               "batch" = batch,
               "celltype" = celltype,
               "sample" = sample,
               "patient" = patient,
               "cont" = cont,
               "k" = k,
               "k_min" = k_min)

saveRDS(params, file = outputfile)


