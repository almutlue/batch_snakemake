#Generate Params file to use characterization scripts:


args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}


## Show arguments
print(outputfile)
print(params)


#libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
})


# Read in data
params <- readRDS(file = params)

# Change params to fit simulation
params[["batch"]] <- "batch_id"
params[["celltype"]] <- "cluster_id"
params[["Sample"]] <- "sample_id"

#save params
saveRDS(params, file = outputfile)
