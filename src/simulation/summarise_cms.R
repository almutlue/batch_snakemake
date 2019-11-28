#Summarize all cms scores into colData of one object:

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

set.seed(1000)

## Show arguments
print(data)
print(outputfile)
print(params)

params_list = list(strsplit(params, ","))
print(params_list[[1]][[1]])


#libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(magrittr)
    library(dplyr)
    library(purrr)
    library(scater)
})


# Read in data
sce <- readRDS(file = data)

new_cms <- lapply(params_list[[1]][[1]], function(sim){
    name <- gsub('.*summarize_cms_sim_', "", outputfile)
    name <- gsub('_sce.rds', "", name)
    filename <- paste0("out/sim_char/", name, "/sim_", name, "_", sim, "_sce.rds")
    new_sce <- readRDS(file = paste0("out/sim_char/", name, "/sim_", name, "_", sim, "_sce.rds"))
    df <- data.frame(name = new_sce$cms)
}) %>% bind_cols() %>% set_colnames(paste0("cms_sim_", params_list[[1]][[1]]))

head(new_cms)

colData(sce) <- cbind(colData(sce), new_cms)
names(colData(sce))

# Add red_dim
new_dim <- lapply(params_list[[1]][[1]], function(sim){
    name <- gsub('.*summarize_cms_sim_', "", outputfile)
    name <- gsub('_sce.rds', "", name)
    filename <- paste0("out/sim_char/", name, "/sim_", name, "_", sim, "_sce.rds")
    new_sce <- readRDS(file = paste0("out/sim_char/", name, "/sim_", name, "_", sim, "_sce.rds"))
    redDim <- reducedDims(sce)[["UMAP"]]
}) %>% set_names(params_list[[1]][[1]])

reducedDims(sce) <- c(reducedDims(sce), new_dim)


### -------------- save sce object ----------------------###
saveRDS(sce, file = outputfile)


