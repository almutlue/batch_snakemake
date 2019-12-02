##Characterize Simulation - Reduced Dim representation :
# Aim: Add reduced dim representations to simulations

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

set.seed(1000)

## Show arguments
print(sce_sim)
print(outputfile)
print(sce_real)



#libraries
suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(countsimQC)
    library(magrittr)
    library(purrr)
    library(dplyr)
})

# Read in data
sample_name <- gsub("src/data/", "", sce_real)
sample_name <- gsub(".rds", "", sample_name)
sim_sce <- readRDS(file = sce_sim)
sce_real <- readRDS(file = sce_real)
out_path <- outputfile
out_path <- gsub("/", "", out_path)

out_file <- paste0("countSimQC_", sample_name, ".html")
paste(out_file)

### --------------Run CountSimQC------------###
count_list <- list("original" = as.matrix(counts(sce_real)), "sim" = as.matrix(counts(sim_sce)))

countsimQCReport(ddsList = count_list, outputFile = out_file,
                 outputDir = out_path, outputFormat = "html_document", 
                 showCode = FALSE, forceOverwrite = TRUE,
                 savePlots = FALSE, description = "Test sbatch simulation", 
                 maxNForCorr = 25, maxNForDisp = 100, 
                 calculateStatistics = TRUE, subsampleSize = 50,
                 kfrac = 0.01, kmin = 5, 
                 permutationPvalues = FALSE, nPermutations = NULL)

