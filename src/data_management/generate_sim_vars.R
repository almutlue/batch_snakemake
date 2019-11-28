##Generate combinations of simulation parameters :
# Call this using 'snakemake generate_sim_vars' 

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

set.seed(1000)

## Show arguments
print(out_path)

#libraries
suppressPackageStartupMessages({
    library(purrr)
})

#Parameter lists
rel_be <- list(1, 0.5, 2, 5, c(0.5, 1), c(0, 1))
rel_be_names <- c("1", "0.5", "2", "5", "0.5_1", "0_1")
rel_be_c <- list(1, 0.5, 2, c(0.1, 0.2, 0.3, 0.5, 0, 0.8, 1, 2))
rel_be_c_names <- c("1", "0.5", "2", "0.1_0.2_0.3_0.5_0_0.8_1_2")
l <- list("rel_be" = rel_be, "rel_be_c" = rel_be_c)
l_names <-list("rel_be" = rel_be_names, "rel_be_c "= rel_be_c_names)

#Combine parameter
sim <- cross(l)
names <- l_names %>% cross() %>% map(lift(paste)) %>% unlist
names <- gsub(" ", "__", names)
names(sim) <- names

#Add negative control
sim[["0__0"]] <- list("rel_be" = 0, "rel_be_c" = 0)

## write each element of sim to own R object
lapply(names(sim), function(sim_name){
    sub <- sim[[sim_name]]
    saveRDS(sub, file = paste0(out_path, sim_name, ".rds"))
})