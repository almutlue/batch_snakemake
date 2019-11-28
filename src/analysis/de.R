# Differential expression:
#Aim: Batch characterization

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

set.seed(1000)

## Show arguments
print(data)
print(outputfile)
print(outputsce)
print(params)


#libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(limma)
  library(purrr)
  library(dplyr)
})


# Read in data
sce <- readRDS(file = data)
params <- readRDS(file = params)

# vars
batch <- params[["batch"]]
celltype <- params[["celltype"]]
sample <- params[["Sample"]]
cont <- params[["cont"]]

### ------------ Differential expression -------------------------###
# define parameter
clust <- as.factor(colData(sce)[,celltype])
kids <- levels(clust)
names(kids) <- kids
group <- as.factor(colData(sce)[,batch])
cs <- names(cont)
names(cs) <- cs
expr <- as.matrix(assays(sce)$logcounts)
ctype <- "contrast"

#prepare res dataframe
res_df <- function(k, tt, ct, c) {
  df <- data.frame(
    gene = rownames(tt), cluster_id = k, tt,
    row.names = NULL, stringsAsFactors = FALSE)
  df[[ct]] <- c
  return(df)
}

#run de
doDE <- function(sce, lfc_cutoff = 0){
  res <- lapply(kids, function (k) {
    cat(k, "..", sep = "")
    n <- clust == k
    es_tmp <- expr[, n]
    grp <- group[n]
    design <- model.matrix( ~ 0 + grp)
    colnames(design) <- levels(group)
    #k1 <- rowSums(es_tmp > 0) >= .2 * min(table(grp))
    #es_tmp <- es_tmp[k1, ]
    f <- lmFit(es_tmp, design)
    f <- eBayes(f, trend = TRUE)
    tt <- lapply(cont, function(c) {
      cc <- names(c)
      fc <- contrasts.fit(f, contrasts = c)
      tr <- treat(fc, lfc = lfc_cutoff)
      tt <- topTreat(tr, n = Inf)
      res_df(k, tt, ctype,cc)
    })
    return(list(tt = tt, data = es_tmp))
  })
  # remove empty clusters
  skipped <- vapply(res, is.null, logical(1))
  if( any(skipped) )
    message(paste("Cluster(s)", dQuote(kids[skipped]), "skipped due to an",
                  "insufficient number of cells in at least 2 samples per group."))
  res <- res[!skipped]
  kids <- kids[names(res)]

  # re-organize by contrast &
  # do global p-value adjustment
  tt <- lapply(res, "[[", "tt")
  tt <- lapply(cs, function(c) map(tt, c))

  # return results
  data <- lapply(res, "[[", "data")
  list(table = tt,
       data = data,
       design = design,
       coef = coef)
}
res <- doDE(sce,lfc_cutoff = 0)

##### ----------- Update rowData --------------- ######

combine_folds <- function(cont_var){
    #extract the contrast of interest and change log2fold colums names to be unique
    B <- res[["table"]][[cont_var]]
    new_name <- function(p){
        colnames(B[[p]])[3] <- paste0(cont_var, "_logFC_", p)
        return(B[[p]][,c(1,3)])
    }
    B_new_names <- lapply(names(B),new_name)
    names(B_new_names) <- names(B)
    #combine log2fold colums
    Folds <- Reduce(function(...){inner_join(..., by="gene")}, B_new_names)
}

all_folds <- lapply(cs, combine_folds)
rd <- Reduce(function(...){inner_join(..., by="gene")}, all_folds)
rd <- rd[match(rownames(sce), rd$gene),]
rd <- rd %>% mutate_all( ~replace(., is.na(.), 0))
rowData(sce)[,colnames(rd)] <- rd

### -------------- save output  ----------------------###
saveRDS(res, file = outputfile)
saveRDS(sce, file = outputsce)


