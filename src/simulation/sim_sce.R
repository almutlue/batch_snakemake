##Simulation :
# Aim: Simulate realistic batch effects
# We want to simulate from real data and thus take the edgeR beta - coefficients 
# and logFoldchanges from real data to add batch effect/de and get back to count 
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
print(sim)
print(out_path)


#libraries
suppressPackageStartupMessages({
    library(edgeR)
    library(scater)
    library(magrittr)
    library(dplyr)
    library(purrr)
    library(data.table)
})


# Read in data
sce <- readRDS(file = data)
params <- readRDS(file = params)
sim <- readRDS(file = sim)

# vars
batch <- params[["batch"]]
celltype <- params[["celltype"]]
sam <- params[["sample"]]

#add cluster_id and sample_id
cd <- data.frame(colData(sce))
cd <- dplyr::mutate_if(cd, is.factor, droplevels)
cd <- dplyr::mutate(cd, 
                    "sample_id" = as.factor(cd[, sam]), 
                    "cluster_id" = as.factor(cd[, celltype]), 
                    "batch_id" = as.factor(cd[, batch]))
colData(sce) <- S4Vectors::DataFrame(cd, row.names = colnames(sce))

colData(sce)[,batch] <- as.factor(colData(sce)[,batch])
colData(sce)[,celltype] <- as.factor(colData(sce)[,celltype])

#update rowData
colnames(rowData(sce)) <- gsub(paste0(".", batch), "", colnames(rowData(sce)))


#simulation input from summary
#dataset parameter
clust_prop <- table(sce$cluster_id)/sum(table(sce$cluster_id))
sample_prop <- table(sce$sample_id)/sum(table(sce$sample_id))
probs <- list("cluster" = table(sce$sample_id, sce$cluster_id)/ncol(sce), 
              "sample" = sample_prop, 
              "group" = c(0.5, 0.5))
n_ct = length(levels(sce$cluster_id))
n_be = length(levels(sce$batch_id))
rel_be <- rep(sim[["rel_be"]], n_be)
rel_be <- rel_be[1:n_be]
rel_be_c <- rep(sim[["rel_be_c"]], n_ct)
rel_be_c <- rel_be_c[1:n_ct]
rel_be
rel_be_c

lib_be <- "real"

if( all(rel_be == 0) ){
    lib_be = rep(1, n_be)
}

lib_be

## Create directory for outputfiles
dir.create(out_path, showWarnings = FALSE)

### ------------ Helper functions -------------------------###
#Predefined variables
cats <- c("ee", "ep", "de", "dp", "dm", "db")
cats <- factor(cats, levels = cats)

#Subsample cells for 2 groups
.sample_n_cells <- function(n, k, s) {
    nk <- length(k)
    ns <- length(s)
    if (length(n) == 1) {
        n <- list(rep(n, 2))
    } else {
        n <- replicate(nk * ns, 
                       list(sample(seq(n[1], n[2]), 2)))
    }
    matrix(n, 
           nrow = ns, ncol = nk, 
           dimnames = list(s, k))
}


#Sample metadata 
.sample_cell_md <- function(n, ids, probs = NULL) {
    ns <- vapply(ids, length, numeric(1))
    #no probs --> assign equally proportions to each celltype, sample, group
    if( is.null(probs) )                                 
        probs <- vector("list", 3)
    probs <- lapply(1:3, function(i) {
        if (!is.null(probs[[i]])) {
            return(probs[[i]])
        } else {
            rep(1 / ns[i], ns[i])
        }
    })
    #probs with celltype specific sample distribution 
    if( length(dim(probs[[1]])) == 2 ){
        if( is.null(names(probs[[2]])) ){
            names(probs[[2]]) <- ids[[2]]
        } 
        ct_s <- lapply(ids[[2]], function(s){
            ct_ids <- data.frame("cluster_id" = sample(ids[[1]], probs[[2]][s] * n, TRUE, probs[[1]][s,]), 
                                 "sample_id" = rep(s, probs[[2]][s] * n))
        }) %>% bind_rows()
        #Make sure that each celltype got at least 10 cells per sample
        if( any(table(ct_s) < 15) ){
            nk <- length(ids[[1]])
            be_tosmall <- ceiling(which(table(ct_s) < 10)/nk)
            ct_tosmall <- which(table(ct_s) < 10)/be_tosmall
            new_ct <- lapply(seq_len(length(be_tosmall)), function(rep){
                samp <- ids[[2]][be_tosmall[rep]]
                cellt <- ids[[1]][ct_tosmall[rep]]
                max_ct <- ids[[1]][table(ct_s[which(ct_s$sample_id %in% samp),]) %>% which.max()]
                start <- rep * 10
                finish <- start + 9
                new_ids <- which(ct_s$cluster_id %in% max_ct & ct_s$sample_id %in% samp)[start:finish]
                new_assign <- data.frame("ids" = new_ids, "new" = c(rep(cellt, 10)))
            }) %>% bind_rows()
            ct_s[new_ct$ids,"cluster_id"] <- new_ct$new
        }
        if( length(ct_s) < n ){
            ct_s[n,] <- c(ids[[1]][[1]], ids[[2]][[1]])
        } 
        ct_s$group_id <- sample(ids[[3]], n, probs[[3]], replace = TRUE)
        ct_s
    }else{
        #independent celltype, sample, group proportions
        vapply(1:3, function(i) 
            sample(ids[[i]], n, TRUE, probs[[i]]), 
            character(n)) %>% data.frame(row.names = NULL) %>% 
            set_colnames(c("cluster_id", "sample_id", "group_id"))
    }
}



#sample non-intersecting genes (typemarker?)
.sample_gene_inds <- function(gs, ns) {
    cluster_ids <- colnames(ns)
    vapply(cluster_ids, function(k)
        split(sample(gs), rep.int(cats, ns[, k])),
        vector("list", length(cats)))
}

### Sample log fold changes
# celltype lfcs
.sim_lfc_ct <- function(lfc_ct, p_type, gs){
    n <- length(gs)
    n_type <- n*p_type
    signs <- sample(c(-1, 1), n, TRUE, prob = c(0.5, 0.5))
    lfcs <- rtgamma(n, lfc_ct[1], lfc_ct[2], a = lfc_ct[3], b = lfc_ct[4]) * signs
    signs <- sample(c(-1, 1), n_type, TRUE, prob = c(0.5, 0.5))
    lfcs_type <- rnorm(n_type, mean = lfc_ct[5], sd = lfc_ct[6]) * signs
    gene_type <- sample(c(1:length(lfcs)), n_type)
    lfcs[gene_type] <- lfcs_type
    names(lfcs) <- gs
    return(lfcs)
}


# batch lfcs
.sim_lfc_be <- function(lfc_be, p_batch, gs){
    n <- length(gs)
    n_batch <- n * p_batch
    signs <- sample(c(-1, 1), n, TRUE)
    lfcs <- rtgamma(n, lfc_be[1], scale = lfc_be[2],
                    a = lfc_be[3],
                    b = lfc_be[4]) * signs
    signs <- sample(c(-1, 1), n_batch, TRUE)
    lfcs_batch <- rnorm(n_batch,
                        mean = lfc_be[5],
                        sd = lfc_be[6]) * signs
    gene_batch <- sample(c(1:length(lfcs)), n_batch)
    lfcs[gene_batch] <- lfcs_batch
    names(lfcs) <- gs
    return(lfcs)
}


#Get beta coefficients
.beta_coef <- function(x, equal, ids){
    n <- length(ids)
    beta_id <- data.frame("beta_ids" = paste("beta", ids, sep = ".")) %>%  
        filter(beta_ids %in% colnames(rowData(x)))
    rd <- as.data.frame(rowData(x)[, beta_id$beta_ids])
    rd <- cbind(rd, mean = apply(rd, 1, mean), 
                sd = abs(0.1 * apply(rd, 1, mean)))
    if( equal ){
        rd <- apply(rd, 1, function(gene){rnorm(n, gene["mean"], sd = 0)}) %>% 
            t() %>% set_colnames(ids)
    }else{
        rd <- apply(rd, 1, function(gene){
            rnorm(n = n, mean = gene["mean"], sd = gene["sd"])
        }) %>% t() %>% set_colnames(ids)
    }
}

#Impute type genes
.impute_type_genes <- function(x, gs_by_k, gs_idx, p_type) {
    kids <- colnames(gs_idx)
    names(kids) <- kids
    # sample gene-classes for genes of categroy EE & EP
    non_de <- c("ee", "ep")
    class <- lapply(kids, function(k) {
        gs <- unlist(gs_idx[non_de, k])
        n <- length(gs)
        data.table(
            stringsAsFactors = FALSE,
            gene = gs, cluster_id = k,
            class = sample(factor(c("state", "type")), n,
                           prob = c(1 - p_type, p_type), replace = TRUE))
    }) %>% map(split, by = "class", flatten = FALSE)
    # sample cluster-specific genes for ea. cluster & type-gene
    for (k in kids) {
        type_gs <- class[[k]]$type$gene
        #gs_by_k <- bind_rows(gs_by_k)
        gs_by_k[type_gs, k] <- apply(as.array(gs_by_k[type_gs, kids != k]), 
                                     1, function(ex) sample(setdiff(rownames(x), ex), 1))
    }
    return(gs_by_k)
}

# helper to sample from a NB across a grid of dispersions ('size') and means ('mu')
.nb <- function(cs, d, m, lfc = NULL) {
    n_gs <- length(d)
    n_cs <- length(cs)
    if (is.null(lfc))
        lfc <- rep(0, n_gs)
    lfc[lfc < 0] <- 0
    fc <- 2 ^ lfc
    fc <- rep(fc, each = n_cs)
    ds <- rep(1/d, each = n_cs)
    ms <- c(t(m[, cs])) * fc 
    rnbinom(n_gs * n_cs, size = ds, mu = ms) %>% 
        matrix(byrow = TRUE,
               nrow = n_gs, ncol = n_cs, 
               dimnames = list(names(d), cs)) %>% 
        list(counts = ., means = split(ms, rep(seq_len(nrow(m)), each = n_cs)))
}

# Simulate differential distributions
.sim <- function(
    cat = c("ee", "ep", "de", "dp", "dm", "db"),
    cs_g1, cs_g2, m_g1, m_g2, d, lfc) {
    
    cat <- match.arg(cat)
    ng1 <- length(cs_g1)
    ng2 <- length(cs_g2)
    
    re <- switch(cat,
                 ee = {
                     list(
                         .nb(cs_g1, d, m_g1),
                         .nb(cs_g2, d, m_g2))
                 },
                 ep = {
                     g1_hi <- sample(ng1, round(ng1 * 0.5))
                     g2_hi <- sample(ng2, round(ng2 * 0.5))
                     list(
                         .nb(cs_g1[-g1_hi], d, m_g1),
                         .nb(cs_g1[ g1_hi], d, m_g1, lfc), # 50% g1 hi
                         .nb(cs_g2[-g2_hi], d, m_g2),
                         .nb(cs_g2[ g2_hi], d, m_g2, lfc)) # 50% g2 hi
                 },
                 de = {
                     list(
                         .nb(cs_g1, d, m_g1, -lfc), # lfc < 0 => all g1 hi
                         .nb(cs_g2, d, m_g2,  lfc)) # lfc > 0 => all g2 hi
                 },
                 dp = {
                     props <- sample(c(0.3, 0.7), 2)
                     g1_hi <- sample(ng1, round(ng1 * props[1]))
                     g2_hi <- sample(ng2, round(ng2 * props[2]))
                     list(                           
                         .nb(cs_g1[-g1_hi], d, m_g1), 
                         .nb(cs_g1[ g1_hi], d, m_g1,  lfc), # lfc > 0 => 30/70% up
                         .nb(cs_g2[-g2_hi], d, m_g2), 
                         .nb(cs_g2[ g2_hi], d, m_g2, -lfc)) # lfc < 0 => 70/30% up
                 },
                 dm = {
                     g1_hi <- sample(ng1, round(ng1 * 0.5))
                     g2_hi <- sample(ng2, round(ng2 * 0.5))
                     list(
                         .nb(cs_g1[-g1_hi], d, m_g1),
                         .nb(cs_g1[ g1_hi], d, m_g1, -lfc), # lfc < 0 => 50% g1 hi
                         .nb(cs_g2[-g2_hi], d, m_g2),
                         .nb(cs_g2[ g2_hi], d, m_g2,  lfc)) # lfc > 0 => 50% g2 hi
                 }, 
                 db = {
                     g2_hi <- sample(ng2, round(ng2 * 0.5))
                     list(
                         .nb(cs_g1, d, m_g1, lfc/2),       # all g1 mi
                         .nb(cs_g2[-g2_hi], d, m_g2),      # 50% g2 lo
                         .nb(cs_g2[ g2_hi], d, m_g2, lfc)) # 50% g2 hi
                 }
    )
    cs <- map(re, "counts")
    cs <- do.call("cbind", cs)
    ms <- map(re, "means") %>%
        map_depth(2, mean) %>% 
        map_depth(1, unlist) %>% 
        bind_cols %>% as.matrix
    ms <- switch(cat, 
                 ee = ms,
                 de = ms,
                 db = cbind(
                     ms[, 1],
                     rowMeans(ms[, 2:3])),
                 cbind(
                     rowMeans(ms[, 1:2]),
                     rowMeans(ms[, 3:4]))) %>% 
        split(col(.)) %>% 
        set_names(c("A", "B"))
    list(cs = cs, ms = ms)
}

# Compute pseudo bulks
.pb <- function(x, cs, assay, fun) {
    fun <- switch(fun,
                  rowMedians = getFromNamespace(fun, "matrixStats"),
                  getFromNamespace(fun, "Matrix"))
    pb <- map_depth(cs, -1, function(i) {
        if (length(i) == 0) return(numeric(nrow(x)))
        fun(assays(x)[[assay]][, i, drop = FALSE])
    })
    map_depth(pb, -2, function(u) 
        data.frame(u, 
                   row.names = rownames(x),
                   check.names = FALSE))
}


#split cells
.split_cells <- function(x, 
                         by = c("cluster_id", "sample_id")) {
    if (is(x, "SingleCellExperiment"))
        x <- colData(x)
    cd <- data.frame(x[by], check.names = FALSE)
    cd <- data.table(cd, cell = rownames(cd)) %>% 
        split(by = by, sorted = TRUE, flatten = FALSE)
    map_depth(cd, length(by), "cell")
}

#Checks
.check_sce <- function(x, req_group = TRUE) {
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(c("cluster_id", "sample_id") %in% colnames(colData(x)))
    if (req_group)
        stopifnot("group_id" %in% colnames(colData(x)))
}


.check_arg_assay <- function(x, y) {
    stopifnot(is.character(y), length(y) == 1, y %in% assayNames(x))
    if (sum(assayNames(x) == y) > 1)
        stop("Argument 'assay' was matched to multiple times.\n ", 
             " Please assure that the input SCE has unique 'assayNames'.")
}

.check_input_param <- function(param, ids, value, th_up = Inf, th_l = 0) {
    if ( is.null(param) ) {
        param <- rep(value, length(ids))
        names(param) <- ids
    } else {
        stopifnot(is.numeric(param), length(param) == length(ids), param <= th_up, param >= th_l)
        if ( is.null(names(param)) ) {
            names(param) <- ids
        } else {
            stopifnot(setequal(names(param), ids))
        }
    }
    return(param)
}


.check_args_simData <- function(u) {
    stopifnot(
        is.numeric(u$n_genes), length(u$n_genes) == 1,
        is.numeric(u$n_cells), length(u$n_cells) == 1 | length(u$n_cells) == 2,
        is.numeric(u$p_dd), length(u$p_dd) == 6, sum(u$p_dd) == 1,
        is.numeric(u$p_type), length(u$p_type) == 1, u$p_type >= 0,
        is.numeric(u$lfc_dd), is.numeric(u$lfc_dd), u$lfc_dd >= 1,
        is.numeric(u$p_batch), length(u$p_batch) == 1, u$p_batch >= 0,
        is.numeric(u$lfc_be), length(u$lfc_be) == 6,
        is.numeric(u$lfc_ct), length(u$lfc_ct) == 6,
        u$be %in% c("real", "sim"), length(u$be) == 1,
        u$ct %in% c("est", "sim"), length(u$ct) == 1,
        all(u$lib_be %in% c("real") | is.numeric(u$lib_be)))
}


### -------------- function for simulation ----------------------###

simData <- function(x, n_genes = 500, n_cells = 300, 
                    probs = NULL, p_dd = diag(6)[1, ], p_type = 0,  p_batch = 0,
                    lfc_dd = 2, lfc_be = c(1, 0.2, -4, 4, 3, 0.5), lfc_ct = c(0.8, 0.5, -5, 5, 3.5, 0.9), 
                    rel_lfc = NULL, rel_be= NULL, p_ct = NULL, rel_be_c = NULL,
                    ct = "est", be = "real", lib_be = "real") {
    
    # throughout this code...
    # k: cluster ID
    # s: sample ID
    # c: gene category
    
    ########################## Start prepare input ##########################################
    
    .check_sce(x, req_group = FALSE)
    .check_args_simData(as.list(environment()))
    
    kids <- levels(x$cluster_id)
    sids <- levels(x$sample_id)
    bids <- levels(x$batch_id)
    gids <- c("A", "B")
    names(kids) <- kids
    names(sids) <- sids
    names(bids) <- bids
    names(gids) <- gids
    nk <- length(kids)
    ns <- length(sids)
    nbe <- length(bids)
    
    rel_lfc <- .check_input_param(rel_lfc, kids, value = 1)
    rel_be  <- .check_input_param(rel_be, bids, value = 1)
    rel_be_c <- .check_input_param(rel_be_c, kids, value = 1)
    p_ct <- .check_input_param(p_ct, kids, value = 0, th_up = 1)
    
    ########################## Finish prepare input ##########################################
    
    ########################## Start set up count matrix ###############################
    
    # initialize count matrix
    gs <- paste0("gene", seq_len(n_genes))
    cs <- paste0("cell", seq_len(n_cells))
    y <- matrix(0, n_genes, n_cells, dimnames = list(gs, cs))
    
    # sample cell metadata (see section Parameter - Detail "probs")
    cd <- .sample_cell_md(
        n = n_cells, probs = probs,
        ids = list(kids, sids, gids)) %>% set_rownames(cs)
    
    new_lab <- as_tibble(colData(x)) %>% group_by(sample_id, batch_id) %>% summarise()
    new_lab <- as.character(new_lab$batch_id) %>% set_names(new_lab$sample_id)
    cd <- cd %>% mutate("batch_id" = recode(as.factor(sample_id), !!!new_lab))
    cs_idx <- .split_cells(cd, by = colnames(cd)[-which(colnames(cd) %in% "batch_id")])
    n_cs <- modify_depth(cs_idx, -1, length)
    
    
    # split input cells by cluster-sample
    cs_by_ks <- .split_cells(x)
    
    # sample nb. of genes to simulate per category & gene indices
    n_dd <- replicate(nk, table(sample(factor(cats, levels = cats), n_genes, TRUE, p_dd))) %>% 
        set_colnames(kids)
    gs_idx <- .sample_gene_inds(gs, n_dd)
    
    # for ea. cluster, sample set of genes to simulate from
    gs_by_k <- setNames(sample(rownames(x), n_genes, FALSE), gs)
    gs_by_k <- replicate(nk, gs_by_k) %>% set_colnames(kids)
    
    
    # split by cluster & category
    gs_by_k <- gs_by_k %>% split(col(.)) %>% 
        set_names(kids) %>% map(set_names, gs)
    gs_by_kc <- lapply(kids, function(k) 
        lapply(cats, function(c)
            gs_by_k[[k]][gs_idx[[c, k]]]) %>% 
            set_names(cats))
    
    ########################## Finish set up count matrix ###############################
    
    ########################## Start sample lfc distributions #############################
    
    #logfold change of celltype
    if( ct %in% "sim" ){
        lfc_ct <- lapply(kids, function(k){
            .sim_lfc_ct(lfc_ct, p_type, gs)
        })
    }else{
        lfc_ct <- lapply(kids, function(k){
            lfcs <- rep(0, length(gs))
            names(lfcs) <- gs
            return(lfcs)
        })
    }
    
    
    #Sample celltype specific batch affected genes
    gs_dep_kc <- lapply(kids, function(k){
        lapply(cats, function(c) {
            n_dep <- n_dd[c, k] * p_ct[k]
            gs_dep <- sample(names(gs_by_kc[[k]][[c]]), n_dep, FALSE)
        }) %>% set_names(cats)
    }) 
    
    #logFC of batch effect
    ## Use celltype specific lfcs determined from the real data. They need to be provided 
    ## within the rowData slot and names strat with logFC
    if( be %in% "real" ){                      
        lfc_all <- names(rowData(x))[grepl("logFC", names(rowData(x)))]
        ref <- gsub('-.*', '', lfc_all[1])
        lfc_ref <- lfc_all[grepl(ref, lfc_all)]
        lfc_be <- lapply(bids, function(b){
            lfcb <- lfc_ref[grepl(b, lfc_ref)]
            lfcb_kc <- vapply(kids, function(k)
                lapply(cats, function(c){
                    if( b %in% ref ){
                        lfc_kc <- rep(0, length(gs)) %>% set_names(gs)
                    }else{
                        lfcbk_name <- lfcb[grepl(paste0('_', k, '$'), lfcb)]
                        lfcbk <- rowData(x)[, lfcbk_name] %>% set_names(gs)
                        lfc_kc <- lfcbk
                        lfc_kc <- lfc_kc * rel_be[[b]] * rel_be_c[[k]]
                    }
                    return(lfc_kc)
                }), vector("list", length(cats))) %>%
                set_rownames(cats)
            lfcb_kc
        })
    }
    
    ## logFold changes can be simulated by sampling from a gamma distribution
    ## Parameter need to be specified in lfc_be. Celltype specificity and 
    ## "batch genes" with more extreme values can also be specified"
    if( be %in% "sim" ){
        lfc_be <- lapply(bids, function(b){
            lfcb <- .sim_lfc_be(lfc_be, p_batch, gs)
            #celltype_specific logFC
            lfcb_kc <- vapply(kids, function(k)
                lapply(cats, function(c) {
                    lfc_kc <- lfcb[names(gs_by_kc[[k]][[c]])]
                    n_dep <- n_dd[c, k] * p_ct[k]
                    gs_dep <- sample(names(gs_by_kc[[k]][[c]]), n_dep, FALSE)
                    if( n_dep > 0 ){
                        lfcb_de <- .sim_lfc_be(lfc_be, p_batch, gs_dep)
                        names(lfcb_de) <- NULL
                    }else{
                        lfcb_de <- numeric(0)
                    }
                    lfc_kc[gs_dep_kc[[k]][[c]]] <- lfcb_de
                    lfc_kc <- lfc_kc * rel_be[[b]] * rel_be_c[[k]]
                    names(lfc_kc) <- names(gs_by_kc[[k]][[c]])
                    lfc_kc
                }), vector("list", length(cats))) %>%
                set_rownames(cats)
            lfcb_kc
        })
    }
    
    
    # sample logFCs group effects
    lfc_dd <- vapply(kids, function(k) 
        lapply(cats, function(c) { 
            n <- n_dd[c, k]
            if (c == "ee") return(rep(NA, n))
            signs <- sample(c(-1, 1), n, TRUE)
            lfcs <- rgamma(n, 4, 4/lfc_dd) * signs
            names(lfcs) <- gs_by_kc[[k]][[c]]
            lfcs * rel_lfc[k]
        }), vector("list", length(cats))) %>% 
        set_rownames(cats)
    lfc_dd <- setNames(replicate(ns, lfc_dd, simplify = FALSE), sids)
    
    
    ########################## Finish sample lfc distributions ###########################
    
    ########################## Start beta coefficients and lib sizes #######################
    
    #Get beta values for each sample
    #Are sample ids equal to batch_ids? (Add a sample specific variance or not)
    equal <- length(unique(names(new_lab))) == length(unique(new_lab))
    rd <- .beta_coef(x, equal = equal, sids)
    
    #Get beta values for celltypes
    if( ct %in% "sim" ){
        rd_ct <- as.data.frame(.beta_coef(x, equal = TRUE, kids)) %>% set_names(kids)
    }else{
        rd_ct <- rowData(x)[, paste("beta", kids, sep = ".")] %>% set_colnames(kids)
    }
    
    # compute NB parameters
    o <- exp(colData(x)$offset)
    
    if( lib_be %in% "real"){
        lib_batch <- tibble("batch_id" = bids,
                            "lib_be" = rep(1, length(bids)))
    }
    if( is.numeric(lib_be) ){
        lib_batch <- tibble("batch_id" = bids,
                            "lib_be" = lib_be)
        o <- sample(o, length(o), replace = FALSE)
    } 
    
    
    #libsizes * beta
    m <- lapply(sids, function(s) {
        b <- exp(rd[,s])
        be <- new_lab[[s]]
        o_be <- o * lib_batch$lib_be[lib_batch$batch_id %in% be]
        vapply(o_be, "*", b, FUN.VALUE = numeric(nrow(x))) %>% 
            set_rownames(rownames(x)) %>% 
            set_colnames(colnames(x))
    })
    
    d <- rowData(x)$dispersion %>% 
        set_names(rownames(x))
    
    
    ########################## Finish beta coefficients and lib sizes #######################
    
    ########################## Start simulate counts #######################################
    
    sim_mean <- lapply(sids, function(s) 
        lapply(kids, function(k)
            lapply(gids, function(g)
                setNames(numeric(n_genes), rownames(y)))))
    for (k in kids) {
        for (s in sids) {
            b <- new_lab[[s]]
            for (c in cats[n_dd[, k] != 0]) {
                gs_kc <- gs_by_kc[[k]][[c]]
                cs_ks <- cs_by_ks[[k]][[s]]
                
                g1 <- cs_idx[[k]][[s]]$A %>% as.numeric()
                g2 <- cs_idx[[k]][[s]]$B %>% as.numeric()
                
                ng1 <- length(g1)
                ng2 <- length(g2) 
                
                if( length(cs_ks) >= sum(ng1, ng2) ){
                    cs_g1 <- sample(cs_ks, ng1, replace = FALSE)
                    cs_g2_all <- cs_ks[-which(cs_ks %in% cs_g1)]
                    cs_g2 <- sample(cs_g2_all, ng2, replace = FALSE)
                }else{
                    cs_g1 <- sample(cs_ks, ng1, replace = TRUE)
                    cs_g2 <- sample(cs_ks, ng2, replace = TRUE)
                }
                
                #celltype specific beta
                b_ct_k <- exp(rd_ct[gs_kc, k]) 
                
                #Add  batch specific logfoldchange 
                lfc_be_s <- lfc_be[[b]][[c, k]][names(gs_kc)]
                lfc_be_s[is.na(lfc_be_s)] <- 0
                lfc_ct_k <- lfc_ct[[k]][names(gs_kc)]
                
                m_g1 <- m[[s]][gs_kc, cs_g1, drop = FALSE] * 2^lfc_be_s * 2^lfc_ct_k * b_ct_k
                m_g2 <- m[[s]][gs_kc, cs_g2, drop = FALSE] * 2^lfc_be_s * 2^lfc_ct_k * b_ct_k
                
                d_kc <- d[gs_kc]
                lfc_dd_s <- lfc_dd[[s]][[c, k]]
                
                gidx <- gs_idx[[c, k]]
                cidx <- c(g1, g2)
                
                re <- .sim(c, cs_g1, cs_g2, m_g1, m_g2, d_kc, lfc_dd_s)
                y[gidx, cidx] <- re$cs
                sim_mean[[s]][[k]]$A[gidx] <- re$ms$A
                sim_mean[[s]][[k]]$B[gidx] <- re$ms$B
            }
        }
    }
    
    
    ########################## Finish simulate counts #######################################
    
    ########################## Prepare results #############################################
    
    sim_mean <- lapply(sids, function(s){
        sim_mean[[s]] %>%
            map(bind_cols) %>% 
            bind_rows(.id = "cluster_id")  %>% 
            dplyr::mutate_at("cluster_id", factor) %>%      
            dplyr::mutate(gene = rep(gs, nk))})%>%
        reduce(inner_join, by = c("gene", "cluster_id"))
    
    colnames(sim_mean)[grepl("A",colnames(sim_mean))] <- paste0("A.",sids)
    colnames(sim_mean)[grepl("B",colnames(sim_mean))] <- paste0("B.",sids)
    
    
    #adjust lfc_be to summary table
    red_lfc_be <- lapply(bids, function(b){
        lfc_red <-lfc_be[[b]][1,] %>% unlist() 
    }) %>% bind_cols() %>% set_colnames(paste0("lfc_be_", colnames(.))) %>% 
        mutate("gene" = rep(gs, nk), "cluster_id" = rep(kids, each = length(gs)))
    
    # construct gene metadata table storing ------------------------------------
    # gene | cluster_id | category | logFC, gene, disp, mean used for sim.
    gi <- data.frame(
        gene = unlist(gs_idx),
        cluster_id = rep.int(rep(kids, each = length(cats)), c(n_dd)),
        category = rep.int(rep(cats, nk), c(n_dd)),
        sim_gene = unlist(gs_by_kc),
        sim_disp = d[unlist(gs_by_kc)],
        logFC_ct = unlist(lfc_ct)) %>% 
        dplyr::mutate_at("gene", as.character)
    # add true simulation means
    gi <- full_join(gi, sim_mean, by = c("gene", "cluster_id"))
    
    #add lfc_be
    gi <- full_join(gi, red_lfc_be, by = c("gene", "cluster_id"))
    
    # reorder
    o <- order(as.numeric(gsub("[a-z]", "", gi$gene)))
    gi <- gi[o, ] %>% set_rownames(NULL)
    
    # construct SCE
    cd$sample_id <- factor(paste(cd$sample_id, cd$group_id, sep = "."))
    m <- match(levels(cd$sample_id), cd$sample_id)
    gids <- cd$group_id[m]
    o <- order(gids)
    sids <- levels(cd$sample_id)[o]
    ei <- data.frame(sample_id = sids, group_id = gids[o])
    cd <- cd %>% dplyr::mutate_at("sample_id", factor, levels = sids)
    
    md <- list(
        experiment_info = ei,
        n_cells = table(cd$sample_id),
        gene_info = gi)
    
    SingleCellExperiment(
        assays = list(counts = as.matrix(y)),
        colData = cd, 
        metadata = md)
}


### -------------- simulate batch effects ----------------###

#Simulation
sim_batch <- simData(sce, n_genes = nrow(sce), n_cells = ncol(sce), 
                     probs = probs, be = "real", ct = "est", 
                     rel_be_c = rel_be_c, rel_be = rel_be, lib_be = lib_be)


### -------------- save sce object ----------------------###
saveRDS(sim_batch, file = outputfile)


