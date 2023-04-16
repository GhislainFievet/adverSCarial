#' Find a one gene min change adversarial attack list. A one gene min change
#' adversarial attack  refers to the modification of a single gene
#' within a cluster, leading to a change in its classification. The function
#' returns a list of genes/new classification.
#'
#' @details This function aims to get all genes that when modified individually
#' can lead to a misclassification. You can refer to the
#' 'advModifications' function documentation to more details on how
#' to define a modification.
#' The function is made as a two step parameter search. The first step is to
#' split the genes in 'firstDichot' sets, 100 by default. Then each set is
#' studied by a dichotomic process in a recursive function.
#' The aim of sarting by a hight value of sets, instead of strating directly by
#' the dichotomic research is to avoid the following scenario: we modify 5000
#' genes, the modification of one gene conpensates the modification of another.
#' The classification remains unchanged, whereas there is a one gene
#' classification modifying inside the 5000.
#' The dichotomic process runs as follow. The function receives a list of
#' genes, make the modification of the whole list and make the classification.
#' Three scenarios possible:
#' - the classification remains the same as the target cluster. The function
#'  returns, and the dichotomic process continues.
#' - the classification is changed. There is only one gene in the list, the
#'  function returns the gene and the new classification.
#' - the classification is changed. There is more than one gene in the list,
#'  the genes list is split in two, and the dichotomic process continues.
#' 
#' @param exprs a matrix, a data.frame or a DataFrame of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param clusters a character vector of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format
#' @param exclGenes a character vector of genes to exclude from the analysis
#' @param genes a character vector of genes in case you want to limit the
#' attack on a subset of genes
#' @param advMethod the name of the method to use
#' @param advFixedValue the numeric value to use in case of
#' advMethod=`fixed`
#' @param advFct the function to use in case advMethod
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param firstDichot the initial number of slices before
#' the dichotomic search
#' @param maxSplitSize max size of dichotomic slices
#' @param returnFirstFound set to TRUE to return result when a
#' the first misclassification is found
#' @param changeType `any` consider each misclassification,
#'  `not_na` consider each misclassification but NA.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a list of genes/new classification tuples
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'    c("T cell", 0.9)
#' }
#' rna_expression <- data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#' 
#' advMinChange(rna_expression, clusters_id,
#' "T cell", MyClassifier, advMethod="perc99")
#' @export
advMinChange <- function(exprs, clusters, target, classifier, exclGenes = c(),
        genes = c(), advMethod = "perc99", advFixedValue = 3,
        advFct = NULL, firstDichot = 100, maxSplitSize = 1,
        returnFirstFound = FALSE, changeType = "any", verbose = FALSE) {
    if ( !is(exprs, 'matrix') && !is(exprs,'data.frame') && !is(exprs,"DFrame")){
        stop("The argument exprs must be a matrix, a data.frame or a DataFrame.")
    }
    if (!is.character(clusters)) {
        stop("The argument clusters must be a vector of character.")
    }
    if (!is.character(target)) {
        stop("The argument target must be character.")
    }
    if (!is.function(classifier)){
        stop("The argument classifier must be a function.")
    }
    if (!is.character(exclGenes) && length(exclGenes)>0) {
        stop("The argument exclGenes must be character or vector of character.")
    }
    if (!is.character(genes) && length(genes)>0) {
        stop("The argument genes must be character or vector of character.")
    }
    if (!is.character(advMethod)) {
        stop("The argument advMethod must be character.")
    }
    if (!is.numeric(firstDichot)) {
        stop("The argument firstDichot must be numeric.")
    }
    if (!is.numeric(maxSplitSize)){
        stop("The argument maxSplitSize must be numeric.")
    }
    if (!is.logical(returnFirstFound)){
        stop("The argument returnFirstFound must be logical.")
    }
    if (!is.character(changeType)) {
        stop("The argument changeType must be a vector of character.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }
    if (is(exprs, "DFrame")){
        exprs <- as.data.frame(exprs)
    }

    genes_index <- seq_len(ncol(exprs))
    if (length(exclGenes) != 0 && length(genes) == 0) {
        genes_to_remove <- unlist(lapply(exclGenes, function(str_gene){
            if (!is.na(match(str_gene, colnames(exprs)))) {
                return(-match(str_gene, colnames(exprs)))
            } else {
                return(NULL)
            }
        }))
        if (length(genes_to_remove) != 0)
            genes_index <- genes_index[genes_to_remove]
    }
    if (length(genes) != 0) {
        exclGenes <- unique(exclGenes)
        genes <- unique(genes)
        genes_to_keep <- genes[!genes %in% exclGenes]
        ind_genes_to_keep <- unlist(lapply(genes_to_keep), function(str_gene){
            if (!is.na(match(str_gene, colnames(exprs)))) {
                return(match(str_gene, colnames(exprs)))
            } else {
                return(NULL)
            }
        })
        if (length(ind_genes_to_keep) != 0)
            genes_index <- genes_index[ind_genes_to_keep]
    }
    c_split_index <- split(genes_index, seq_len(firstDichot))
    l_splits_results <- c(); i <- 1
    while(i <= firstDichot) {
        prev_time <- Sys.time()
        message(paste0("Split number: ", i, "/", firstDichot))
        if (length(unlist(c_split_index[i])) != 0) {
            l_splits_results <- c( l_splits_results,
                .predictionDichotMinChange(exprs,
                    colnames(exprs)[unlist(c_split_index[i])],
                    clusters, target, classifier, advMethod = advMethod,
                    advFixedValue = advFixedValue, advFct = advFct,
                    maxSplitSize = maxSplitSize,
                    returnFirstFound = returnFirstFound,
                    changeType=changeType,verbose = verbose))
        }
        message(paste0("Split time: ", Sys.time() - prev_time))
        ifelse(length(l_splits_results) > 0 && returnFirstFound,
            i <- firstDichot + 1, i <- i + 1)
    }
    l_splits_results
}

.dichotMinSplit <- function(lResults, cGeneSplitValue, exprs,
            genes, clusters, target, classifier, advMethod, advFixedValue,
            advFct, maxSplitSize, changeType, returnFirstFound, verbose){
    # if ( length(lResults) > 0 && returnFirstFound){ return(lResults)}
    if (length(cGeneSplitValue) != 0) {
        if (verbose) {
            message("before predictWithNewValue")
            message(length(cGeneSplitValue))}
        mod_obj <- predictWithNewValue(exprs, cGeneSplitValue,
            clusters, target, classifier,
            advMethod = advMethod,
            advFixedValue = advFixedValue, advFct = advFct,
            verbose = TRUE)
        cell_type <- mod_obj[1]; celltype_score <- mod_obj[2]
        message(cell_type); message(celltype_score)
        if (cell_type != target) {
            if (length(cGeneSplitValue) <= maxSplitSize &&
                (changeType == "any" || cell_type != "NA" )) {
                message("Store gene:")
                message(cGeneSplitValue)
                lResults[[paste(cGeneSplitValue, collapse = "__")]] <-
                    c(cell_type, celltype_score)
            } else {
                if ( length(cGeneSplitValue) > maxSplitSize){
                    if (verbose) { message("Before predictionDichot:")
                        message(length(cGeneSplitValue))}
                    lResults <- .predictionDichotMinChange(exprs,
                        cGeneSplitValue, clusters, target, classifier,
                        lResults = lResults, advMethod = advMethod,
                        advFixedValue = advFixedValue, advFct = advFct,
                        maxSplitSize = maxSplitSize,
                        returnFirstFound=returnFirstFound,
                        changeType = changeType, verbose = verbose)
                }
            }
        }
    }
    lResults
}

.predictionDichotMinChange <- function(exprs, genes, clusters, target,
        classifier, lResults = list(), advMethod = "perc99",
        advFixedValue = 3, advFct = NULL, maxSplitSize = 1,
        returnFirstFound = FALSE, changeType = "any",
        verbose = TRUE) {
    if ( returnFirstFound && length(lResults) > 0) return (lResults)
    c_gene_split <- split(unlist(genes), c(1,2))
    message(paste0("genes size: ", length(unlist(genes))))
    lResults <- .dichotMinSplit(lResults,
                    unname(unlist(c_gene_split[1])), exprs,
                    genes, clusters, target, classifier, advMethod,
                    advFixedValue, advFct,
                    maxSplitSize, changeType,returnFirstFound,
                    verbose)
    if ( returnFirstFound && length(lResults) > 0) return (lResults)
    lResults <- .dichotMinSplit(lResults,
                    unname(unlist(c_gene_split[2])), exprs,
                    genes, clusters, target, classifier,
                    advMethod, advFixedValue,
                    advFct, maxSplitSize, changeType,returnFirstFound,
                    verbose)
    lResults
}

.maxOverListModifs <- function(exprs, clusters, classifier, exclGenes,
                            genes, modifications, maxSplitSize, verbose){
    df_result <- data.frame(todel = unique(clusters))
    rownames(df_result) <- unique(clusters)
    df_names <- c()
    for (modif_ind in seq_len(length(modifications))) {
        mod1 <- modifications[[modif_ind]][[1]]
        attacks_length <- c()
        for (cell_type in unique(clusters)) {
            if (verbose) {
                message(paste0("Running maxChange attack on ",
                    cell_type, ", with a maxSplitSize of: ", maxSplitSize))
                message(paste0("The smaller the maxSplitSize,",
                        " the more precise the result will be,",
                        " but it will take longer."))
                message("Modification: ",
                    paste(modifications[[modif_ind]], collapse = " "))
            }
            if (length(modifications[[modif_ind]]) == 1) {
                max_change_genes <- advMaxChange(exprs, clusters, cell_type,
                    classifier, advMethod = mod1,
                    maxSplitSize = maxSplitSize, exclGenes = exclGenes,
                    genes = genes, verbose = verbose)
            } else {
                mod2 <- modifications[[modif_ind]][[2]]
                max_change_genes <- advMaxChange(exprs, clusters,
                    cell_type, classifier, advMethod = mod1,
                    advFixedValue = mod2, advFct = mod2,
                    maxSplitSize = maxSplitSize, exclGenes = exclGenes,
                    genes = genes, verbose = verbose)
            }
            result_length <- length(max_change_genes)
            attacks_length <- c(attacks_length, result_length)
            if (verbose) {
                message(paste0("At least ", result_length,
                    " genes can be modified with the ",
                    paste(modifications[[modif_ind]], collapse = " "),
                    " method, and the cluster will still",
                    " be classified as ", cell_type))
            }
        }
        df_names <- c(df_names, paste(modifications[[modif_ind]],
            collapse = "_"))
        df_result <- cbind(df_result, attacks_length)
        df_result$todel <- NULL
        colnames(df_result) <- df_names
    }
    df_result
}

.maxOverArgModifs <- function(exprs, clusters, classifier, exclGenes,
                            genes, advMethod, advFixedValue,
                            advFct, maxSplitSize, verbose){
    attacks_length <- vapply(unique(clusters), function(cell_type){
        if (verbose) {
            message(paste0(
                "Running maxChange attack on ",
                cell_type, ", with a maxSplitSize of: ", maxSplitSize
            ))
            message(paste0("The smaller the maxSplitSize, the more",
                " precise the result will be, but it will take longer."))
        }
        max_change_genes <- advMaxChange(exprs, clusters, cell_type,
            classifier, advMethod = advMethod,
            advFixedValue = advFixedValue, advFct = advFct,
            maxSplitSize = maxSplitSize, exclGenes = exclGenes,
            genes = genes, verbose = verbose
        )
        result_length <- length(max_change_genes)
        if (verbose) {
            message(paste0(
                "At least ", result_length,
                " genes can be modified with the ", advMethod,
                " method, and the cluster will still be classified as ",
                cell_type
            ))
        }
        return(result_length)
    }, numeric(1))
    names(attacks_length) <- unique(clusters)
    df_result <- data.frame(at_least_gene_number = unlist(attacks_length))
    rownames(df_result) <- names(attacks_length)
    df_result
}


