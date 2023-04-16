#' Returns a modified RNA expression matrix, for a given cluster,
#' for a given modification.
#'
#' @details The motivation for this function is to standardize the modifications
#' we want to study in the attacks. We give as argument a matrix of the RNA
#' expression, the gene and the target cells we want to modify. Then we have
#' three arguments allowing to specify what modification we want to apply on
#' these cells. The advMethod contains, a specific prerecorded modification
#' or an indication on how to use the other two arguments.
#' The prerecorded modifications available for the advMethod argument are:
#' - 'perc1', replace the value by the whole matrix 1 percentile value of
#'  the gene. It is as if we biologically switched off the gene.
#' - 'perc99', replace the value by the whole matrix 99 percentile value of
#'  the gene. It is as if we biologically switched on the gene to the maximum.
#' The value of the advMethod argument can also be 'fixed', in this case the
#' modification would be to replace the value of the gene of the wanted cells
#' by the value of the argument 'advFixedValue'. This can be useful to test
#' aberrant values like negative integer, absurdly high values of character
#' values.
#' The value of the advMethod argument can also be 'full_row_fct',
#' 'target_row_fct', 'target_matrix_fct' or 'full_matrix_fct'. They are
#' used when we want to use a custom modification function, with the 'advFct'
#' argument:
#' - 'full_row_fct' indicate that the 'advFct' function takes the whole gene
#'  values as input.
#' - 'target_row_fct' indicate that the 'advFct' function takes target cells
#'  gene values as input.
#' - 'full_matrix_fct' indicate that the 'advFct' function takes the whole gene
#'  expression values as input.
#' - 'target_matrix_fct' indicate that the 'advFct' function takes target cells
#'  all genes values as input.
#' 
#' @param exprs a matrix, a data.frame or a DataFrame of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param genes the character vector of genes to modify
#' @param clusters a character vector of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param advMethod the name of the method to use
#' @param advFixedValue the numeric value to use in case of
#' advMethod=`fixed`
#' @param advFct the function to use in case advMethod
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return the matrix or a dataframe exprs modified on asked genes
#' with the specified modification
#' @examples
#' rna_expression <- data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#'
#' advModifications(rna_expression, genes, clusters_id,
#' "T cell", advMethod="perc99")
#' @export
advModifications <- function(exprs, genes, clusters,
                            target, advMethod = "perc99",
                            advFixedValue = 3, advFct = NULL,
                            verbose = FALSE) {
    if ( !is.matrix(exprs) && !is.data.frame(exprs) && !is(exprs,"DFrame")){
        stop("The argument exprs must be a matrix, a data.frame or a DataFrame.")
    }
    if (!is.character(genes)) {
        stop("The argument genes must be character or vector of character.")
    }
    if (!is.character(clusters)) {
        stop("The argument clusters must be a vector of character.")
    }
    if (!is.character(target)) {
        stop("The argument target must be character.")
    }
    if (!is.character(advMethod)) {
        stop("The argument advMethod must be character.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }

    if (verbose) {
        message("Modify data for ", length(genes),
            " genes for cluster ", target)
    }
    if (is(exprs, "DFrame")){
        exprs <- as.data.frame(exprs)
    }

    if (!is.function(advFct)) {
        exprs <- .advModificationsNotFunction(exprs,
                                    genes, clusters,
                                    target, advMethod = advMethod,
                                    advFixedValue = advFixedValue)
    } else {
        exprs <- .advModificationsFunction(exprs, genes, clusters,
                            target, advMethod = advMethod,
                            advFct = advFct)
    }

    exprs
}

.advModificationsNotFunction <- function(exprs,
                            genes, clusters,
                            target, advMethod = "perc99",
                            advFixedValue = 3) {
        cell_mask <- clusters == target
        if (advMethod == "fixed") {
            num_genes <- c()
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    num_genes <- append(num_genes, my_gene)
                }
            }
            exprs[cell_mask, num_genes] <- advFixedValue
        }
        if (advMethod == "perc99") {
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    exprs[cell_mask, my_gene] <- stats::quantile(
                        exprs[, my_gene],
                        0.99
                    )
                }
            }
        }
        if (advMethod == "perc1") {
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    exprs[cell_mask, my_gene] <- stats::quantile(
                        exprs[, my_gene],
                        0.01
                    )
                }
            }
        }
        exprs
}

.advModificationsFunction <- function(exprs, genes, clusters,
                            target, advMethod = "perc99",
                            advFct = NULL) {
        cell_mask <- clusters == target
        if (advMethod == "full_row_fct" || advMethod == "fixed") {
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    exprs[cell_mask, my_gene] <- advFct(exprs[, my_gene])
                }
            }
        }
        if (advMethod == "target_row_fct" || advMethod == "fixed") {
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    exprs[cell_mask, my_gene] <- advFct(exprs[
                        cell_mask,
                        my_gene
                    ])
                }
            }
        }
        if (advMethod == "target_matrix_fct") {
            num_genes <- unlist(lapply(unique(genes),function(my_gene){
                if (is(exprs[, my_gene], "numeric")) {
                    return(my_gene)
                } else {
                    return(NULL)
                }
            }))
            exprs[cell_mask, num_genes] <- advFct(exprs[cell_mask, num_genes])
        }
        if (advMethod == "full_matrix_fct") {
            num_genes <- unlist(lapply(unique(genes),function(my_gene){
                if (is(exprs[, my_gene], "numeric")) {
                    return(my_gene)
                } else {
                    return(NULL)
                }
            }))
            exprs[cell_mask, num_genes] <- advFct(exprs[, num_genes])
        }
    exprs
}

#' Returns a classification and an odd value from a
#' RNA expression matrix, for given genes, for a given cluster,
#' for a given modification.
#'
#' @details This function aims to concatenate the following actions:
#' - modify the RNA gene expression
#' - classify the result
#' This is a function widely used in the other functions of the package.
#' @param exprs a matrix, a data.frame or a DataFrame of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param genes the character vector of genes to modify
#' @param clusters a character vector of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format
#' @param advMethod the name of the method to use
#' @param advFixedValue the numeric value to use in case of
#' advMethod=`fixed`
#' @param advFct the function to use in case advMethod
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a vector of the classification, and the associated odd
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'    c("T cell", 0.9)
#' }
#' rna_expression <- data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#'
#' predictWithNewValue(rna_expression, genes, clusters_id,
#' "T cell", MyClassifier, advMethod="perc99")
#' @export
predictWithNewValue <- function(exprs, genes, clusters, target,
                                classifier, advMethod = "perc99",
                                advFixedValue = 3,
                                advFct = NULL, verbose = FALSE) {
    if ( !is.matrix(exprs) && !is.data.frame(exprs) && !is(exprs,"DFrame")){
        stop("The argument exprs must be a matrix, a data.frame or a DataFrame.")
    }
    if (!is.character(genes)) {
        stop("The argument genes must be character or vector of character.")
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
    if (!is.character(advMethod)) {
        stop("The argument advMethod must be character.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }

    if (verbose) {
        message("Modify data for ", length(genes),
            " genes for cluster ", target)
    }
    if (is(exprs, "DFrame")){
        exprs <- as.data.frame(exprs)
    }

    modif_exprs <- advModifications(exprs, genes, clusters,
        target,
        advMethod = advMethod,
        advFixedValue = advFixedValue,
        advFct = advFct
    )

    classifier(modif_exprs, clusters, target)
}

#' Find a max change adversarial attack. It finds the longer
#' list of genes you can modify on a cluster without changing its
#' classification.
#' 
#' @details This function aims to get the largest part of the genes
#' that can be modified without altering the classification,
#' considering a given modification. You can refer to the
#' 'advModifications' function documentation to more details on how
#' to define a modification.
#' The search is made by a dichotomic process, on a reccursive function.
#' At each iteration the function splits the genes in two groups. It
#' proceeds to the modification of the RNA gene value of the first group,
#' makes its classification. Then three possible scenarios.
#' - the classification is the same as the target cluster. We concat
#'  the genes list to the previous one, make the classification, and it still
#'  gives same classification. Then we return the genes list.
#' - the classification is the same as the target cluster. We concat
#'  the genes list to the previous one, make the classification, and it
#'  gives a different classification. This happens often, you can modify
#'  the gene A with a classification of T cell, or modify the gene B with a
#'  classification of T cell, but modifying A and B returns another
#'  classification. In this case we split the genes list in two and try again.
#' - the classification is not the same as the target cluster. In this case we
#'  split the genes list in two and try again.
#' The iteration process stops when the length of the genes list is lower
#' than the value of the 'maxSplitSize' argument. So you should set it to 1
#' to have the maximum number of genes for the max change attack. This function
#' is used by the 'overMaxChange' function with a default argument value of 100
#' to increase speed, and still returns significant results.
#' @param exprs a matrix, a data.frame or a DataFrame of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param clusters a character vector of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format
#' @param exclGenes a list of genes to exclude from the analysis
#' @param genes a list of genes in case you want to limit the
#' attack on a subset of genes
#' @param advMethod the name of the method to use
#' @param advFixedValue the numeric value to use in case of
#' advMethod=`fixed`
#' @param advFct the function to use in case advMethod
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param maxSplitSize max size of dichotomic slices.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a character vector of genes you can modify on a cluster without
#' modifying its classification
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'    c("T cell", 0.9)
#' }
#' rna_expression <- data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#' 
#' advMaxChange(rna_expression, clusters_id,
#' "T cell", MyClassifier, advMethod="perc99")
#' @export
advMaxChange <- function(exprs, clusters, target, classifier,
                        exclGenes = c(), genes = c(), advMethod = "perc99",
                        advFixedValue = 3, advFct = NULL,
                        maxSplitSize = 1, verbose = FALSE) {
    if ( !is.matrix(exprs) && !is.data.frame(exprs) && !is(exprs,"DFrame")){
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
    if (!is.numeric(maxSplitSize)){
        stop("The argument maxSplitSize must be numeric.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }
    if (is(exprs, "DFrame")){
        exprs <- as.data.frame(exprs)
    }

    if (maxSplitSize < 1){ maxSplitSize <- 1 }
    genes_to_keep <- unlist(lapply(unique(genes), function(str_gene){
        if (!is.na(match(str_gene, colnames(exprs)))) {
            return(str_gene)
        } else {
            return(NULL)
        }
    }))
    if (length(genes) == 0) {
        genes_to_keep <- colnames(exprs)
    }
    for (str_gene in unique(exclGenes)) {
        if (!is.na(match(str_gene, genes_to_keep))) {
            genes_to_keep <- genes_to_keep[-match(str_gene, genes_to_keep)]
        }
    }
    lResults <- .predictionDichotMaxChange(exprs, genes_to_keep,
        clusters, target, classifier,
        advMethod = advMethod,
        advFixedValue = advFixedValue, advFct = advFct,
        maxSplitSize = maxSplitSize, verbose = verbose
    )
    lResults
}

.dichotMaxSameType <- function(lResults, cGeneSplitValue, exprs,
                clusters, target, classifier, advMethod,
                advFixedValue, advFct, maxSplitSize, verbose = TRUE){
        if (verbose) { message("same cell_type") }
        if (length(lResults) == 0) {
            lResults <- cGeneSplitValue
        } else {
            if (verbose){
                message("check if concat lists still gives target")
            }
            concat_genes_list <- c(lResults, cGeneSplitValue)
            if (verbose) {
                message(paste0(
                    "length of tested gene list: ",
                    length(unique(concat_genes_list))
                ))
            }
            # check if concat lists still gives target
            mod_obj <- predictWithNewValue(exprs, concat_genes_list, clusters,
                target, classifier,
                advMethod = advMethod,
                advFixedValue = advFixedValue, advFct = advFct,
                verbose = verbose
            )
            concat_cell_type <- mod_obj[1]
            if (concat_cell_type == target) {
                if (verbose) { message("YES: merge results") }
                # Merge results
                lResults <- c(lResults, cGeneSplitValue)
                if (verbose) { message(paste0( "results length after merge: ",
                        length(lResults)
                    ))
                }
            } else {
                if (verbose) { message("NO: split and retry") }
                if (length(cGeneSplitValue) > maxSplitSize) {
                    lResults <- .predictionDichotMaxChange(exprs,
                        cGeneSplitValue, clusters, target,
                        classifier,
                        lResults = lResults,
                        advMethod = advMethod,
                        advFixedValue = advFixedValue, advFct = advFct,
                        maxSplitSize = maxSplitSize, verbose = verbose
                    )
                }
            }
        }
    lResults
}

.dichotMaxSplit <- function(lResults, cGeneSplitValue, exprs,
                    genes, clusters, target,
                                    classifier,
                                    advMethod,
                                    advFixedValue,
                                    advFct, maxSplitSize,
                                    verbose){
    if (length(cGeneSplitValue) != 0) {
        if (verbose) {
            message("before predictWithNewValue 1")
            message(length(cGeneSplitValue))
        }
        mod_obj <- predictWithNewValue(exprs, cGeneSplitValue,
            clusters, target, classifier,
            advMethod = advMethod,
            advFixedValue = advFixedValue, advFct = advFct,
            verbose = TRUE
        )
        cell_type <- mod_obj[1]
        celltype_score <- mod_obj[2]
        message(cell_type)
        message(celltype_score)
    }
    if (cell_type == target) {
        lResults <- .dichotMaxSameType(lResults,
                    cGeneSplitValue,
                    exprs, clusters,
                    target, classifier,
                    advMethod,
                    advFixedValue, advFct,
                    maxSplitSize,
                    verbose = TRUE)
    } else {
        if (verbose) { message("NOT same cell_type") }
        if (length(cGeneSplitValue) > maxSplitSize) {
            lResults <- .predictionDichotMaxChange(exprs,
                unname(cGeneSplitValue), clusters, target,
                classifier, lResults = lResults, advMethod = advMethod,
                advFixedValue = advFixedValue, advFct = advFct,
                maxSplitSize = maxSplitSize, verbose = verbose
            )
        }
    }
    lResults
}


.predictionDichotMaxChange <- function(exprs, genes, clusters, target,
                                    classifier, lResults = c(),
                                    advMethod = "perc99",
                                    advFixedValue = 3,
                                    advFct = NULL, maxSplitSize = 1,
                                    verbose = TRUE) {
    c_gene_split <- split(unlist(genes), c(1,2))
    message(paste0("genes size: ", length(unlist(genes))))
    message(paste0("current gene results length: ", length(lResults)))
    # Random order for each reccursion
    random_index <- sample(c(1,2), 2)
    c_gene_split <- c_gene_split[random_index]
    lResults <- .dichotMaxSplit(lResults, unname(unlist(c_gene_split[1])),
                    exprs, genes, clusters, target, classifier, advMethod,
                    advFixedValue, advFct, maxSplitSize, verbose)
    lResults <- .dichotMaxSplit(lResults, unname(unlist(c_gene_split[2])),
                exprs, genes, clusters, target, classifier, advMethod,
                advFixedValue, advFct, maxSplitSize, verbose)
    lResults
}

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
    if ( !is.matrix(exprs) && !is.data.frame(exprs) && !is(exprs,"DFrame")){
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


#' Gives an overview of the susceptibility to max change
#' attacks, for each cell type, for a given list of modifications.
#'
#' @details Running the advMaxChange function for each cell type
#' to see which ones are more vulerable can take a long time. The
#' aim of the maxChangeOverview function is to make this process faster.
#' It uses a default value of 100 for the 'maxSplitSize' parameter. So,
#' the dichotomic process of the advMaxChange function stops as soon
#' as the fold length is lower than 100. You can have more accurate
#' results with maxSplitSize=1, but it will take longer.
#' This function aims also to run the advMaxChange for several given
#' modifications. You can specify a list of modifications as so - each
#' item of the list should be 1 or 2 length size.
#' The 1 length vector must contain the prerecorded modifications, 'perc1'
#' or 'perc99'.
#' The 2 length vector must have as first item:
#'  - 'fixed', in this case the second item should be the value to be
#'  replaced by.
#'  - 'full_row_fct', 'target_row_fct', 'target_matrix_fct' or
#'  'full_matrix_fct'. In this case the second item should be a function.
#' Let's say we want to analysis the susceptibility to max change attack
#' for 3 modifications: "perc1", the modification of each value of the
#' cluster by 1000, and a custom modification stored inside a function myFct.
#' Then the 'modification' parameter should be:
#' my_modifications = list(c("perc1"),
#'                         c("fixed", 1000),
#'                         c("full_matrix_fct", myFct))
#' 
#' The function returns a dataframe with the number of genes of the max
#' change attack for each modification in columns, for each cell type in rows.
#' 
#' @param exprs a matrix, a data.frame or a DataFrame of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param clusters a character vector of the clusters to which the cells belong
#' @param classifier a classifier in the suitable format
#' @param exclGenes a character vector of genes to exclude from the analysis
#' @param genes a character vector of genes in case you want to limit the
#' analysis on a subset of genes
#' @param modifications the list of the modifications to study
#' @param advMethod the name of the method to use
#' @param advFixedValue the numeric value to use in case of
#' advMethod=`fixed`
#' @param advFct the function to use in case advMethod
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param maxSplitSize max size of dichotomic slices.
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
#' maxChangeOverview(rna_expression, clusters_id,
#' MyClassifier, modifications = list(c("perc1"), c("perc99")))
#' 
#' myModif = function(x){
#'    return(sample(1:10,1))
#' }
#' 
#' my_modifications = list(c("perc1"),
#'                         c("fixed", 1000),
#'                         c("full_matrix_fct", myModif))
#' maxChangeOverview(rna_expression, clusters_id,
#'  MyClassifier, modifications = my_modifications)
#' @export
maxChangeOverview <- function(exprs, clusters, classifier, exclGenes = c(),
                            genes = c(),
                            modifications = list(c("perc1"), c("perc99")),
                            advMethod = "perc99", advFixedValue = 3,
                            advFct = NULL, maxSplitSize = 100,
                            verbose = FALSE) {
    if ( !is.matrix(exprs) && !is.data.frame(exprs) && !is(exprs,"DFrame")){
        stop("The argument exprs must be a matrix, a data.frame or a DataFrame.")
    }
    if (!is.character(clusters)) {
        stop("The argument clusters must be a vector of character.")
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
    if (!is.list(modifications)) {
        stop("The argument modifications must be a list.")
    }
    if (!is.character(advMethod)) {
        stop("The argument advMethod must be character.")
    }
    if (!is.numeric(maxSplitSize)){
        stop("The argument maxSplitSize must be numeric.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }
    if (is(exprs, "DFrame")){
        exprs <- as.data.frame(exprs)
    }

    if (length(modifications) == 0) {
        df_result <- .maxOverArgModifs(exprs, clusters, classifier,
            exclGenes, genes, advMethod, advFixedValue,
            advFct, maxSplitSize, verbose)
    } else {
        df_result <- .maxOverListModifs(exprs, clusters, classifier, exclGenes,
            genes, modifications, maxSplitSize, verbose)
    }
    df_result
}


.minOverListModifs <- function(exprs, clusters, classifier, exclGenes,
                            genes, modifications, firstDichot, maxSplitSize,
                            changeType, verbose){
    df_result <- data.frame(todel = unique(clusters))
    rownames(df_result) <- unique(clusters); df_names <- c()
    for (modif_ind in seq_len(length(modifications))) {
        mod1 <- modifications[[modif_ind]][[1]]
        attacks_length <- vapply(unique(clusters), function(cell_type){
            if (verbose) {
                message("Running minChange attack on ", cell_type,
                    ", with a maxSplitSize of: ", maxSplitSize)
                message("The smaller the maxSplitSize, the more precise",
                    " the result will be, but it will take longer.")
                message("Modification: ",
                    paste(modifications[[modif_ind]], collapse = " "))
            }
            if (length(modifications[[modif_ind]]) == 1) {
                min_change_genes <- advMinChange(exprs, clusters, cell_type,
                    classifier, exclGenes = exclGenes, genes = genes,
                    advMethod = mod1, maxSplitSize = maxSplitSize,
                    firstDichot = firstDichot, changeType = changeType,
                    verbose = verbose)
            } else {
                mod2 <- modifications[[modif_ind]][[2]]
                min_change_genes <- advMinChange(exprs, clusters, cell_type,
                    classifier, exclGenes = exclGenes, genes = genes,
                    advMethod = mod1, advFixedValue = mod2, advFct = mod2,
                    maxSplitSize = maxSplitSize,
                    firstDichot = firstDichot, changeType = changeType,
                    verbose = verbose)
            }
            result_length <- length(min_change_genes)
            if (verbose) {
                message("An approximation gives about ", result_length,
                    " genes can cause a one gene min change attack on the ",
                    cell_type, " cell type for the modification ",
                    paste(modifications[[modif_ind]], collapse = " "))
            }
            return(result_length)
        }, numeric(1))
        df_names <- c(df_names,
            paste(modifications[[modif_ind]], collapse = "_"))
        df_result <- cbind(df_result, attacks_length)
        df_result$todel <- NULL
        colnames(df_result) <- df_names
    }
    df_result
}


.minOverArgModifs <- function(exprs, clusters, classifier, exclGenes,
                    genes, advMethod, advFixedValue, advFct,
                    firstDichot, maxSplitSize, changeType, verbose) {
    attacks_length <- vapply(unique(clusters), function(cell_type){
        if (verbose) {
            message(paste0(
                "Running minChange attack on ", cell_type,
                ", with a maxSplitSize of: ", maxSplitSize
            ))
            message(
                "The smaller the maxSplitSize, the more",
                " precise the result will be, but it will take longer."
            )
        }
        min_change_genes <- advMinChange(exprs, clusters, cell_type,
            classifier,
            exclGenes = exclGenes,
            genes = genes, advMethod = advMethod,
            advFixedValue = advFixedValue,
            advFct = advFct,
            maxSplitSize = maxSplitSize,
            firstDichot = firstDichot, changeType = changeType,
            verbose = verbose
        )
        result_length <- length(min_change_genes)
        if (verbose) {
            message(paste0(
                "An approximation gives about ",
                result_length,
                " genes can cause a one gene min change attack on the ",
                cell_type, " cell type"
            ))
        }
        return(result_length)
    }, numeric(1))
    names(attacks_length) <- unique(clusters)
    df_result <- data.frame(about_gene_number = unlist(attacks_length))
    rownames(df_result) <- names(attacks_length)
    df_result
}

#' Gives an overview of the susceptibility to min change
#' attacks, for each cell type, for a given list of modifications.
#'
#' @details Running the advMinChange function for each cell type
#' to see which ones are more vulerable can take a long time. The
#' aim of the minChangeOverview function is to make this process faster.
#' It uses a default value of 100 for the 'maxSplitSize' parameter. So,
#' the dichotomic process of the advMinChange function stops as soon
#' as the fold length is lower than 100. You can have more accurate
#' results with maxSplitSize=1, but it will take longer.
#' This function aims also to run the advMinChange for several given
#' modifications. You can specify a list of modifications as so - each
#' item of the list should be 1 or 2 length size.
#' The 1 length vector must contain the prerecorded modifications, 'perc1'
#' or 'perc99'.
#' The 2 length vector must have as first item:
#'  - 'fixed', in this case the second item should be the value to be
#'  replaced by.
#'  - 'full_row_fct', 'target_row_fct', 'target_matrix_fct' or
#'  'full_matrix_fct'. In this case the second item should be a function.
#' Let's say we want to analysis the susceptibility to min change attack
#' for 3 modifications: "perc1", the modification of each value of the
#' cluster by 1000, and a custom modification stored inside a function myFct.
#' Then the 'modification' parameter should be:
#' my_modifications = list(c("perc1"),
#'                         c("fixed", 1000),
#'                         c("full_matrix_fct", myFct))
#' 
#' The function returns a dataframe with the number of genes of the max
#' change attack for each modification in columns, for each cell type in rows.
#' @param exprs a matrix, a data.frame or a DataFrame of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param clusters a character vector of the clusters to which the cells belong
#' @param classifier a classifier in the suitable format
#' @param exclGenes a character vector of genes to exclude from the analysis
#' @param genes a character vector of genes in case you want to limit the
#' analysis on a subset of genes
#' @param modifications the list of the modifications to study
#' @param advMethod the name of the method to use
#' @param advFixedValue the numeric value to use in case of
#' advMethod=`fixed`
#' @param advFct the function to use in case advMethod
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param firstDichot the initial number of slices before
#' the dichotomic search
#' @param maxSplitSize max size of dichotomic slices.
#' @param changeType `any` consider each misclassification,
#'  `not_na` consider each misclassification but NA.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a data.frame storing the number of min change attack 
#' possible for each cell type and each modification.
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'    c("T cell", 0.9)
#' }
#' rna_expression <- data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#'
#' minChangeOverview(rna_expression, clusters_id,
#' MyClassifier, modifications = list(c("perc1"), c("perc99")))
#' 
#' myModif = function(x){
#'    return(sample(1:10,1))
#' }
#' 
#' my_modifications = list(c("perc1"),
#'                         c("fixed", 1000),
#'                         c("full_matrix_fct", myModif))
#' minChangeOverview(rna_expression, clusters_id,
#'  MyClassifier, modifications = my_modifications)
#' 
#' @export
minChangeOverview <- function(exprs, clusters, classifier, exclGenes = c(),
            genes = c(), modifications = list(c("perc1"), c("perc99")),
            advMethod = "perc99", advFixedValue = 3, advFct = NULL,
            firstDichot = 100, maxSplitSize = 100, changeType = "any",
            verbose = FALSE) {
    if ( !is.matrix(exprs) && !is.data.frame(exprs) && !is(exprs,"DFrame")){
        stop("The argument exprs must be a matrix, a data.frame or a DataFrame.")
    }
    if (!is.character(clusters)) {
        stop("The argument clusters must be a vector of character.")
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
    if (!is.list(modifications)) {
        stop("The argument modifications must be a list.")
    }
    if (!is.character(advMethod)) {
        stop("The argument advMethod must be character.")
    }
    if (!is.numeric(firstDichot)){
        stop("The argument firstDichot must be numeric.")
    }
    if (!is.numeric(maxSplitSize)){
        stop("The argument maxSplitSize must be numeric.")
    }
    if (!is.character(changeType)) {
        stop("The argument changeType must be character.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }
    if (is(exprs, "DFrame")){
        exprs <- as.data.frame(exprs)
    }

    if (length(modifications) == 0) {
        df_result <- .minOverArgModifs(exprs, clusters, classifier, exclGenes,
                            genes, advMethod, advFixedValue, advFct,
                            firstDichot, maxSplitSize, changeType, verbose)
    } else {
        df_result <- .minOverListModifs(exprs, clusters, classifier, exclGenes,
                            genes, modifications, firstDichot, maxSplitSize,
                            changeType, verbose)
    }
    df_result
}

.gridWarning <- function(modifications, genes, iamsure){
    if ((length(modifications) + 1)^length(genes) > 100000 & !iamsure) {
        message("Exit because of too many combinations to test: ",
            (length(modifications) + 1)^length(genes))
        message("This will probably make your computer freeze.")
        message("You should lower the number of genes and/or",
            " of modifications. For example 5 genes and 2 modifications",
            " gives 243 combinations to test")
        message("You can use the iamsure=TRUE option",
            "to run the function anyway.")
        return(TRUE)
    } else {
        return (FALSE)
    }
}

#' Grid search of min change adversarial attack. Tries each
#' combination on a cluster, given a list of genes and a list of modifications.
#'
#' @details This function aims to find the shortest combination of genes
#' allowing to make a min change attack. It will test every possible combination
#' for a given gene list. This function can take a long time to run, and we
#' recommand to use the random walk search advRandWalkMinChange function instead
#' for lists above 10 genes.
#' 
#' You can specify a list of modifications as so - each item of the list should
#' be 1 or 2 length size.
#' The 1 length vector must contain the prerecorded modifications, 'perc1'
#' or 'perc99'.
#' The 2 length vector must have as first item:
#'  - 'fixed', in this case the second item should be the value to be
#'  replaced by.
#'  - 'full_row_fct', 'target_row_fct', 'target_matrix_fct' or
#'  'full_matrix_fct'. In this case the second item should be a function.
#' Let's say we want to analysis the susceptibility to min change attack
#' for 3 modifications: "perc1", the modification of each value of the
#' cluster by 1000, and a custom modification stored inside a function myFct.
#' Then the 'modification' parameter should be:
#' my_modifications = list(c("perc1"),
#'                         c("fixed", 1000),
#'                         c("full_matrix_fct", myFct))
#' @param exprs a matrix, a data.frame or a DataFrame of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param clusters a character vector of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format
#' @param genes the character vector of genes to study
#' @param modifications the list of the modifications to study
#' @param returnFirstFound set to TRUE to return result when a
#' the first misclassification is found
#' @param verbose logical, set to TRUE to activate verbose mode
#' @param iamsure logical, prevents from expansive calculations
#' when `genes` list is too long, set to `TRUE` to run anyway.
#' @return data.frame results of the classification of all the grid combinations
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'    c("T cell", 0.9)
#' }
#' rna_expression <- data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#' 
#' advGridMinChange(rna_expression, clusters_id, "T cell",
#' MyClassifier, genes=genes,
#' modifications = list(c("perc1"), c("perc99")))
#' 
#' myModif = function(x){
#'    return(sample(1:10,1))
#' }
#' 
#' my_modifications = list(c("perc1"),
#'                         c("fixed", 1000),
#'                         c("full_matrix_fct", myModif))
#' advGridMinChange(rna_expression, clusters_id, "T cell",
#' MyClassifier, genes=genes, modifications = my_modifications)
#' 
#' @export
advGridMinChange <- function(exprs, clusters, target, classifier,
                genes, modifications = list(c("perc1"), c("perc99")),
                returnFirstFound = FALSE, verbose = FALSE, iamsure = FALSE) {
    if ( !is.matrix(exprs) && !is.data.frame(exprs) && !is(exprs,"DFrame")){
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
    if (!is.character(genes)) {
        stop("The argument genes must be character or vector of character.")
    }
    if (!is.list(modifications)) {
        stop("The argument modifications must be a list.")
    }
    if (!is.logical(returnFirstFound)){
        stop("The argument returnFirstFound must be logical.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }
    if (!is.logical(iamsure)){
        stop("The argument iamsure must be logical.")
    }
    if (is(exprs, "DFrame")){
        exprs <- as.data.frame(exprs)
    }

    function_results <- data.frame(matrix(ncol = length(genes) + 5, nrow = 0))
    colnames(function_results) <- c("prediction", "odd", "genes_modified",
        "type_modified", "iteration", genes)
    if (.gridWarning(modifications, genes, iamsure))
        return(function_results)
    tests_grid <- gtools::permutations(n = length(modifications) + 1,
        r = length(genes), repeats.allowed = TRUE)
    results <- data.frame(matrix(ncol = length(genes) + 2, nrow = 0))
    colnames(results) <- c(genes, "prediction", "odd"); i <- 1
    while (i <= nrow(tests_grid)) {
        message("Running combination: ", i, " on ", nrow(tests_grid))
        exprs_temp <- exprs; row_results <- c()
        for (gene_ind in seq_len(length(genes))) {
            modif_ind <- tests_grid[i, gene_ind]
            if (modif_ind != length(modifications) + 1) {
                mod1 <- modifications[[modif_ind]][[1]]
                if (length(modifications[[modif_ind]]) == 1) {
                    exprs_temp <- advModifications(exprs_temp, genes[gene_ind],
                        clusters, target, advMethod = mod1, verbose = verbose)
                } else { mod2 <- modifications[[modif_ind]][[2]]
                    exprs_temp <- advModifications(exprs_temp, genes[gene_ind],
                        clusters, target, advMethod = mod1, advFct = mod2,
                        advFixedValue = mod2, verbose = verbose)}
                row_results <- c(row_results,
                    paste(modifications[[modif_ind]], collapse = " "))
            } else { row_results <- c(row_results, "NA")}
        }
        class_results <- classifier(exprs_temp, clusters, target)
        row_results <- c(row_results, class_results[1], class_results[2])
        results <- rbind(results, row_results)
        colnames(results) <- c(genes, "prediction", "odd")
        ifelse(class_results[1] == target, i <- i + 1, i <- nrow(tests_grid)+1)
    }
    results$genes_modified <-
        length(genes) - apply(results, 1, function(x) sum(x == "NA"))
    results$type_modified <- results$prediction != target
    results <- results[, c( "prediction", "odd",
            "genes_modified", "type_modified", genes)]
    function_results <- rbind(results, function_results)
    function_results <-
        function_results[order(as.numeric(function_results$genes_modified)), ]
    function_results <- function_results[order(function_results$type_modified,
            decreasing = TRUE), ]
    function_results
}

.randWalkTryNewVector <- function(exprs, genes, newWalkParams, modifications,
        clusters, target, classifier, bestAttackParams, i,
        function_results, changeType, verbose){
    exprs_temp <- exprs
    row_results <- c()
    rowResultsInt <- c()
    # Try the new params
    for (gene_ind in seq_len(length(genes))) {
        modif_ind <- as.numeric(newWalkParams[gene_ind])
        if (modif_ind != length(modifications) + 1) {
            mod1 <- modifications[[modif_ind]][[1]]
            if (length(modifications[[modif_ind]]) == 1) {
                exprs_temp <- advModifications(exprs_temp,
                    genes[gene_ind], clusters, target,
                    advMethod = mod1, verbose = verbose)
            } else {
                mod2 <- modifications[[modif_ind]][[2]]
                exprs_temp <- advModifications(exprs_temp, genes[gene_ind],
                    clusters, target, advMethod = mod1,
                    advFixedValue = mod2, advFct = mod2,
                    verbose = verbose)
            }
            row_results <- c(row_results, paste(modifications[[modif_ind]],
                        collapse = " "))
            rowResultsInt <- c(rowResultsInt, modif_ind)
        } else {
            row_results <- c(row_results, "NA")
            rowResultsInt <- c(rowResultsInt, modif_ind)
        }
    }
    genes_modified <- length(genes) - sum(row_results == "NA")
    previous_genes_modified <- length(genes) -
        sum(bestAttackParams == (length(modifications) + 1))
    class_results <- classifier(exprs_temp, clusters, target)
    type_modified <- class_results[1] != target
    if ( changeType == "not_na")
        type_modified <- (class_results[1] != target) &&
            class_results[1] != "NA"
    if (type_modified && genes_modified <= previous_genes_modified) {
        bestAttackParams <- rowResultsInt
        message( "Better attack with only ", genes_modified,
            " genes modified")
    }
    row_results <- c(class_results[1], class_results[2],
        genes_modified, type_modified, (i + 1), row_results)
    function_results <- rbind(row_results, function_results)
    function_results
}

.randWalkNewVector <- function(previousStrComb, modifications,
    bestAttackParams, stepChangeRatio, whileMaxCount, function_results){
    newWalkParams <- c()
    new_walk_step <- -1
    kill_count <- 0
    while (new_walk_step %in% previousStrComb ||
        length(newWalkParams[
            newWalkParams != (length(modifications) + 1)
        ]) >=
            length(bestAttackParams[bestAttackParams !=
                (length(modifications) + 1)])) {
        newWalkParams <- c()
        for (param in bestAttackParams) {
            if (param != length(modifications) + 1) {
                if (sample(seq_len(round(1 / stepChangeRatio)), 1) == 1) {
                    newWalkParams <- c( newWalkParams,
                        sample(seq_len(length(modifications) + 1), 1))
                } else {
                    newWalkParams <- c(newWalkParams, param)
                }
            } else {
                if (sample(seq_len(round(sum(newWalkParams ==
                    (length(modifications) + 1))^(1 / 2) /
                    stepChangeRatio)), 1) == 1) {
                    newWalkParams <- c( newWalkParams,
                        sample(seq_len(length(modifications) + 1), 1))
                } else {
                    newWalkParams <- c(newWalkParams, param)
                }
            }
        }
        new_walk_step <- paste(newWalkParams, collapse = " ")
        # If inifite loop kill function and return results
        if (kill_count > whileMaxCount) {
            function_results <- function_results[order(
                    as.numeric(function_results$genes_modified)), ]
            function_results <- function_results[order(
                function_results$type_modified, decreasing = TRUE), ]
            message("Inifite loop, kill function and return results")
            return(list(FALSE, function_results))
        }
        kill_count <- kill_count + 1
    }
    return(list(TRUE, newWalkParams, new_walk_step))
}

.randWalkBeforeWalk <- function(exprs, genes, modifications,
            clusters, target, classifier, firstBatch, verbose){
    init_objs <- .initRandWalk(modifications, genes, firstBatch)
    previousStrComb <- init_objs[[1]]; tests_grid <- init_objs[[2]]
    function_results <- init_objs[[3]]; functionResultsInt <- init_objs[[4]]
    results <- init_objs[[5]]; resultsInt <- init_objs[[6]]
    rw_seed <- .randWalkGetSeed(tests_grid, exprs, genes, modifications,
        clusters, target, classifier, results, resultsInt,
        previousStrComb, verbose)
    results <- rw_seed[[1]]; resultsInt <- rw_seed[[2]]
    previousStrComb <- rw_seed[[3]]

    results$genes_modified <-
        length(genes) - apply(results, 1, function(x) sum(x == "NA"))
    results$type_modified <- results$prediction != target
    results$iteration <- 1
    results <- results[, c("prediction", "odd", "genes_modified",
            "type_modified", "iteration", genes)]
    resultsInt$genes_modified <- length(genes) - apply( resultsInt, 1,
                        function(x) sum(x == (length(modifications) + 1)))
    resultsInt$type_modified <- resultsInt$prediction != target
    resultsInt$iteration <- 1
    resultsInt <- resultsInt[, c("prediction", "odd", "genes_modified",
        "type_modified", "iteration", genes)]
    function_results <- rbind(results, function_results)
    function_results <-
        function_results[order(function_results$genes_modified), ]
    function_results <- function_results[order(function_results$type_modified,
            decreasing = TRUE), ]
    functionResultsInt <- rbind(resultsInt, functionResultsInt)
    functionResultsInt <-
        functionResultsInt[order(functionResultsInt$genes_modified), ]
    functionResultsInt <- functionResultsInt[
        order(functionResultsInt$type_modified, decreasing = TRUE), ]
    bestAttackParams <- functionResultsInt[1, 6:ncol(functionResultsInt)]
    list(previousStrComb, function_results, bestAttackParams)
}

.randWalkGetSeed <- function(tests_grid, exprs, genes, modifications, clusters,
    target, classifier, results, resultsInt, previousStrComb, verbose){
    i<- 1; newPreviousStrComb <- c()
    while (i <= nrow(tests_grid)) {
        message("Running first batch to determine walk seed: ",
            i, " on ", nrow(tests_grid))
        exprs_temp <- exprs
        newPreviousStrComb <- c(newPreviousStrComb, previousStrComb[i])
        row_results <- c()
        rowResultsInt <- c()
        for (gene_ind in seq_len(length(genes))) {
            modif_ind <- tests_grid[i, gene_ind]
            if (modif_ind != length(modifications) + 1) {
                mod1 <- modifications[[modif_ind]][[1]]
                if (length(modifications[[modif_ind]]) == 1) {
                    exprs_temp <- advModifications(exprs_temp,
                        genes[gene_ind], clusters, target,
                        advMethod = mod1, verbose = verbose)
                } else {
                    mod2 <- modifications[[modif_ind]][[2]]
                    exprs_temp <- advModifications(exprs_temp,
                        genes[gene_ind], clusters, target,
                        advMethod = mod1,advFixedValue = mod2,
                        advFct = mod2, verbose = verbose)
                }
                row_results <- c( row_results,
                    paste(modifications[[modif_ind]], collapse = " "))
                rowResultsInt <- c(rowResultsInt, modif_ind)
            } else {
                row_results <- c(row_results, "NA")
                rowResultsInt <- c(rowResultsInt, modif_ind)
            }
        }
        class_results <- classifier(exprs_temp, clusters, target)
        row_results <- c(row_results, class_results[1], class_results[2])
        rowResultsInt <- c(unlist(rowResultsInt),
                    class_results[1], class_results[2])
        results <- rbind(results, row_results)
        resultsInt <- rbind(resultsInt, rowResultsInt)
        colnames(results) <- c(genes, "prediction", "odd")
        colnames(resultsInt) <- c(genes, "prediction", "odd")
        if (class_results[1] == target) {
            i <- i + 1
        } else {
            i <- nrow(tests_grid) + 1
        }
    }
    list(results, resultsInt, newPreviousStrComb)
}

.initRandWalk <- function(modifications, genes, firstBatch){
    previousStrComb <- c(-1)
    if ((length(modifications) + 1)^length(genes) > firstBatch) {
        tests_grid <- data.frame(matrix(ncol = length(genes), nrow = 0))
        while (nrow(tests_grid) < (firstBatch - length(modifications))) {
            rand_sample <- sample(seq_len(length(modifications) + 1),
                length(genes),
                replace = TRUE
            )
            str_rand_sample <- paste(rand_sample, collapse = " ")
            if (!str_rand_sample %in% previousStrComb) {
                tests_grid <- rbind(tests_grid, rand_sample)
            }
            previousStrComb <- c(previousStrComb, str_rand_sample)
        }
        for (i in seq_len(length(modifications))) {
            tests_grid <- rbind(rep(i, length(genes)), tests_grid)
        }
        tests_grid <- as.matrix(tests_grid)
    } else {
        tests_grid <- gtools::permutations(n = length(modifications) + 1,
            r = length(genes), repeats.allowed = TRUE)
    }
    function_results <- data.frame(matrix(ncol = length(genes) + 5, nrow = 0))
    colnames(function_results) <- c("prediction", "odd", "genes_modified",
        "type_modified", "iteration", genes)
    functionResultsInt <- data.frame(
        matrix(ncol = length(genes) + 5, nrow = 0)
    )
    colnames(functionResultsInt) <- c( "prediction", "odd",
        "genes_modified", "type_modified", "iteration", genes)
    results <- data.frame(matrix(ncol = length(genes) + 2, nrow = 0))
    colnames(results) <- c(genes, "prediction", "odd")
    resultsInt <- data.frame(matrix(ncol = length(genes) + 2, nrow = 0))
    colnames(resultsInt) <- c(genes, "prediction", "odd")
    return (list(previousStrComb, tests_grid, function_results,
        functionResultsInt, results, resultsInt))
}

#' Random walk search of min change adversarial attack.
#' 
#' @details The parameter search by grid can take a long time,
#' this function aims to make a parameter search faster.
#' We have a function, advMinChange, looking for one gene attacks.
#' The advRandWalkMinChange function aims to find a multiple genes
#' attack, with the fewer genes possible.
#' At first the user have to provide a list of genes to test, for
#' example by running differential statistics between two cell clusters.
#' The user should also provide a list of modifications to test, to
#' define as so - each item of the list should be 1 or 2 length size.
#' The 1 length vector must contain the prerecorded modifications, 'perc1'
#' or 'perc99'.
#' The 2 length vector must have as first item:
#'  - 'fixed', in this case the second item should be the value to be
#'  replaced by.
#'  - 'full_row_fct', 'target_row_fct', 'target_matrix_fct' or
#'  'full_matrix_fct'. In this case the second item should be a function.
#' Let's say we want to analysis the susceptibility to min change attack
#' for 3 modifications: "perc1", the modification of each value of the
#' cluster by 1000, and a custom modification stored inside a function myFct.
#' Then the 'modification' parameter should be:
#' my_modifications = list(c("perc1"),
#'                         c("fixed", 1000),
#'                         c("full_matrix_fct", myFct))
#' 
#' Then the function will try to find the best combination of these genes
#' and modifications to make the min change attack.
#' Step 1 is to find a seed by trying random combinations of
#' genes and modifications on a cluster until the classification is altered.
#' Step 2 is to perform a random walk search to reduce the number of genes
#' needed to change the classification.
#' The 
#' @param exprs a matrix, a data.frame or a DataFrame of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param clusters a character vector of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format
#' @param genes the character vector of genes to study
#' @param modifications the list of the modifications to study
#' @param firstBatch the maximum number of try in step 1
#' @param walkLength the maximum number of try in step 2
#' @param stepChangeRatio ratio of parameters change in new walk step
#' @param whileMaxCount the maximum number of try when looking
#' for new combination of parameters
#' @param changeType `any` consider each misclassification,
#'  `not_na` consider each misclassification but NA.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return data.frame results of the classification of all the grid combinations
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'    c("T cell", 0.9)
#' }
#' rna_expression <- data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#'
#' advRandWalkMinChange(rna_expression, clusters_id, "T cell",
#' MyClassifier, genes=genes,
#' modifications = list(c("perc1"), c("perc99")))
#'
#' # Stop at first attack discovery, whitout going into the walk
#' # parameter search.
#' advRandWalkMinChange(rna_expression, clusters_id, "T cell",
#'  MyClassifier, genes=genes,
#'  modifications = list(c("perc1"), c("perc99")), walkLength=0)
#' 
#' myModif = function(x){
#'    return(sample(1:10,1))
#' }
#' 
#' my_modifications = list(c("perc1"),
#'                         c("fixed", 1000),
#'                         c("full_matrix_fct", myModif))
#' advRandWalkMinChange(rna_expression, clusters_id, "T cell",
#'  MyClassifier, genes=genes, modifications = my_modifications)
#' @export
advRandWalkMinChange <- function(exprs, clusters, target, classifier, genes,
        modifications = list(c("perc1"), c("perc99")), firstBatch = 100,
        walkLength = 100, stepChangeRatio = 0.2, whileMaxCount = 10000,
        changeType = "any", verbose = FALSE) {
    if ( !is.matrix(exprs) && !is.data.frame(exprs) && !is(exprs,"DFrame")){
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
    if (!is.character(genes)) {
        stop("The argument genes must be character or vector of character.")
    }
    if (!is.list(modifications)) {
        stop("The argument modifications must be a list.")
    }
    if (!is.numeric(firstBatch)) {
        stop("The argument firstBatch must be numeric.")
    }
    if (!is.numeric(walkLength)) {
        stop("The argument walkLength must be numeric.")
    }
    if (!is.numeric(stepChangeRatio)) {
        stop("The argument stepChangeRatio must be numeric.")
    }
    if (!is.numeric(whileMaxCount)) {
        stop("The argument whileMaxCount must be numeric.")
    }
    if (!is.character(changeType)) {
        stop("The argument changeType must be character.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }
    if (is(exprs, "DFrame")){
        exprs <- as.data.frame(exprs)
    }

    before_walk <- .randWalkBeforeWalk(exprs, genes, modifications,
        clusters, target, classifier, firstBatch, verbose)
    previousStrComb <- before_walk[[1]]
    function_results <- before_walk[[2]]
    bestAttackParams <- before_walk[[3]]
    if (function_results[1, "type_modified"] == FALSE) {
        message("No modified type, try with a higher firstBatch argument")
        return(function_results)
    }
    for (i in seq_len(walkLength)) {
        message("Walk step ", i, " on ", walkLength)
        if (verbose) {
            message("Current best attack length: ",
                length(bestAttackParams), " genes")
        }
        # Create a new vector of parameters to test
        rw_new_vector <- .randWalkNewVector(previousStrComb, modifications,
            bestAttackParams, stepChangeRatio, whileMaxCount,
            function_results)
        # If infinite loop kill function and return result
        if ( !rw_new_vector[[1]] ) {return(rw_new_vector[[2]])}
        newWalkParams <- rw_new_vector[[2]]
        new_walk_step <- rw_new_vector[[3]]

        previousStrComb <- c(previousStrComb, new_walk_step)
        function_results <- .randWalkTryNewVector(exprs, genes,
            newWalkParams, modifications, clusters, target,
            classifier, bestAttackParams, i, function_results,
            changeType, verbose)
    }
    function_results <- function_results[order(
            as.numeric(function_results$genes_modified)), ]
    function_results <- function_results[order(
            function_results$type_modified, decreasing = TRUE), ]
    function_results
}
