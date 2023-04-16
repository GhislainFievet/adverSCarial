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

