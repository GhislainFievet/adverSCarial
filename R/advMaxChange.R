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
#' makes its classification. Then three possible scenarios:
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
#' @param exprs DelayedMatrix of numeric RNA expression, cells are rows and genes
#' are columns - or a SingleCellExperiment object, a matrix or a data.frame. By default,
#' these are converted to a data.frame to increase speed performance during modifications.
#' However, this conversion can consume a significant amount of memory, see 'argForModif'
#' argument for options.
#' @param clusters a character vector of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format.
#' A classifier function should be formated as follow:
#' classifier = function(expr, clusters, target){
#'      # Making the classification
#'      c("cell type", score)
#' }
#' `score` should be numeric between 0 and 1, 1 being the highest confidance
#' into the cell type classification.
#' The matrix `expr` contains RNA expression values, the vector `clusters`
#' consists of the cluster IDs for each cell in `expr`, and `target` is the
#' ID of the cluster for which we want to have a classification.
#' The function returns a vector with the classification result, and a score.
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
#' @param argForClassif the type of the first argument to feed to the
#' classifier function. 'DelayedMatrix' by default, can be 'SingleCellExperiment'
#' or 'data.frame'.
#' @param argForModif type of matrix during for the modification, 'DelayedMatrix'
#' by default. Can be 'data.frame', which is faster, but need more memory.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a character vector of genes you can modify on a cluster without
#' modifying its classification
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'    c("T cell", 0.9)
#' }
#' rna_expression <- DelayedArray(data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3)))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#' 
#' advMaxChange(rna_expression, clusters_id,
#' "T cell", MyClassifier, advMethod="perc99")
#' @export
advMaxChange <- function(exprs, clusters, target, classifier,
                        exclGenes = c(), genes = c(), advMethod = "perc99",
                        advFixedValue = 3, advFct = NULL,
                        maxSplitSize = 1, argForClassif = 'DelayedMatrix',
                        argForModif = 'data.frame',
                        verbose = FALSE) {
    if (!is(exprs, 'matrix') && !is(exprs,'data.frame') &&
        !is(exprs,'SingleCellExperiment') && !is(exprs,'DelayedMatrix')){
        stop("The argument exprs must be a DelayedMatrix, a SingleCellExperiment, a matrix or a data.frame")
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
    if (!is.character(argForClassif)) {
        stop("The argument argForClassif must be character: 'data.frame' or 'SingleCellExperiment'.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }
    if (is(exprs,'SingleCellExperiment') ){
        exprs <- t(counts(exprs))
    }
    if (!is(exprs,'DelayedMatrix') && argForModif=="DelayedMatrix"){
        message("Converting exprs object to a DelayedArray object")
        exprs <- DelayedArray::DelayedArray(exprs)
    }
    if (!is(exprs,'data.frame') && argForModif=="data.frame"){
        message("Converting exprs object to a data.frame object")
        exprs <- as.data.frame(exprs)
    }

    if (maxSplitSize < 1){ maxSplitSize <- 1 }
    genesToKeep <- unlist(lapply(unique(genes), function(strGene){
        if (!is.na(match(strGene, colnames(exprs)))) {
            return(strGene)
        } else {
            return(NULL)
        }
    }))
    if (length(genes) == 0) {
        genesToKeep <- colnames(exprs)
    }
    for (strGene in unique(exclGenes)) {
        if (!is.na(match(strGene, genesToKeep))) {
            genesToKeep <- genesToKeep[-match(strGene, genesToKeep)]
        }
    }

    lResults <- .predictionDichotMaxChange(exprs, genesToKeep,
        clusters, target, classifier,
        advMethod = advMethod,
        advFixedValue = advFixedValue, advFct = advFct,
        maxSplitSize = maxSplitSize, argForClassif = argForClassif,
        argForModif=argForModif,
        verbose = verbose
    )
    if ( is.null(lResults)){
        NULL
    } else {
        new("advChar", values=lResults)
    }
}

.dichotMaxSameType <- function(lResults, cGeneSplitValue, exprs,
                clusters, target, classifier, advMethod,
                advFixedValue, advFct, maxSplitSize, argForClassif, argForModif, verbose){
        if (verbose) { message("same cellType") }
        if (length(lResults) == 0) {
            lResults <- cGeneSplitValue
        } else {
            if (verbose){
                message("check if concat lists still gives target")
            }
            concatGenesList <- c(lResults, cGeneSplitValue)
            if (verbose) {
                message(paste0(
                    "length of tested gene list: ",
                    length(unique(concatGenesList))
                ))
            }
            # check if concat lists still gives target
            modObj <- predictWithNewValue(exprs, concatGenesList, clusters,
                target, classifier,
                advMethod = advMethod,
                advFixedValue = advFixedValue, advFct = advFct,
                argForClassif = argForClassif, argForModif=argForModif, verbose = verbose
            )
            concatCellType <- modObj[1]
            if (concatCellType == target) {
                if (verbose) { message("YES: merge results") }
                # Merge results
                lResults <- c(lResults, cGeneSplitValue)
                if (verbose) {
                    message(paste0( "result length after merge: ", length(lResults)))
                } else {
                    message(paste0( "result length: ", length(lResults)))
                }
            } else {
                if (verbose) { message("NO: split and retry") }
                if (length(cGeneSplitValue) > maxSplitSize) {
                    lResults <- .predictionDichotMaxChange(exprs,
                        cGeneSplitValue, clusters, target,
                        classifier,
                        lResults = lResults,
                        advMethod = advMethod,
                        advFixedValue = advFixedValue,
                        advFct = advFct,
                        maxSplitSize = maxSplitSize,
                        argForClassif = argForClassif,
                        argForModif=argForModif,
                        verbose = verbose
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
                                    argForClassif,
                                    argForModif,
                                    verbose){
    if (length(cGeneSplitValue) != 0) {
        if (verbose) {
            message("before predictWithNewValue 1")
            message(length(cGeneSplitValue))
        }
        modObj <- predictWithNewValue(exprs, cGeneSplitValue,
            clusters, target, classifier,
            advMethod = advMethod,
            advFixedValue = advFixedValue, advFct = advFct,
            argForClassif = argForClassif,
            argForModif=argForModif,
            verbose = verbose
        )
        cellType <- modObj[1]
        celltypeScore <- modObj[2]
        if (verbose){
            message(cellType)
            message(celltypeScore)
        }
    }
    if (cellType == target) {
        lResults <- .dichotMaxSameType(lResults,
                    cGeneSplitValue,
                    exprs, clusters,
                    target, classifier,
                    advMethod,
                    advFixedValue, advFct,
                    maxSplitSize,
                    argForClassif,
                    argForModif,
                    verbose)
    } else {
        if (verbose) { message("NOT same cellType") }
        if (length(cGeneSplitValue) > maxSplitSize) {
            lResults <- .predictionDichotMaxChange(exprs,
                unname(cGeneSplitValue), clusters, target,
                classifier, lResults = lResults, advMethod = advMethod,
                advFixedValue = advFixedValue, advFct = advFct,
                maxSplitSize = maxSplitSize, argForClassif = argForClassif,
                argForModif=argForModif, verbose = verbose
            )
        }
    }
    lResults
}

.custSplit <- function(v1, v2){
    result <- list()
    for ( i in seq_along(v1)){
        strKey = as.character(v2[(i %% length(v2)) + 1])
        if ( !strKey %in% names(result)){
            result[[strKey]] = v1[i]
        } else {
            result[[strKey]] = c(result[[strKey]], v1[i])
        }
    }
    result
}

.predictionDichotMaxChange <- function(exprs, genes, clusters, target,
                                    classifier, lResults = c(),
                                    advMethod = "perc99",
                                    advFixedValue = 3,
                                    advFct = NULL,
                                    maxSplitSize = 1,
                                    argForClassif,
                                    argForModif,
                                    verbose) {
    cGeneSplit <- .custSplit(unlist(genes), c(1,2))
    if (verbose){
        message(paste0("genes size: ", length(unlist(genes))))
        message(paste0("current gene results length: ", length(lResults)))
    }
    # Random order for each reccursion
    randomIndex <- sample(c(1,2), 2)
    cGeneSplit <- cGeneSplit[randomIndex]
    lResults <- .dichotMaxSplit(lResults, unname(unlist(cGeneSplit[1])),
                    exprs, genes, clusters, target, classifier, advMethod,
                    advFixedValue, advFct, maxSplitSize, argForClassif, argForModif, verbose)
    lResults <- .dichotMaxSplit(lResults, unname(unlist(cGeneSplit[2])),
                exprs, genes, clusters, target, classifier, advMethod,
                advFixedValue, advFct, maxSplitSize, argForClassif, argForModif, verbose)
    lResults
}

