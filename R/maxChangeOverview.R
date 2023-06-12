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
#' @param exprs DelayedMatrix of numeric RNA expression, cells are rows and genes
#' are columns - or a SingleCellExperiment object, a matrix or a data.frame. By default,
#' these are converted to a data.frame to increase speed performance during modifications.
#' However, this conversion can consume a significant amount of memory, see 'argForModif'
#' argument for options.
#' @param clusters a character vector of the clusters to which the cells belong
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
#' @param argForClassif the type of the first argument to feed to the
#' classifier function. 'DelayedMatrix' by default, can be 'SingleCellExperiment'
#' or 'data.frame'.
#' @param argForModif type of matrix during for the modification, 'DelayedMatrix'
#' by default. Can be 'data.frame', which is faster, but need more memory.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a DataFrame storing the number of possible max change attacks
#' for each cell type and each modification.
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'    c("T cell", 0.9)
#' }
#' rna_expression <- DelayedArray(data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3)))
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
                            argForClassif = 'DelayedMatrix', argForModif = 'data.frame',
                            verbose = FALSE) {
    if (!is(exprs, 'matrix') && !is(exprs,'data.frame') &&
        !is(exprs,'SingleCellExperiment') && !is(exprs,'DelayedMatrix')){
        stop("The argument exprs must be a DelayedMatrix, a SingleCellExperiment, a matrix or a data.frame")
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

    if (length(modifications) == 0) {
        dfResult <- .maxOverArgModifs(exprs, clusters, classifier,
            exclGenes, genes, advMethod, advFixedValue,
            advFct, maxSplitSize, argForClassif, argForModif, verbose)
    } else {
        dfResult <- .maxOverListModifs(exprs, clusters, classifier, exclGenes,
            genes, modifications, maxSplitSize, argForClassif, argForModif, verbose)
    }
    S4Vectors::DataFrame(dfResult)
}


.maxOverListModifs <- function(exprs, clusters, classifier, exclGenes,
                            genes, modifications, maxSplitSize, argForClassif, argForModif, verbose){
    dfResult <- data.frame(todel = unique(clusters))
    rownames(dfResult) <- unique(clusters)
    dfNames <- unlist(lapply(seq_along(modifications), function(modifInd){
        paste(modifications[[modifInd]], collapse = "_")
    }))
    vecsAttacksLength <- lapply(seq_along(modifications), function(modifInd){
        mod1 <- modifications[[modifInd]][[1]]
        attacksLength <- unname(vapply(unique(clusters), function(cellType){
            if (verbose) {
                message(paste0("Running maxChange attack on ", cellType, ", with a maxSplitSize of: ", maxSplitSize))
                message("The smaller the maxSplitSize, the more precise the result will be, but it will take longer.")
                message("Modification: ", paste(modifications[[modifInd]], collapse = " "))
            }
            if (length(modifications[[modifInd]]) == 1) {
                maxChangeGenes <- advMaxChange(exprs, clusters, cellType,
                    classifier, advMethod = mod1,
                    maxSplitSize = maxSplitSize, exclGenes = exclGenes,
                    genes = genes, argForClassif = argForClassif, argForModif=argForModif, verbose = verbose)
            } else {
                mod2 <- modifications[[modifInd]][[2]]
                maxChangeGenes <- advMaxChange(exprs, clusters,
                    cellType, classifier, advMethod = mod1,
                    advFixedValue = mod2, advFct = mod2,
                    maxSplitSize = maxSplitSize, exclGenes = exclGenes,
                    genes = genes, argForClassif = argForClassif, argForModif=argForModif, verbose = verbose)
                if (!is.null(maxChangeGenes)){
                    maxChangeGenes <- maxChangeGenes@values
                }
            }
            resultLength <- length(maxChangeGenes)
            if (verbose) {
                message(paste0("At least ", resultLength, " genes can be modified with the ",
                    paste(modifications[[modifInd]], collapse = " "), " method, and the cluster will still",
                    " be classified as ", cellType))
            }
            return(resultLength)
        }, numeric(1)))
        return(attacksLength)
    })

    dfResult <- do.call(cbind, list(dfResult, vecsAttacksLength))

    dfResult$todel <- NULL
    colnames(dfResult) <- dfNames
    dfResult
}

.maxOverArgModifs <- function(exprs, clusters, classifier, exclGenes,
                            genes, advMethod, advFixedValue,
                            advFct, maxSplitSize, argForClassif, argForModif, verbose){
    attacksLength <- vapply(unique(clusters), function(cellType){
        if (verbose) {
            message("Running maxChange attack on ", cellType, ", with a maxSplitSize of: ", maxSplitSize)
            message("The smaller the maxSplitSize, the more precise the result will be, but it will take longer.")
        }
        maxChangeGenes <- advMaxChange(exprs, clusters, cellType,
            classifier, advMethod = advMethod,
            advFixedValue = advFixedValue, advFct = advFct,
            maxSplitSize = maxSplitSize, exclGenes = exclGenes,
            genes = genes, argForClassif = argForClassif, argForModif=argForModif, verbose = verbose
        )
        if (!is.null(maxChangeGenes)){
            maxChangeGenes <- maxChangeGenes@values
        }
        resultLength <- length(maxChangeGenes)
        if (verbose) {
            message( "At least ", resultLength, " genes can be modified with the ", advMethod,
                " method, and the cluster will still be classified as ", cellType)
        }
        return(resultLength)
    }, numeric(1))
    names(attacksLength) <- unique(clusters)
    dfResult <- data.frame(atLeastGeneNumber = unlist(attacksLength))
    rownames(dfResult) <- names(attacksLength)
    dfResult
}

