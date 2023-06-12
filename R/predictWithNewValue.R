#' Returns a classification and an odd value from a
#' RNA expression DelayedMatrix or a SingleCellExperiment object,
#' for given genes, for a given cluster, for a given modification.
#'
#' @details This function aims to concatenate the following actions:
#' - modify the RNA gene expression
#' - classify the result
#' This is a widely used function in the other functions of the package.
#' @param exprs DelayedMatrix of numeric RNA expression, cells are rows and genes
#' are columns - or a SingleCellExperiment object, a matrix or a data.frame.
#' @param genes the character vector of genes to modify
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
#' @param advMethod the name of the method to use
#' @param advFixedValue the numeric value to use in case of
#' advMethod=`fixed`
#' @param advFct the function to use in case advMethod
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param argForClassif the type of the first argument to feed to the
#' classifier function. 'DelayedMatrix' by default, can be 'SingleCellExperiment'
#' or 'data.frame'.
#' @param argForModif type of matrix during for the modification, 'DelayedMatrix'
#' by default. Can be 'data.frame', which is faster, but need more memory.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a vector of the classification, and the associated odd
#' @examples
#' library(DelayedArray)
#' 
#' MyClassifier <- function(expr, clusters, target) {
#'    c("T cell", 0.9)
#' }
#' rna_expression <- DelayedArray(data.frame(CD4=c(0,0,0,0),
#'      CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3)))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#'
#' predictWithNewValue(rna_expression, genes, clusters_id,
#' "T cell", MyClassifier, advMethod="perc99")
#' @export
predictWithNewValue <- function(exprs, genes, clusters, target,
                                classifier, advMethod = "perc99",
                                advFixedValue = 3,
                                advFct = NULL,
                                argForClassif = 'DelayedMatrix',
                                argForModif = 'DelayedMatrix',
                                verbose = FALSE) {
    if (!is(exprs, 'matrix') && !is(exprs,'data.frame') &&
        !is(exprs,'SingleCellExperiment') && !is(exprs,'DelayedMatrix')){
        stop("The argument exprs must be a DelayedMatrix, a SingleCellExperiment, a matrix or a data.frame")
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
    if (!is.character(argForClassif)) {
        stop("The argument argForClassif must be character: 'data.frame' or 'SingleCellExperiment'.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }

    if (verbose) {
        message("Modify data for ", length(genes), " genes for cluster ", target)
    }

    modifExprs <- advModifications(exprs,
        genes,
        clusters,
        target,
        advMethod = advMethod,
        advFixedValue = advFixedValue,
        advFct = advFct,
        argForClassif = argForClassif,
        argForModif = argForModif
    )

    classifier(modifExprs, clusters, target)
}

