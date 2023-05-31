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
#' @param exprs can be a matrix or a data.frame of numeric RNA expression,
#' cells are rows and genes are columns. Or can be a SingleCellExperiment object.
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
#' @param argForClassif the type of the first argument to feed to the
#' classifier function. 'data.frame' by default, can be 'SingleCellExperiment'
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a list of genes/new classification tuples
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'    c("B cell", 0.9)
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
        returnFirstFound = FALSE, changeType = "any",
        argForClassif = 'data.frame',
        verbose = FALSE) {
    if (!is(exprs, 'matrix') && !is(exprs,'data.frame') && !is(exprs,'SingleCellExperiment')){
        stop("The argument exprs must be a matrix, a data.frame or a SingleCellExperiment")
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
    if (!is.character(argForClassif)) {
        stop("The argument argForClassif must be character: 'data.frame' or 'SingleCellExperiment'.")
    }
    if (!is.character(changeType)) {
        stop("The argument changeType must be a vector of character.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }

    if (is(exprs,'SingleCellExperiment') ){
        exprs <- as.matrix(t(counts(exprs)))
    }

    genesIndex <- seq_len(ncol(exprs))
    if (length(exclGenes) != 0 && length(genes) == 0) {
        genesToRemove <- unlist(lapply(exclGenes, function(strGene){
            if (!is.na(match(strGene, colnames(exprs)))) {
                return(-match(strGene, colnames(exprs)))
            } else {
                return(NULL)
            }
        }))
        if (length(genesToRemove) != 0)
            genesIndex <- genesIndex[genesToRemove]
    }
    if (length(genes) != 0) {
        exclGenes <- unique(exclGenes)
        genes <- unique(genes)
        genesToKeep <- genes[!genes %in% exclGenes]
        indGenesToKeep <- unlist(lapply(genesToKeep), function(strGene){
            if (!is.na(match(strGene, colnames(exprs)))) {
                return(match(strGene, colnames(exprs)))
            } else {
                return(NULL)
            }
        })
        if (length(indGenesToKeep) != 0)
            genesIndex <- genesIndex[indGenesToKeep]
    }
    cSplitIndex <- .custSplit(genesIndex, seq_len(firstDichot))
    lastResLength <<- 0
    lSplitsResults <- unlist(lapply(seq_len(firstDichot), function(i){
        prevTime <- Sys.time()
        if (lastResLength > 0 && returnFirstFound){
            return(NULL)
        }
        message(paste0("Split number: ", i, "/", firstDichot))
        if (length(unlist(cSplitIndex[i])) != 0) {
            temp_result <- .predictionDichotMinChange(exprs,
                    colnames(exprs)[unlist(cSplitIndex[i])],
                    clusters, target, classifier, advMethod = advMethod,
                    advFixedValue = advFixedValue, advFct = advFct,
                    maxSplitSize = maxSplitSize,
                    returnFirstFound = returnFirstFound,
                    changeType=changeType,
                    argForClassif=argForClassif,
                    verbose = verbose)
            if ( lastResLength == 0 ){
                lastResLength <<- length(temp_result)
            }
            return(temp_result)
        }
        message(paste0("Split time: ", Sys.time() - prevTime))
    }), recursive = FALSE)
    if ( is.null(lSplitsResults)){
        NULL
    } else {
        new("advList", values=lSplitsResults)
    }
}

.dichotMinSplit <- function(lResults, cGeneSplitValue, exprs,
            genes, clusters, target, classifier, advMethod, advFixedValue,
            advFct, maxSplitSize, changeType, returnFirstFound, argForClassif, 
            verbose){
    if (length(cGeneSplitValue) != 0) {
        if (verbose) {
            message("before predictWithNewValue")
            message(length(cGeneSplitValue))}
        modObj <- predictWithNewValue(exprs, cGeneSplitValue,
            clusters, target, classifier,
            advMethod = advMethod,
            advFixedValue = advFixedValue, advFct = advFct,
            argForClassif = argForClassif, verbose = TRUE)
        cellType <- modObj[1]
        celltypeScore <- modObj[2]
        if (verbose){
            message(cellType)
            message(celltypeScore)
        }
        if (cellType != target) {
            if (length(cGeneSplitValue) <= maxSplitSize &&
                (changeType == "any" || cellType != "NA" )) {
                message("Store gene:")
                message(cGeneSplitValue)
                res_key <- paste(cGeneSplitValue, collapse = "__")
                lResults[[res_key]] <- c(cellType, celltypeScore)
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
                        changeType = changeType,
                        argForClassif = argForClassif,
                        verbose = verbose)
                }
            }
        }
    }
    lResults
    
}

.predictionDichotMinChange <- function(exprs, genes, clusters, target,
        classifier, lResults = list(), advMethod = "perc99",
        advFixedValue = 3, advFct = NULL, maxSplitSize = 1,
        returnFirstFound = FALSE, changeType = "any", argForClassif,
        verbose = TRUE) {
    if ( returnFirstFound && length(lResults) > 0){
        return (lResults)
    }
    cGeneSplit <- .custSplit(unlist(genes), c(1,2))
    message(paste0("genes size: ", length(unlist(genes))))
    lResults <- .dichotMinSplit(lResults,
                    unname(unlist(cGeneSplit[1])), exprs,
                    genes, clusters, target, classifier, advMethod,
                    advFixedValue, advFct,
                    maxSplitSize, changeType,returnFirstFound, argForClassif,
                    verbose)
    if ( returnFirstFound && length(lResults) > 0){
        return (lResults)
    } 
    lResults <- .dichotMinSplit(lResults,
                    unname(unlist(cGeneSplit[2])), exprs,
                    genes, clusters, target, classifier,
                    advMethod, advFixedValue,
                    advFct, maxSplitSize, changeType, returnFirstFound, argForClassif,
                    verbose)
    lResults
}




