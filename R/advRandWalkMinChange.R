
#' Random walk search of min change adversarial attack.
#' 
#' @details The parameter search by grid can take a long time,
#' this function aims to make a parameter search faster.
#' We have a function, advSingleGene, looking for one gene attacks.
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
#' @param exprs DelayedMatrix of numeric RNA expression, cells are rows and genes
#' are columns - or a SingleCellExperiment object, a matrix or a data.frame. By
#' default matrix and data.frame are converted to DelayedMatrix for memory
#' performance, see 'argForModif' argument for options.
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
#' @param genes the character vector of genes to study
#' @param modifications the list of the modifications to study
#' @param firstBatch the maximum number of try in step 1
#' @param walkLength the maximum number of try in step 2
#' @param stepChangeRatio ratio of parameters change in new walk step
#' @param whileMaxCount the maximum number of try when looking
#' for new combination of parameters
#' @param changeType `any` consider each misclassification,
#'  `not_na` consider each misclassification but NA.
#' @param argForClassif the type of the first argument to feed to the
#' classifier function. 'DelayedMatrix' by default, can be 'SingleCellExperiment'
#' or 'data.frame'.
#' @param argForModif type of matrix during for the modification, 'DelayedMatrix'
#' by default. Can be 'data.frame', which is faster, but need more memory.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return DataFrame results of the classification of all the grid combinations
#' @examples
#' library(DelayedArray)
#' 
#' MyClassifier <- function(expr, clusters, target) {
#'    c("T cell", 0.9)
#' }
#' rna_expression <- DelayedArray(data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3)))
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
#' myModif = function(x, y){
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
        changeType = "any", argForClassif = 'DelayedMatrix',
        argForModif = 'DelayedMatrix',
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

    beforeWalk <- .randWalkBeforeWalk(exprs, genes, modifications,
        clusters, target, classifier, firstBatch, argForClassif, argForModif, verbose)
    previousStrComb <- beforeWalk[[1]]
    functionResults <- beforeWalk[[2]]
    bestAttackParams <- beforeWalk[[3]]
    if (functionResults[1, "typeModified"] == FALSE) {
        message("No modified type, try with a higher firstBatch argument")
        return(functionResults)
    }
    for (i in seq_len(walkLength)) {
        if (verbose || log2(i) == round(log2(i)) || i == walkLength){
            message("Walk step ", i, " on ", walkLength)
            message("Current best attack length: ",
                length(bestAttackParams), " genes")
        }
        # Create a new vector of parameters to test
        rwNewVector <- .randWalkNewVector(previousStrComb, modifications,
            bestAttackParams, stepChangeRatio, whileMaxCount,
            functionResults)
        # If infinite loop kill function and return result
        if ( !rwNewVector[[1]] ) {return(rwNewVector[[2]])}
        newWalkParams <- rwNewVector[[2]]
        newWalkStep <- rwNewVector[[3]]

        previousStrComb <- c(previousStrComb, newWalkStep)
        lrwtnv <- .randWalkTryNewVector(exprs, genes,
            newWalkParams, modifications, clusters, target,
            classifier, bestAttackParams, i, functionResults,
            changeType, argForClassif, argForModif, verbose)
        functionResults <- lrwtnv[['functionResults']]
        bestAttackParams <- lrwtnv[['bestAttackParams']]
    }
    functionResults <- functionResults[order(as.numeric(functionResults$genesModified)), ]
    functionResults <- functionResults[order(functionResults$typeModified, decreasing = TRUE), ]
    S4Vectors::DataFrame(functionResults)
}

.randWalkTryNewVector <- function(exprs, genes, newWalkParams, modifications,
        clusters, target, classifier, bestAttackParams, i,
        functionResults, changeType, argForClassif, argForModif, verbose){
    exprsTemp <- exprs
    rowResults <- c()
    rowResultsInt <- c()
    # Try the new params
    for (geneInd in seq_along(genes)) {
        modifInd <- as.numeric(newWalkParams[geneInd])
        if (modifInd != length(modifications) + 1) {
            mod1 <- modifications[[modifInd]][[1]]
            if (length(modifications[[modifInd]]) == 1) {
                exprsTemp <- advModifications(exprsTemp,
                    genes[geneInd], clusters, target,
                    advMethod = mod1, argForModif=argForModif, verbose = verbose)
            } else {
                mod2 <- modifications[[modifInd]][[2]]
                exprsTemp <- advModifications(exprsTemp, genes[geneInd],
                    clusters, target, advMethod = mod1,
                    advFixedValue = mod2, advFct = mod2, argForModif=argForModif,
                    verbose = verbose)
            }
            rowResults <- c(rowResults, paste(modifications[[modifInd]],
                        collapse = " "))
            rowResultsInt <- c(rowResultsInt, modifInd)
        } else {
            rowResults <- c(rowResults, "NA")
            rowResultsInt <- c(rowResultsInt, modifInd)
        }
    }
    genesModified <- length(genes) - sum(rowResults == "NA")
    previousGenesModified <- length(genes) - sum(bestAttackParams == (length(modifications) + 1))
    if ( argForClassif == 'SingleCellExperiment'){
        exprsTemp <- SingleCellExperiment(assays = list(counts = t(exprsTemp)))
    }
    if (!is(exprsTemp, 'data.frame') && argForClassif=='data.frame'){
        exprsTemp <- as.data.frame(exprsTemp)
    }
    if (!is(exprsTemp,'DelayedMatrix') && argForClassif=="DelayedMatrix"){
        exprsTemp <- DelayedArray::DelayedArray(exprsTemp)
    }
    classResults <- classifier(exprsTemp, clusters, target)
    typeModified <- classResults[1] != target
    if ( changeType == "not_na")
        typeModified <- (classResults[1] != target) &&
            classResults[1] != "NA"
    if (typeModified && genesModified <= previousGenesModified) {
        tempBestAttackParams <- as.data.frame(t(rowResultsInt))
        colnames(tempBestAttackParams) <- colnames(bestAttackParams)
        bestAttackParams <- tempBestAttackParams
        message( "Better attack with only ", genesModified, " genes modified")
    }
    rowResults <- c(classResults[1], classResults[2],
        genesModified, typeModified, (i + 1), rowResults)
    functionResults <- rbind(rowResults, functionResults)
    list(functionResults=functionResults, bestAttackParams=bestAttackParams)
}

.randWalkNewVector <- function(previousStrComb, modifications,
    bestAttackParams, stepChangeRatio, whileMaxCount, functionResults){
    newWalkParams <- c()
    newWalkStep <- -1
    killCount <- 0
    while (newWalkStep %in% previousStrComb ||
        length(newWalkParams[newWalkParams != (length(modifications) + 1)]) >=
            length(bestAttackParams[bestAttackParams != (length(modifications) + 1)])) {
        newWalkParamsb <- vapply(as.numeric(bestAttackParams), function(param){
            if (param != length(modifications) + 1) {
                if (sample(seq_len(round(1 / stepChangeRatio)), 1) == 1) {
                    return(sample(seq_len(length(modifications) + 1), 1))
                } else {
                    return(param)
                }
            } else {
                if (round(sum(newWalkParams ==
                        (length(modifications) + 1))^(1 / 2) /
                        stepChangeRatio) == 0){
                    return(param)
                } else {
                    if (sample(seq_len(round(sum(newWalkParams ==
                        (length(modifications) + 1))^(1 / 2) /
                        stepChangeRatio)), 1) == 1) {
                        return(sample(seq_len(length(modifications) + 1), 1))
                    } else {
                        return(param)
                    }
                }
            }
        }, numeric(1))
        newWalkParams <- newWalkParamsb
        newWalkStep <- paste(newWalkParams, collapse = " ")
        # If inifite loop kill function and return results
        if (killCount > whileMaxCount) {
            functionResults <- functionResults[order(
                    as.numeric(functionResults$genesModified)), ]
            functionResults <- functionResults[order(
                functionResults$typeModified, decreasing = TRUE), ]
            message("Inifite loop, kill function and return results")
            return(list(FALSE, functionResults))
        }
        killCount <- killCount + 1
    }
    return(list(TRUE, newWalkParams, newWalkStep))
}

.randWalkBeforeWalk <- function(exprs, genes, modifications,
            clusters, target, classifier, firstBatch, argForClassif, argForModif, verbose){
    initObjs <- .initRandWalk(modifications, genes, firstBatch)
    previousStrComb <- initObjs[[1]]
    testsGrid <- initObjs[[2]]
    functionResults <- initObjs[[3]]
    functionResultsInt <- initObjs[[4]]
    results <- initObjs[[5]]
    resultsInt <- initObjs[[6]]
    rwSeed <- .randWalkGetSeed(testsGrid, exprs, genes, modifications,
        clusters, target, classifier, results, resultsInt,
        previousStrComb, argForClassif, argForModif, verbose)
    results <- rwSeed[[1]]
    resultsInt <- rwSeed[[2]]
    previousStrComb <- rwSeed[[3]]
    results$genesModified <- length(genes) - apply(results, 1, function(x) sum(x == "NA"))
    results$typeModified <- results$prediction != target
    results$iteration <- 1
    results <- results[, c("prediction", "odd", "genesModified",
            "typeModified", "iteration", genes)]
    resultsInt$genesModified <- length(genes) - apply( resultsInt, 1,
                        function(x) sum(x == (length(modifications) + 1)))
    resultsInt$typeModified <- resultsInt$prediction != target
    resultsInt$iteration <- 1
    resultsInt <- resultsInt[, c("prediction", "odd", "genesModified",
        "typeModified", "iteration", genes)]
    functionResults <- rbind(results, functionResults)
    functionResults <- functionResults[order(functionResults$genesModified), ]
    functionResults <- functionResults[order(functionResults$typeModified, decreasing = TRUE), ]
    functionResultsInt <- rbind(resultsInt, functionResultsInt)
    functionResultsInt <- functionResultsInt[order(functionResultsInt$genesModified), ]
    functionResultsInt <- functionResultsInt[order(functionResultsInt$typeModified, decreasing = TRUE), ]
    bestAttackParams <- functionResultsInt[1, 6:ncol(functionResultsInt)]
    list(previousStrComb, functionResults, bestAttackParams)
}

.randWalkGetSeed <- function(testsGrid, exprs, genes, modifications, clusters,
    target, classifier, results, resultsInt, previousStrComb, argForClassif, argForModif, verbose){
    i<- 1
    newPreviousStrComb <- c()
    while (i <= nrow(testsGrid)) {
        if (verbose || log2(i) == round(log2(i)) || i == nrow(testsGrid)){
            message("Running first batch to determine walk seed: ", i, " on ", nrow(testsGrid))
        }
        exprsTemp <- exprs
        newPreviousStrComb <- c(newPreviousStrComb, previousStrComb[i])
        rowResults <- c()
        rowResultsInt <- c()
        for (geneInd in seq_along(genes)) {
            modifInd <- testsGrid[i, geneInd]
            if (modifInd != length(modifications) + 1) {
                mod1 <- modifications[[modifInd]][[1]]
                if (length(modifications[[modifInd]]) == 1) {
                    exprsTemp <- advModifications(exprsTemp,
                        genes[geneInd], clusters, target,
                        advMethod = mod1, verbose = verbose)
                } else {
                    mod2 <- modifications[[modifInd]][[2]]
                    exprsTemp <- advModifications(exprsTemp,
                        genes[geneInd], clusters, target,
                        advMethod = mod1,advFixedValue = mod2,
                        advFct = mod2, verbose = verbose)
                }
                rowResults <- c( rowResults,
                    paste(modifications[[modifInd]], collapse = " "))
                rowResultsInt <- c(rowResultsInt, modifInd)
            } else {
                rowResults <- c(rowResults, "NA")
                rowResultsInt <- c(rowResultsInt, modifInd)
            }
        }
        if ( argForClassif == 'SingleCellExperiment'){
            exprsTemp <- SingleCellExperiment(assays = list(counts = t(exprsTemp)))
        }
        if (!is(exprsTemp, 'data.frame') && argForClassif=='data.frame'){
            exprsTemp <- as.data.frame(exprsTemp)
        }
        if (!is(exprsTemp,'DelayedMatrix') && argForClassif=="DelayedMatrix"){
            exprsTemp <- DelayedArray::DelayedArray(exprsTemp)
        }
        classResults <- classifier(exprsTemp, clusters, target)
        rowResults <- c(rowResults, classResults[1], classResults[2])
        rowResultsInt <- c(unlist(rowResultsInt),
                    classResults[1], classResults[2])
        results <- rbind(results, rowResults)
        resultsInt <- rbind(resultsInt, rowResultsInt)
        colnames(results) <- c(genes, "prediction", "odd")
        colnames(resultsInt) <- c(genes, "prediction", "odd")
        if (classResults[1] == target) {
            i <- i + 1
        } else {
            i <- nrow(testsGrid) + 1
        }
    }
    list(results, resultsInt, newPreviousStrComb)
}

.initRandWalk <- function(modifications, genes, firstBatch){
    previousStrComb <- c(-1)
    if ((length(modifications) + 1)^length(genes) > firstBatch) {
        testsGrid <- data.frame(matrix(ncol = length(genes), nrow = 0))
        while (nrow(testsGrid) < (firstBatch - length(modifications))) {
            randSample <- sample(seq_len(length(modifications) + 1),
                length(genes),
                replace = TRUE
            )
            strRandSample <- paste(randSample, collapse = " ")
            if (!strRandSample %in% previousStrComb) {
                testsGrid <- rbind(testsGrid, randSample)
            }
            previousStrComb <- c(previousStrComb, strRandSample)
        }
        for (i in seq_along(modifications)) {
            testsGrid <- rbind(rep(i, length(genes)), testsGrid)
        }
        testsGrid <- as.matrix(testsGrid)
    } else {
        testsGrid <- gtools::permutations(n = length(modifications) + 1,
            r = length(genes), repeats.allowed = TRUE)
    }
    functionResults <- data.frame(matrix(ncol = length(genes) + 5, nrow = 0))
    colnames(functionResults) <- c("prediction", "odd", "genesModified",
        "typeModified", "iteration", genes)
    functionResultsInt <- data.frame(
        matrix(ncol = length(genes) + 5, nrow = 0)
    )
    colnames(functionResultsInt) <- c( "prediction", "odd",
        "genesModified", "typeModified", "iteration", genes)
    results <- data.frame(matrix(ncol = length(genes) + 2, nrow = 0))
    colnames(results) <- c(genes, "prediction", "odd")
    resultsInt <- data.frame(matrix(ncol = length(genes) + 2, nrow = 0))
    colnames(resultsInt) <- c(genes, "prediction", "odd")
    return (list(previousStrComb, testsGrid, functionResults,
        functionResultsInt, results, resultsInt))
}
