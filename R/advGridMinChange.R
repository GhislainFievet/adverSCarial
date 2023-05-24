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
#' @param genes the character vector of genes to study
#' @param modifications the list of the modifications to study
#' @param returnFirstFound set to TRUE to return result when a
#' the first misclassification is found
#' @param argForClassif the type of the first argument to feed to the
#' classifier function. 'data.frame' by default, can be 'SingleCellExperiment'
#' @param verbose logical, set to TRUE to activate verbose mode
#' @param iamsure logical, prevents from expansive calculations
#' when `genes` list is too long, set to `TRUE` to run anyway.
#' @return DataFrame results of the classification of all the grid combinations
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
                returnFirstFound = FALSE, argForClassif = 'data.frame',
                verbose = FALSE, iamsure = FALSE) {
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
    if (!is.character(genes)) {
        stop("The argument genes must be character or vector of character.")
    }
    if (!is.list(modifications)) {
        stop("The argument modifications must be a list.")
    }
    if (!is.logical(returnFirstFound)){
        stop("The argument returnFirstFound must be logical.")
    }
    if (!is.character(argForClassif)) {
        stop("The argument argForClassif must be character: 'data.frame' or 'SingleCellExperiment'.")
    }
    if (!is.logical(verbose)){
        stop("The argument verbose must be logical.")
    }
    if (!is.logical(iamsure)){
        stop("The argument iamsure must be logical.")
    }

    if (is(exprs,'SingleCellExperiment') ){
        exprs <- as.matrix(t(counts(exprs)))
    }

    functionResults <- data.frame(matrix(ncol = length(genes) + 5, nrow = 0))
    colnames(functionResults) <- c("prediction", "odd", "genesModified",
        "typeModified", "iteration", genes)
    if (.gridWarning(modifications, genes, iamsure))
        return(functionResults)
    testsGrid <- gtools::permutations(n = length(modifications) + 1,
        r = length(genes), repeats.allowed = TRUE)
    results <- data.frame(matrix(ncol = length(genes) + 2, nrow = 0))
    colnames(results) <- c(genes, "prediction", "odd"); i <- 1
    while (i <= nrow(testsGrid)) {
        message("Running combination: ", i, " on ", nrow(testsGrid))
        exprsTemp <- exprs
        rowResults <- c()
        for (geneInd in seq_along(genes)) {
            modifInd <- testsGrid[i, geneInd]
            if (modifInd != length(modifications) + 1) {
                mod1 <- modifications[[modifInd]][[1]]
                if (length(modifications[[modifInd]]) == 1) {
                    exprsTemp <- advModifications(exprsTemp, genes[geneInd],
                        clusters, target, advMethod = mod1, verbose = verbose)
                } else {
                    mod2 <- modifications[[modifInd]][[2]]
                    exprsTemp <- advModifications(exprsTemp, genes[geneInd],
                        clusters, target, advMethod = mod1, advFct = mod2,
                        advFixedValue = mod2, verbose = verbose)}
                    rowResults <- c(rowResults,
                        paste(modifications[[modifInd]], collapse = " "))
            } else { rowResults <- c(rowResults, "NA")}
        }
        if ( argForClassif == 'SingleCellExperiment'){
            exprsTemp <- SingleCellExperiment(assays = list(counts = t(exprsTemp)))
        }
        classResults <- classifier(exprsTemp, clusters, target)
        rowResults <- c(rowResults, classResults[1], classResults[2])
        results <- rbind(results, rowResults)
        colnames(results) <- c(genes, "prediction", "odd")
        ifelse(classResults[1] == target, i <- i + 1, i <- nrow(testsGrid)+1)
    }
    results$genesModified <- length(genes) - apply(results, 1, function(x) sum(x == "NA"))
    results$typeModified <- results$prediction != target
    results <- results[, c( "prediction", "odd",
            "genesModified", "typeModified", genes)]
    functionResults <- rbind(results, functionResults)
    functionResults <-
        functionResults[order(as.numeric(functionResults$genesModified)), ]
    functionResults <- functionResults[order(functionResults$typeModified,
            decreasing = TRUE), ]
    S4Vectors::DataFrame(functionResults)
}
