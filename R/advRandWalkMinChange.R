
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

