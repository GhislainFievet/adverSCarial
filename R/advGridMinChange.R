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
