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
    if ( !is(exprs, 'matrix') && !is(exprs,'data.frame') && !is(exprs,"DFrame")){
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

