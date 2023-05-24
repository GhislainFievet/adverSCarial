#' Returns a modified RNA expression matrix, or a modified SingleCellExperiment,
#' for a given cluster, for a given modification.
#'
#' @details The motivation for this function is to standardize the modifications
#' we want to study in the attacks. We give as argument a matrix of the RNA
#' expression, the gene and the target cells we want to modify. Then we have
#' three arguments allowing to specify what modification we want to apply on
#' these cells. The advMethod contains, a specific prerecorded modification
#' or an indication on how to use the other two arguments.
#' The prerecorded modifications available for the advMethod argument are:
#' - 'perc1', replace the value by the whole matrix 1 percentile value of
#'  the gene. It is as if we biologically switched off the gene.
#' - 'perc99', replace the value by the whole matrix 99 percentile value of
#'  the gene. It is as if we biologically switched on the gene to the maximum.
#' The value of the advMethod argument can also be 'fixed', in this case the
#' modification would be to replace the value of the gene of the wanted cells
#' by the value of the argument 'advFixedValue'. This can be useful to test
#' aberrant values like negative integer, absurdly high values of character
#' values.
#' The value of the advMethod argument can also be 'full_row_fct',
#' 'target_row_fct', 'target_matrix_fct' or 'full_matrix_fct'. They are
#' used when we want to use a custom modification function, with the 'advFct'
#' argument:
#' - 'full_row_fct' indicate that the 'advFct' function takes the whole gene
#'  values as input.
#' - 'target_row_fct' indicate that the 'advFct' function takes target cells
#'  gene values as input.
#' - 'full_matrix_fct' indicate that the 'advFct' function takes the whole gene
#'  expression values as input.
#' - 'target_matrix_fct' indicate that the 'advFct' function takes target cells
#'  all genes values as input.
#' 
#' @param exprs can be a matrix or a data.frame of numeric RNA expression,
#' cells are rows and genes are columns. Or can be a SingleCellExperiment object.
#' @param genes the character vector of genes to modify
#' @param clusters a character vector of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param advMethod the name of the method to use
#' @param advFixedValue the numeric value to use in case of
#' advMethod=`fixed`
#' @param advFct the function to use in case advMethod
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param argForClassif the type of the first argument to feed to the
#' classifier function. 'data.frame' by default, can be 'SingleCellExperiment'
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return the matrix or a data.frame exprs modified on asked genes
#' with the specified modification
#' @examples
#' rna_expression <- data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#'
#' advModifications(rna_expression, genes, clusters_id,
#' "T cell", advMethod="perc99")
#' @export
advModifications <- function(exprs, genes, clusters,
                            target, advMethod = "perc99",
                            advFixedValue = 3, advFct = NULL,
                            argForClassif = 'data.frame',
                            verbose = FALSE) {
    if (!is(exprs, 'matrix') && !is(exprs,'data.frame') && !is(exprs,'SingleCellExperiment')){
        stop("The argument exprs must be a matrix, a data.frame or a SingleCellExperiment")
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
        message("Modify data for ", length(genes),
            " genes for cluster ", target)
    }
    if (is(exprs,'SingleCellExperiment') ){
        exprs <- as.matrix(t(counts(exprs)))
    }
    if (!is.function(advFct)) {
        exprs <- .advModificationsNotFunction(exprs,
                                    genes, clusters,
                                    target, advMethod = advMethod,
                                    advFixedValue = advFixedValue)
    } else {
        exprs <- .advModificationsFunction(exprs, genes, clusters,
                            target, advMethod = advMethod,
                            advFct = advFct)
    }
    if ( argForClassif == 'SingleCellExperiment'){
        exprs <- SingleCellExperiment(assays = list(counts = t(exprs)))
    }
    exprs
}

.advModificationsNotFunction <- function(exprs,
                            genes, clusters,
                            target, advMethod = "perc99",
                            advFixedValue = 3) {
    if (is(exprs,'SingleCellExperiment')){
        exprs <- as.matrix(t(counts(exprs)))
    }
    cellMask <- clusters == target
    if (advMethod == "fixed") {
        numGenes <- unlist(lapply(unique(genes), function(myGene){
            if (is(exprs[, myGene], "numeric")) {
                return(myGene)
            }
        }))
        exprs[cellMask, numGenes] <- advFixedValue
    }
    if (advMethod == "perc99") {
        for (myGene in unique(genes)) {
            if (is(exprs[, myGene], "numeric")) {
                exprs[cellMask, myGene] <- stats::quantile(
                    exprs[, myGene],
                    0.99
                )
            }
        }
    }
    if (advMethod == "perc1") {
        for (myGene in unique(genes)) {
            if (is(exprs[, myGene], "numeric")) {
                exprs[cellMask, myGene] <- stats::quantile(
                    exprs[, myGene],
                    0.01
                )
            }
        }
    }
    exprs
}

.advModificationsFunction <- function(exprs, genes, clusters,
                            target, advMethod = "perc99",
                            advFct = NULL) {
    cellMask <- clusters == target
    if (advMethod == "full_row_fct" || advMethod == "fixed") {
        for (myGene in unique(genes)) {
            if (is(exprs[, myGene], "numeric")) {
                exprs[cellMask, myGene] <- advFct(exprs[, myGene])
            }
        }
    }
    if (advMethod == "target_row_fct" || advMethod == "fixed") {
        for (myGene in unique(genes)) {
            if (is(exprs[, myGene], "numeric")) {
                exprs[cellMask, myGene] <- advFct(exprs[
                    cellMask,
                    myGene
                ])
            }
        }
    }
    if (advMethod == "target_matrix_fct") {
        numGenes <- unlist(lapply(unique(genes),function(myGene){
            if (is(exprs[, myGene], "numeric")) {
                return(myGene)
            } else {
                return(NULL)
            }
        }))
        exprs[cellMask, numGenes] <- advFct(exprs[cellMask, numGenes])
    }
    if (advMethod == "full_matrix_fct") {
        numGenes <- unlist(lapply(unique(genes),function(myGene){
            if (is(exprs[, myGene], "numeric")) {
                return(myGene)
            } else {
                return(NULL)
            }
        }))
        exprs[cellMask, numGenes] <- advFct(exprs[, numGenes])
    }
    exprs
}

