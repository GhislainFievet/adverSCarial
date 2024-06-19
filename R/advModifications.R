#' Returns a modified RNA expression DelayedMatrix, or a modified SingleCellExperiment,
#' for a given cluster, for a given modification.
#'
#' @details The motivation for this function is to standardize the modifications
#' we want to study in the attacks. We give as argument a DelayedMatrix of the RNA
#' expression, the gene and the target cells we want to modify. Then we have
#' three arguments allowing to specify what modification we want to apply on
#' these cells. The advMethod contains, a specific prerecorded modification
#' or an indication on how to use the other two arguments.
#' The prerecorded modifications available for the advMethod argument are:
#' - 'perc1', replace the value by the whole matrix 1 percentile value of
#'  the gene. It is as if we biologically switched off the gene.
#' - 'perc99', replace the value by the whole matrix 99 percentile value of
#'  the gene. It is as if we biologically switched on the gene to the maximum.
#' - 'random', replace the value by from a uniform distribution between min 
#' and max of the gene on the dataset
#' - 'positive_aberrant' replace value by 10,000 times the max value of the
#' gene on the dataset
#' - 'negative_aberrant' replace value by -10,000 times the max value of the
#' gene on the dataset
#' - 'decile+X', shifts the gene value by + X deciles.
#' - 'decile-X', shifts the gene value by - X deciles.
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
#' @param exprs DelayedMatrix of numeric RNA expression, cells are rows and genes
#' are columns - or a SingleCellExperiment object, a matrix or a data.frame. By default,
#' these are converted to a data.frame to increase speed performance during modifications.
#' However, this conversion can consume a significant amount of memory, see 'argForModif'
#' argument for options.
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
#' classifier function. 'DelayedMatrix' by default, can be 'SingleCellExperiment'
#' or 'data.frame'.
#' @param argForModif type of matrix during for the modification, 'DelayedMatrix'
#' by default. Can be 'data.frame', which is faster, but need more memory.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return the matrix or a data.frame exprs modified on asked genes
#' with the specified modification
#' @examples
#' library(DelayedArray)
#' 
#' rna_expression <- DelayedArray(data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3)))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#'
#' advModifications(rna_expression, genes, clusters_id,
#' "T cell", advMethod="perc99")
#' @export
advModifications <- function(exprs, genes, clusters,
                            target, advMethod = "perc99",
                            advFixedValue = 3, advFct = NULL,
                            argForClassif = 'DelayedMatrix',
                            argForModif = 'data.frame',
                            slot=NULL,
                            verbose = FALSE) {
    if (!is(exprs, 'matrix') && !is(exprs,'data.frame') &&
        !is(exprs,'SingleCellExperiment') && !is(exprs,'DelayedMatrix') &&
        !is(exprs,'Seurat')){
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
    # if (is(exprs,'Seurat') && tolower(argForModif)=="seurat"){
    if (is(exprs,'Seurat')){
        so_exprs <- exprs
        exprs <- as.data.frame(so_exprs@assays$RNA@layers[[slot]])
        colnames(exprs) <- rownames(so_exprs@assays$RNA@cells)
        rownames(exprs) <- rownames(so_exprs@assays$RNA@features)
        exprs <- as.data.frame(t(exprs))
    }
    if (is(exprs,'SingleCellExperiment') ){
        exprs <- t(as.matrix(counts(exprs)))
    }
    if (!is(exprs,'DelayedMatrix') && argForModif=="DelayedMatrix"){
        if (verbose){
            message("Converting exprs object to a DelayedArray object")
        }
        exprs <- DelayedArray::DelayedArray(exprs)
    }


    if (!is(exprs,'data.frame') && argForModif=="data.frame"){
        if (verbose){
            message("Converting exprs object to a data.frame object")
        }
        exprs <- as.data.frame(exprs)
    }

    
    if (advMethod != "none"){
        if (advMethod %in% c("random", "positive_aberrant", "negative_aberrant")){
            if (advMethod == "random"){
                exprs <- .advModificationsFunction(exprs, genes, clusters,
                                    target, advMethod = "full_matrix_fct",
                                    advFct = .fctRand)
            }
            if (advMethod == "positive_aberrant"){
                exprs <- .advModificationsFunction(exprs, genes, clusters,
                                    target, advMethod = "full_row_fct",
                                    advFct = .fctAbbPos)
            }
            if (advMethod == "negative_aberrant"){
                exprs <- .advModificationsFunction(exprs, genes, clusters,
                                    target, advMethod = "full_row_fct",
                                    advFct = .fctAbbNeg)
            }
        } else {
            if (startsWith(advMethod, "decile")){
                if ( unlist(strsplit(advMethod, "\\+"))[1] == "decile" ){
                    decilNumber <- as.numeric(unlist(strsplit(advMethod, "\\+")))[2]
                    exprs <- .advModificationsFunction(exprs, genes, clusters,
                                target, advMethod = "full_row_fct",
                                advFct = function(x,y){.advDecil(x,y,decilNumber=decilNumber)})
                }
                if ( unlist(strsplit(advMethod, "\\-"))[1] == "decile" ){
                    decilNumber <- - as.numeric(unlist(strsplit(advMethod, "\\-")))[2]
                    exprs <- .advModificationsFunction(exprs, genes, clusters,
                                target, advMethod = "full_row_fct",
                                advFct = function(x,y){.advDecil(x,y,decilNumber=decilNumber)})
                }
            } else {
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
            }
        }


        if(tolower(argForClassif)=="seurat"){
            so_exprs@assays$RNA@layers[[slot]] <- t(exprs)
            exprs <- so_exprs
        }
    } else {
        if(tolower(argForClassif)=="seurat"){
            exprs <- so_exprs
        }
    }

    if ( argForClassif == 'SingleCellExperiment'){
        exprs <- SingleCellExperiment(assays = list(counts = t(exprs)))
    }
    if (!is(exprs, 'data.frame') && argForClassif=='data.frame'){
        exprs <- as.data.frame(exprs)
    }
    if (!is(exprs,'DelayedMatrix') && argForClassif=="DelayedMatrix"){
        exprs <- DelayedArray::DelayedArray(exprs)
    }

    exprs
}

# Function for the default "random" modification
.fctRand <- function(x, cMask){
    sample(min(x):max(x), sum(cMask), replace=T)
}

# Function for the default "positive_aberrant" modification
.fctAbbPos <- function(x, y){10000*max(x)}

# Function for the default "negative_aberrant" modification
.fctAbbNeg <- function(x, y){-10000*max(x)}

.advModificationsNotFunction <- function(exprs,
                            genes, clusters,
                            target, advMethod = "perc99",
                            advFixedValue = 3) {
    if (is(exprs,'SingleCellExperiment')){
        exprs <- t(as.matrix(counts(exprs)))
    }
    cellMask <- clusters == target
    if (advMethod == "fixed") {
        numGenes <- unlist(lapply(unique(genes), function(myGene){
            if (is(exprs[, myGene], "numeric")) {
                return(myGene)
            }
        }))
        exprs[cellMask, numGenes] <- as.numeric(advFixedValue)
    }
    if (advMethod == "perc99") {
        myGenes <- unlist(lapply(unique(genes), function(myGene){
            if(is(exprs[, myGene], "numeric")){
                return(myGene)
            } else {
                return(NULL)
            }
        }))
        myQuantiles <- unlist(lapply(myGenes, function(myGene){
            stats::quantile(exprs[, myGene], 0.99)
        }))
        exprs[cellMask, myGenes] <- rep(myQuantiles, each=sum(cellMask))
    }
    if (advMethod == "perc1") {
        myGenes <- unlist(lapply(unique(genes), function(myGene){
            if(is(exprs[, myGene], "numeric")){
                return(myGene)
            } else {
                return(NULL)
            }
        }))
        myQuantiles <- unlist(lapply(myGenes, function(myGene){
            stats::quantile(exprs[, myGene], 0.01)
        }))
        exprs[cellMask, myGenes] <- rep(myQuantiles, each=sum(cellMask))
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
                exprs[cellMask, myGene] <- advFct(exprs[, myGene], cellMask)
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
        exprs[cellMask, numGenes] <- advFct(exprs[, numGenes], cellMask)
    }
    exprs
}

.advDecil <- function(x, y, decilNumber=1){
    orderedX <- x[order(x)]
    modValues <- unlist(lapply(x[y], function(val){
        xRank <- sum(orderedX < val)
        if ( decilNumber > 0){
            if ( xRank + round(decilNumber*length(x) / 10) > length(x) ){
                return(max(x))
            } else {
                return(orderedX[xRank + round(decilNumber*length(x) / 10)])
            }
        } else {
            if ( xRank + round(decilNumber*length(x) / 10) < 1 ){
                return(min(x))
            } else {
                return(orderedX[xRank + round(decilNumber*length(x) / 10)])
            }
        }
    }))
    modValues
}

