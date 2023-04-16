#' Example cell type classifier for pbmc clustered datasets.
#'
#' @details This classifier aims at testing the adverSCarial
#' package of real pbmc data. It is a simple marker based
#' classifier. It looks at the average value of a few genes
#' inside a cluster, and returns the associated cell type.
#' Markers where found by differential expressions.
#' @param expr a matrix, a data.frame or a DataFrame of numeric
#' RNA expression, cells are rows and genes are columns.
#' @param clusters vector of clusters to which each cell belongs
#' @param target name of the cell cluster to classify
#' @return a vector with the classification, and the odd
#' @examples
#' library(TENxPBMCData)
#'
#' pbmc <- TENxPBMCData(dataset = "pbmc3k")
#' mat_rna <- matrixFromSCE(pbmc)
#' cell_types <- system.file("extdata",
#'     "pbmc3k_cell_types.tsv",
#'     package = "adverSCarial"
#' )
#' cell_types <- read.table(cell_types, sep = "\t")$cell_type
#'
#' MClassifier(mat_rna, cell_types, "DC")
#'
#' @export
MClassifier <- function(expr, clusters, target) {
    if ( !is(expr, 'matrix') && !is(expr,'data.frame') && !is(expr,"DFrame")){
        stop("The argument expr must be a matrix, a data.frame or a DataFrame.")
    }
    if (!is(target,"character")) {
        stop("The argument target must be character.")
    }
    if (!is(clusters,"character")) {
        stop("The argument clusters must be a vector of character.")
    }
    if (is(expr, "DFrame")){
        expr <- as.data.frame(expr)
    }
    if (mean(expr[clusters == target, "LTB"]) > 7) {
        return(c("Memory CD4 T", 1))
    }
    if (mean(expr[clusters == target, "CD79A"]) > 2) {
        return(c("B", 1))
    }
    if (mean(expr[clusters == target, "S100A9"]) > 10) {
        return(c("CD14+ Mono", 1))
    }
    if (mean(expr[clusters == target, "GZMB"]) > 5) {
        return(c("NK", 1))
    }
    if (mean(expr[clusters == target, "GZMK"]) > 2) {
        return(c("CD8 T", 1))
    }
    if (mean(expr[clusters == target, "LDHB"]) > 3 &&
        mean(expr[clusters == target, "LDHB"]) < 3.3) {
        return(c("Naive CD4 T", 1))
    }
    if (mean(expr[clusters == target, "CCR7"]) > 0.5) {
        return(c("Naive CD4 T", 1))
    }

    if (mean(expr[clusters == target,"LST1"]) > 10) {
        return (c("FCGR3A+ Mono", 1))
    }
    if (mean(expr[clusters == target, "CD74"]) > 50) {
        return(c("DC", 1))
    }
    if (mean(expr[clusters == target, "PF4"]) > 10) {
        return(c("Platelet", 1))
    }
    c("UNDETERMINED", 1)
}

.get_unique_index <- function(cArray) {
    c_prev <- c()
    c_result <- c()
    for (i in seq_len(length(cArray))) {
        if (!is.na(cArray[i])) {
            if (!cArray[i] %in% c_prev) {
                c_result <- c(c_result, i)
                c_prev <- c(c_prev, cArray[i])
            }
        }
    }
    c_result
}

#' Returns the RNA expression matrix from a SingleCellExperiment
#' with unique hgnc gene names in columns
#'
#' @details This function retrieves from a SingleCellExperiment
#' object the raw RNA expression value corresponding to the hgnc
#' genes. The resulting matrix can then be used with adverSCarial
#' packages.
#' @param sce SingleCellExperiment object to convert
#' @return the RNA expression matrix from a SingleCellExperiment
#' with unique hgnc gene names in columns
#' @examples
#' library(TENxPBMCData)
#'
#' pbmc <- TENxPBMCData(dataset = "pbmc3k")
#' mat_rna <- matrixFromSCE(pbmc)
#'
#' @export
matrixFromSCE <- function(sce) {
    if ( !is(sce,'SingleCellExperiment')){
        stop("The argument sce must be a SingleCellExperiment.")
    }
    ind2keep <- .get_unique_index(
        #sce@rowRanges@elementMetadata@listData$Symbol_TENx
        sce@rowRanges@elementMetadata@listData$Symbol
    )
    df <- data.frame(sce@assays@data$counts)
    df <- df[ind2keep, ]
    rownames(df) <- sce@ rowRanges@ elementMetadata@listData$Symbol[ind2keep]
    colnames(df) <- colData(sce)$Barcode
    t(df)
}

#' Returns a SingleCellExperiment object keeping unique HGNC gene
#'
#' @details Sometimes classifiers need HGNC instead of ensemble genes
#' to run. This function allows to make the conversion.
#' @param sce SingleCellExperiment object to convert
#' @return the SingleCellExperiment object keeping unique HGNC gene
#' @examples
#' library(TENxPBMCData)
#'
#' pbmc <- TENxPBMCData(dataset = "pbmc3k")
#' hgnc_pbmc <- sceConvertToHGNC(pbmc)
#'
#' @export
sceConvertToHGNC <- function(sce){
    if ( !is(sce,'SingleCellExperiment')){
        stop("The argument sce must be a SingleCellExperiment.")
    }
    mat_sce <- matrixFromSCE(sce)
    SingleCellExperiment(assays = list(counts = t(mat_sce)),
        colData = colData(sce))
}