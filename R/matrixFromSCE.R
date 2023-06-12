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
    # symbols <- sce@rowRanges@elementMetadata@listData$Symbol_TENx
    symbols <- sce@rowRanges@elementMetadata@listData$Symbol
    ind2keep <- setdiff(match(unique(symbols), symbols), which(is.na(symbols)))

    df <- data.frame(sce@assays@data$counts)
    df <- df[ind2keep, ]
    rownames(df) <- sce@ rowRanges@ elementMetadata@listData$Symbol[ind2keep]
    colnames(df) <- colData(sce)$Barcode
    t(df)
}