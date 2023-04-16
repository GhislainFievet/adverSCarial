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
