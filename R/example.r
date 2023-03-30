#' Example cell type classifier for the pbmc3k dataset
#'
#' @param expr RNA expression matrix
#' @param clusters list of clusters to which each cell belongs
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
    if (mean(expr[clusters == target, "LDHB"]) > 17 &&
        mean(expr[clusters == target, "LDHB"]) < 20) {
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
    c("undetermined", 1)
}

.get_unique_index <- function(c_array) {
    c_prev <- c()
    c_result <- c()
    for (i in seq_len(length(c_array))) {
        if (!is.na(c_array[i])) {
            if (!c_array[i] %in% c_prev) {
                c_result <- c(c_result, i)
                c_prev <- c(c_prev, c_array[i])
            }
        }
    }
    c_result
}

#' Returns the RNA expression matrix from a SingleCellExperiment
#' with unique hgnc gene names in columns
#'
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
    mat_sce <- matrixFromSCE(sce)
    SingleCellExperiment(assays = list(counts = t(mat_sce)),
        colData = colData(sce))
}