#' Example cell type classifier for pbmc clustered datasets.
#'
#' @details This classifier aims at testing the adverSCarial
#' package of real pbmc data. It is a simple marker based
#' classifier. It looks at the average value of a few genes
#' inside a cluster, and returns the associated cell type.
#' Markers where found by differential expressions.
#' @param expr a matrix or a data.frame of numeric
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
    if ( !is(expr, 'matrix') && !is(expr,'data.frame')){
        stop("The argument expr must be a matrix or a data.frame.")
    }
    if (!is(target,"character")) {
        stop("The argument target must be character.")
    }
    if (!is(clusters,"character")) {
        stop("The argument clusters must be a vector of character.")
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


