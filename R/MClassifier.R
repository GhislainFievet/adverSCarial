#' Example cell type classifier for pbmc clustered datasets.
#'
#' @details This classifier aims at testing the adverSCarial
#' package of real pbmc data. It is a simple marker based
#' classifier. It looks at the average value of a few genes
#' inside a cluster, and returns the associated cell type.
#' Markers where found by differential expressions.
#' @param exprs DelayedMatrix of numeric RNA expression, cells are rows and genes
#' are columns - or a SingleCellExperiment object, a matrix or a data.frame.
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
MClassifier <- function(exprs, clusters, target) {
    if (!is(exprs, 'matrix') && !is(exprs,'data.frame') &&
        !is(exprs,'SingleCellExperiment') && !is(exprs,'DelayedMatrix')){
        stop("The argument exprs must be a DelayedMatrix, a SingleCellExperiment, a matrix or a data.frame")
    }
    if (is(exprs,'SingleCellExperiment') ){
        exprs <- t(counts(exprs))
    }

    res_cell_type <- "UNDETERMINED"
    if (!is(target,"character")) {
        stop("The argument target must be character.")
    }
    if (!is(clusters,"character")) {
        stop("The argument clusters must be a vector of character.")
    }
    if (mean(exprs[clusters == target, "LTB"]) > 7) {
        res_cell_type <- "Memory CD4 T"
    }
    if (mean(exprs[clusters == target, "CD79A"]) > 2) {
        res_cell_type <- "B"
    }
    if (mean(exprs[clusters == target, "S100A9"]) > 10) {
        res_cell_type <- "CD14+ Mono"
    }
    if (mean(exprs[clusters == target, "GZMB"]) > 5) {
        res_cell_type <- "NK"
    }
    if (mean(exprs[clusters == target, "GZMK"]) > 2) {
        res_cell_type <- "CD8 T"
    }
    if (mean(exprs[clusters == target, "LDHB"]) > 3 &&
        mean(exprs[clusters == target, "LDHB"]) < 3.3) {
        res_cell_type <- "Naive CD4 T"
    }
    if (mean(exprs[clusters == target, "CCR7"]) > 0.5) {
        res_cell_type <- "Naive CD4 T"
    }

    if (mean(exprs[clusters == target,"LST1"]) > 10) {
        res_cell_type <- "FCGR3A+ Mono"
    }
    if (mean(exprs[clusters == target, "CD74"]) > 50) {
        res_cell_type <- "DC"
    }
    if (mean(exprs[clusters == target, "PF4"]) > 10) {
        res_cell_type <- "Platelet"
    }

    df_tp <- as.data.frame(matrix(0,
        nrow=nrow(exprs), ncol=length(unique(clusters))))
    colnames(df_tp) <- unique(clusters)
    rownames(df_tp) <- rownames(exprs)
    df_tp[,res_cell_type] <- 1
    df_tp <- t(df_tp)

    resClass <- list(prediction=res_cell_type, odd=1,
                      typePredictions=df_tp,
                      cellTypes=rep(res_cell_type, sum(clusters == target)))
    return(resClass)
}


