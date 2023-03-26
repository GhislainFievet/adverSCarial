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
#' RFClassifier(mat_rna, cell_types, "DC")
#'
#' @export
RFClassifier <- function(expr, clusters, target) {
    set.seed(20)
    colnames(expr) <- stringr::str_replace_all(colnames(expr), "-", "_")
    colnames(expr) <- stringr::str_replace(colnames(expr), "^", "g_")

    if (!exists("rf_scrnaseq")) {
        load(system.file("extdata", "rf_scrnaseq", package = "adverSCarial"))
        rf_scrnaseq <<- rf
    }

    rf_features <- names(rf_scrnaseq$forest$xlevels)
    c_diff_genes <- setdiff(rf_features, colnames(expr))
    expr <- as.data.frame(expr)
    expr[, c_diff_genes] <- 0

    final_predictions <- predict(
                    rf_scrnaseq, expr[clusters == target, ])
    ratio <- as.numeric(sort(table(final_predictions),
        decreasing = T
    )[1]) /
        sum(as.numeric(sort(table(final_predictions), decreasing = T)))
    predicted_class <- names(sort(table(final_predictions), decreasing = T)[1])
    if (ratio < 0.5) {
        predicted_class <- "NA"
    }
    c(predicted_class, ratio)
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
    #rownames(df) <- sce@ rowRanges@ elementMetadata@listData$Symbol_TENx[ind2keep]
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