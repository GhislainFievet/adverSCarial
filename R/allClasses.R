#' adverSCarial class
#' 
#' \code{advChar} is a class used to store the output values
#' of the \code{advMaxChange} function. The result can be a vector
#' of few thousands genes, so a specific *show* method is associated to
#' this class to avoid overflooding the R scripts outputs.
#' @rdname advChar
#' @return A advChar object
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'   c("T cell", 0.9)
#' }
#'
#' genes <- paste0("gene_",1:10000)
#' rna_expression <- data.frame(lapply(genes, function(x) numeric(0)))
#' rna_expression <- rbind(rna_expression, rep(1,10000))
#' rna_expression <- rbind(rna_expression, rep(2,10000))
#' colnames(rna_expression) <- genes
#' clusters_id <- c("B cell","T cell")
#'
#' max_change_genes <- advMaxChange(rna_expression, clusters_id,
#' "T cell", MyClassifier, advMethod="perc99")
#'                                    
#' max_change_genes
#' @export
advChar <- methods::setClass("advChar",
        slots = list(
        values = "character"
        ))
setMethod("show", "advChar",
        function(object) {
            if ( length(object@values)<10){
                print(object@values)
            } else {
                print(paste0("Vector with ", length(object@values), " values:"))
                if (is.null(names(object@values))){
                    print(paste0(paste(object@values[1:10], collapse=", "), " ..."))
                } else {
                    print(paste0(paste(paste0(names(object@values[1:5]),": ",object@values[1:5]),
                                        collapse=", "), " ..."))
                }
            }
        })

#' adverSCarial class
#' 
#' \code{advList} is a class used to store the output values
#' of the \code{advSingleGene} function. The result can be a list
#' of few thousands genes:cell_type key/values, so a specific *show*
#' method is associated to this class to avoid overflooding the R
#' scripts outputs.
#' @rdname advList
#' @return A advList object
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'    c("B cell", 0.9)
#' }
#' rna_expression <- data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'      CD8B=c(2,2,3,3))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","B cell","T cell","T cell")
#' 
#' adv_min_change <- advSingleGene(rna_expression, clusters_id,
#' "T cell", MyClassifier, advMethod="perc99")
#' 
#' adv_min_change
#' 
#' @export
advList <- methods::setClass("advList",
        slots = list(
        values = "list"
        ))
setMethod("show", "advList",
        function(object) {
            if ( length(object@values)<10){
                print(object@values)
            } else {
                print(paste0("Vector with ", length(object@values), " values:"))
                if (is.null(names(object@values))){
                    print(paste0(paste(object@values[1:10], collapse=", "), " ..."))
                } else {
                    print(paste0(paste(paste0(names(object@values[1:5]),": ",object@values[1:5]),
                                        collapse=", "), " ..."))
                }
            }
        })
