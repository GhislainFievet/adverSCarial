#' Returns a modified RNA expression matrix, for a given cluster,
#' for a given modification.
#'
#' @param exprs a matrix or dataframe of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param genes the list of genes to modify
#' @param clusters a list of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param adv_method the name of the method to use
#' @param adv_fixed_value the numeric value to use in case of
#' adv_method=`fixed`
#' @param adv_fct the function to use in case adv_method
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return the matrix or a dataframe exprs modified on asked genes
#' with the specified modification
#' @examples
#' advModifications(rna_expression, c("CD4", "CD8"), clusters_id,
#' "Memory CD4+", adv_method="perc99")
#' @export
advModifications <- function(exprs, genes, clusters,
                            target, adv_method = "perc99",
                            adv_fixed_value = 3, adv_fct = NULL,
                            verbose = FALSE) {
    if (verbose) {
        message("Modify data for ", length(genes),
            " genes for cluster ", target)
    }

    if (is.null(adv_fct)) {
        exprs <- .advModificationsNotFunction(exprs,
                                    genes, clusters,
                                    target, adv_method = adv_method,
                                    adv_fixed_value = adv_fixed_value)
    } else {
        exprs <- .advModificationsFunction(exprs, genes, clusters,
                            target, adv_method = adv_method,
                            adv_fct = adv_fct)
    }

    exprs
}

.advModificationsNotFunction <- function(exprs,
                            genes, clusters,
                            target, adv_method = "perc99",
                            adv_fixed_value = 3) {
        cell_mask <- clusters == target
        if (adv_method == "fixed") {
            num_genes <- c()
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    num_genes <- append(num_genes, my_gene)
                }
            }
            exprs[cell_mask, num_genes] <- adv_fixed_value
        }
        if (adv_method == "perc99") {
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    exprs[cell_mask, my_gene] <- stats::quantile(
                        exprs[, my_gene],
                        0.99
                    )
                }
            }
        }
        if (adv_method == "perc1") {
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    exprs[cell_mask, my_gene] <- stats::quantile(
                        exprs[, my_gene],
                        0.01
                    )
                }
            }
        }
        exprs
}

.advModificationsFunction <- function(exprs, genes, clusters,
                            target, adv_method = "perc99",
                            adv_fct = NULL) {
        cell_mask <- clusters == target
        if (adv_method == "full_row_fct" || adv_method == "fixed") {
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    exprs[cell_mask, my_gene] <- adv_fct(exprs[, my_gene])
                }
            }
        }
        if (adv_method == "target_row_fct" || adv_method == "fixed") {
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    exprs[cell_mask, my_gene] <- adv_fct(exprs[
                        cell_mask,
                        my_gene
                    ])
                }
            }
        }
        if (adv_method == "target_matrix_fct") {
            num_genes <- c()
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    num_genes <- append(num_genes, my_gene)
                }
            }
            exprs[cell_mask, num_genes] <- adv_fct(exprs[cell_mask, num_genes])
        }
        if (adv_method == "full_matrix_fct") {
            num_genes <- c()
            for (my_gene in unique(genes)) {
                if (is(exprs[, my_gene], "numeric")) {
                    num_genes <- append(num_genes, my_gene)
                }
            }
            exprs[cell_mask, num_genes] <- adv_fct(exprs[, num_genes])
        }
    exprs
}

#' Returns a classification and an odd value from a
#' RNA expression matrix, for given genes, for a given cluster,
#' for a given modification.
#'
#' @param exprs a matrix or dataframe of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param genes the list of genes to modify
#' @param clusters a list of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format
#' @param adv_method the name of the method to use
#' @param adv_fixed_value the numeric value to use in case of
#' adv_method=`fixed`
#' @param adv_fct the function to use in case adv_method
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a vector of the classification, and the associated odd
#' @examples
#' predictWithNewValue(rna_expression, c("CD4", "CD8"), clusters_id,
#' "Memory CD4+", myClassifier, adv_method="perc99")
#' @export
predictWithNewValue <- function(exprs, genes, clusters, target,
                                classifier, adv_method = "perc99",
                                adv_fixed_value = 3,
                                adv_fct = NULL, verbose = FALSE) {
    modif_exprs <- advModifications(exprs, genes, clusters,
        target,
        adv_method = adv_method,
        adv_fixed_value = adv_fixed_value,
        adv_fct = adv_fct
    )

    classifier(modif_exprs, clusters, target)
}

#' Find a max change adversarial attack. It finds the longer
#' list of genes you can modify on a cluster without changing its
#' classification.
#' 
#' @param exprs a matrix or dataframe of numeric RNA expression,
#' cells are rows and genes are columns.

#' @param clusters a list of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format
#' @param excl_genes a list of genes to exclude from the analysis
#' @param genes a list of genes in case you want to limit the
#' attack on a subset of genes
#' @param adv_method the name of the method to use
#' @param adv_fixed_value the numeric value to use in case of
#' adv_method=`fixed`
#' @param adv_fct the function to use in case adv_method
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param max_split_size max size of dichotomic slices.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a list of genes you can modify on a cluster without
#' modifying its classification
#' @examples
#' advMaxChange(rna_expression, clusters_id,
#' "Memory CD4+", myClassifier, adv_method="perc99")
#' @export
advMaxChange <- function(exprs, clusters, target, classifier,
                        excl_genes = c(), genes = c(), adv_method = "perc99",
                        adv_fixed_value = 3, adv_fct = NULL,
                        max_split_size = 1, verbose = FALSE) {
    genes_to_keep <- c()
    for (str_gene in unique(genes)) {
        if (!is.na(match(str_gene, colnames(exprs)))) {
            genes_to_keep <- append(genes_to_keep, str_gene)
        }
    }
    if (length(genes) == 0) {
        genes_to_keep <- colnames(exprs)
    }
    for (str_gene in unique(excl_genes)) {
        if (!is.na(match(str_gene, genes_to_keep))) {
            genes_to_keep <- genes_to_keep[-match(str_gene, genes_to_keep)]
        }
    }
    l_results <- .predictionDichotMaxChange(exprs, genes_to_keep,
        clusters, target, classifier,
        adv_method = adv_method,
        adv_fixed_value = adv_fixed_value, adv_fct = adv_fct,
        max_split_size = max_split_size, verbose = verbose
    )
    l_results
}

.dichotMaxSameType <- function(l_results, c_gene_split_value, exprs,
                clusters, target, classifier, adv_method,
                adv_fixed_value, adv_fct, max_split_size, verbose = TRUE){
        if (verbose) { message("same cell_type") }
        if (length(l_results) == 0) {
            l_results <- c_gene_split_value
        } else {
            if (verbose){
                message("check if concat lists still gives target")
            }
            concat_genes_list <- c(l_results, c_gene_split_value)
            if (verbose) {
                message(paste0(
                    "length of tested gene list: ",
                    length(unique(concat_genes_list))
                ))
            }
            # check if concat lists still gives target
            mod_obj <- predictWithNewValue(exprs, concat_genes_list, clusters,
                target, classifier,
                adv_method = adv_method,
                adv_fixed_value = adv_fixed_value, adv_fct = adv_fct,
                verbose = verbose
            )
            concat_cell_type <- mod_obj[1]
            if (concat_cell_type == target) {
                if (verbose) { message("YES: merge results") }
                # Merge results
                l_results <- c(l_results, c_gene_split_value)
                if (verbose) { message(paste0( "results length after merge: ",
                        length(l_results)
                    ))
                }
            } else {
                if (verbose) { message("NO: split and retry") }
                if (length(c_gene_split_value) > max_split_size) {
                    l_results <- .predictionDichotMaxChange(exprs,
                        c_gene_split_value, clusters, target,
                        classifier,
                        l_results = l_results,
                        adv_method = adv_method,
                        adv_fixed_value = adv_fixed_value, adv_fct = adv_fct,
                        max_split_size = max_split_size, verbose = verbose
                    )
                }
            }
        }
    l_results
}

.dichotMaxSplit <- function(l_results, c_gene_split_value, exprs,
                    genes, clusters, target,
                                    classifier,
                                    adv_method,
                                    adv_fixed_value,
                                    adv_fct, max_split_size,
                                    verbose){
    if (length(c_gene_split_value) != 0) {
        if (verbose) {
            message("before predictWithNewValue 1")
            message(length(c_gene_split_value))
        }
        mod_obj <- predictWithNewValue(exprs, c_gene_split_value,
            clusters, target, classifier,
            adv_method = adv_method,
            adv_fixed_value = adv_fixed_value, adv_fct = adv_fct,
            verbose = TRUE
        )
        cell_type <- mod_obj[1]
        celltype_score <- mod_obj[2]
        message(cell_type)
        message(celltype_score)
    }
    if (cell_type == target) {
        l_results <- .dichotMaxSameType(l_results,
                    c_gene_split_value,
                    exprs, clusters,
                    target, classifier,
                    adv_method,
                    adv_fixed_value, adv_fct,
                    max_split_size,
                    verbose = TRUE)
    } else {
        if (verbose) { message("NOT same cell_type") }
        if (length(c_gene_split_value) > max_split_size) {
            l_results <- .predictionDichotMaxChange(exprs,
                unname(c_gene_split_value), clusters, target,
                classifier, l_results = l_results, adv_method = adv_method,
                adv_fixed_value = adv_fixed_value, adv_fct = adv_fct,
                max_split_size = max_split_size, verbose = verbose
            )
        }
    }
    l_results
}


.predictionDichotMaxChange <- function(exprs, genes, clusters, target,
                                    classifier, l_results = c(),
                                    adv_method = "perc99",
                                    adv_fixed_value = 3,
                                    adv_fct = NULL, max_split_size = 1,
                                    verbose = TRUE) {
    c_gene_split <- split(unlist(genes), 1:2)
    message(paste0("genes size: ", length(unlist(genes))))
    message(paste0("current gene results length: ", length(l_results)))
    # Random order for each reccursion
    random_index <- sample(1:2, 2)
    c_gene_split <- c_gene_split[random_index]
    l_results <- .dichotMaxSplit(l_results, unname(unlist(c_gene_split[1])),
                    exprs, genes, clusters, target, classifier, adv_method,
                    adv_fixed_value, adv_fct, max_split_size, verbose)
    l_results <- .dichotMaxSplit(l_results, unname(unlist(c_gene_split[2])),
                exprs, genes, clusters, target, classifier, adv_method,
                adv_fixed_value, adv_fct, max_split_size, verbose)
    l_results
}

#' Find a one gene min change adversarial attack list. A one gene min change
#' adversarial attack  refers to the modification of a single gene
#' within a cluster, leading to a change in its classification. The function
#' returns a list of genes/new classification.
#'
#' @param exprs a matrix or dataframe of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param clusters a list of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format
#' @param excl_genes a list of genes to exclude from the analysis
#' @param genes a list of genes in case you want to limit the
#' attack on a subset of genes
#' @param adv_method the name of the method to use
#' @param adv_fixed_value the numeric value to use in case of
#' adv_method=`fixed`
#' @param adv_fct the function to use in case adv_method
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param first_dichot the initial number of slices before
#' the dichotomic search
#' @param max_split_size max size of dichotomic slices
#' @param return_first_found set to TRUE to return result when a
#' the first misclassification is found
#' @param change_type `any` consider each misclassification,
#'  `not_na` consider each misclassification but NA.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a list of genes/new classification tuples
#' @examples
#' advMinChange(rna_expression, clusters_id,
#' "Memory CD4+", myClassifier, adv_method="perc99")
#' @export
advMinChange <- function(exprs, clusters, target, classifier, excl_genes = c(),
        genes = c(), adv_method = "perc99", adv_fixed_value = 3,
        adv_fct = NULL, first_dichot = 100, max_split_size = 1,
        return_first_found = FALSE, change_type = "any", verbose = FALSE) {
    genes_index <- seq_len(ncol(exprs))
    if (length(excl_genes) != 0 && length(genes) == 0) {
        genes_to_remove <- c()
        for (str_gene in excl_genes) {
            if (!is.na(match(str_gene, colnames(exprs)))) {
                genes_to_remove <- append(genes_to_remove,
                    -match(str_gene, colnames(exprs)))}
        }
        if (length(genes_to_remove) != 0)
            genes_index <- genes_index[genes_to_remove]
    }
    if (length(genes) != 0) {
        excl_genes <- unique(excl_genes)
        genes <- unique(genes)
        genes_to_keep <- genes[!genes %in% excl_genes]
        ind_genes_to_keep <- c()
        for (str_gene in genes_to_keep) {
            if (!is.na(match(str_gene, colnames(exprs)))) {
                ind_genes_to_keep <- append(ind_genes_to_keep,
                    match(str_gene, colnames(exprs)))
            }
        }
        if (length(ind_genes_to_keep) != 0)
            genes_index <- genes_index[ind_genes_to_keep]
    }
    c_split_index <- split(genes_index, 1:first_dichot)
    l_splits_results <- c(); i <- 1
    while(i <= first_dichot) {
        prev_time <- Sys.time()
        message(paste0("Split number: ", i, "/", first_dichot))
        if (length(unlist(c_split_index[i])) != 0) {
            l_splits_results <- c( l_splits_results,
                .predictionDichotMinChange(exprs,
                    colnames(exprs)[unlist(c_split_index[i])],
                    clusters, target, classifier, adv_method = adv_method,
                    adv_fixed_value = adv_fixed_value, adv_fct = adv_fct,
                    max_split_size = max_split_size,
                    return_first_found = return_first_found,
                    change_type=change_type,verbose = verbose))
        }
        message(paste0("Split time: ", Sys.time() - prev_time))
        ifelse(length(l_splits_results) > 0 && return_first_found,
            i <- first_dichot + 1, i <- i + 1)
    }
    l_splits_results
}

.dichotMinSplit <- function(l_results, c_gene_split_value, exprs,
            genes, clusters, target, classifier, adv_method, adv_fixed_value,
            adv_fct, max_split_size, change_type, return_first_found, verbose){
    # if ( length(l_results) > 0 && return_first_found){ return(l_results)}
    if (length(c_gene_split_value) != 0) {
        if (verbose) {
            message("before predictWithNewValue")
            message(length(c_gene_split_value))}
        mod_obj <- predictWithNewValue(exprs, c_gene_split_value,
            clusters, target, classifier,
            adv_method = adv_method,
            adv_fixed_value = adv_fixed_value, adv_fct = adv_fct,
            verbose = TRUE)
        cell_type <- mod_obj[1]; celltype_score <- mod_obj[2]
        message(cell_type); message(celltype_score)
        if (cell_type != target) {
            if (length(c_gene_split_value) <= max_split_size &&
                (change_type == "any" || cell_type != "NA" )) {
                message("Store gene:")
                message(c_gene_split_value)
                l_results[[paste(c_gene_split_value, collapse = "__")]] <-
                    c(cell_type, celltype_score)
            } else {
                if ( length(c_gene_split_value) > max_split_size){
                    if (verbose) { message("Before predictionDichot:")
                        message(length(c_gene_split_value))}
                    l_results <- .predictionDichotMinChange(exprs,
                        c_gene_split_value, clusters, target, classifier,
                        l_results = l_results, adv_method = adv_method,
                        adv_fixed_value = adv_fixed_value, adv_fct = adv_fct,
                        max_split_size = max_split_size,
                        return_first_found=return_first_found,
                        change_type = change_type, verbose = verbose)
                }
            }
        }
    }
    l_results
}

.predictionDichotMinChange <- function(exprs, genes, clusters, target,
        classifier, l_results = list(), adv_method = "perc99",
        adv_fixed_value = 3, adv_fct = NULL, max_split_size = 1,
        return_first_found = FALSE, change_type = "any",
        verbose = TRUE) {
    if ( return_first_found && length(l_results) > 0) return (l_results)
    c_gene_split <- split(unlist(genes), 1:2)
    message(paste0("genes size: ", length(unlist(genes))))
    l_results <- .dichotMinSplit(l_results,
                    unname(unlist(c_gene_split[1])), exprs,
                    genes, clusters, target, classifier, adv_method,
                    adv_fixed_value, adv_fct,
                    max_split_size, change_type,return_first_found,
                    verbose)
    if ( return_first_found && length(l_results) > 0) return (l_results)
    l_results <- .dichotMinSplit(l_results,
                    unname(unlist(c_gene_split[2])), exprs,
                    genes, clusters, target, classifier,
                    adv_method, adv_fixed_value,
                    adv_fct, max_split_size, change_type,return_first_found,
                    verbose)
    l_results
}

.maxOverListModifs <- function(exprs, clusters, classifier, excl_genes,
                            genes, modifications, max_split_size, verbose){
    df_result <- data.frame(todel = unique(clusters))
    rownames(df_result) <- unique(clusters)
    df_names <- c()
    for (modif_ind in seq_len(length(modifications))) {
        mod1 <- modifications[[modif_ind]][[1]]
        attacks_length <- c()
        for (cell_type in unique(clusters)) {
            if (verbose) {
                message(paste0("Running maxChange attack on ",
                    cell_type, ", with a max_split_size of: ", max_split_size))
                message(paste0("The smaller the max_split_size,",
                        " the more precise the result will be,",
                        " but it will take longer."))
                message("Modification: ",
                    paste(modifications[[modif_ind]], collapse = " "))
            }
            if (length(modifications[[modif_ind]]) == 1) {
                max_change_genes <- advMaxChange(exprs, clusters, cell_type,
                    classifier, adv_method = mod1,
                    max_split_size = max_split_size, excl_genes = excl_genes,
                    genes = genes, verbose = verbose)
            } else {
                mod2 <- modifications[[modif_ind]][[2]]
                max_change_genes <- advMaxChange(exprs, clusters,
                    cell_type, classifier, adv_method = mod1,
                    adv_fixed_value = mod2, adv_fct = mod2,
                    max_split_size = max_split_size, excl_genes = excl_genes,
                    genes = genes, verbose = verbose)
            }
            result_length <- length(max_change_genes)
            attacks_length <- c(attacks_length, result_length)
            if (verbose) {
                message(paste0("At least ", result_length,
                    " genes can be modified with the ",
                    paste(modifications[[modif_ind]], collapse = " "),
                    " method, and the cluster will still",
                    " be classified as ", cell_type))
            }
        }
        df_names <- c(df_names, paste(modifications[[modif_ind]],
            collapse = "_"))
        df_result <- cbind(df_result, attacks_length)
        df_result$todel <- NULL
        colnames(df_result) <- df_names
    }
    df_result
}

.maxOverArgModifs <- function(exprs, clusters, classifier, excl_genes,
                            genes, adv_method, adv_fixed_value,
                            adv_fct, max_split_size, verbose){
    attacks_length <- c()
    for (cell_type in unique(clusters)) {
        if (verbose) {
            message(paste0(
                "Running maxChange attack on ",
                cell_type, ", with a max_split_size of: ", max_split_size
            ))
            message(paste0("The smaller the max_split_size, the more",
                " precise the result will be, but it will take longer."))
        }
        max_change_genes <- advMaxChange(exprs, clusters, cell_type,
            classifier, adv_method = adv_method,
            adv_fixed_value = adv_fixed_value, adv_fct = adv_fct,
            max_split_size = max_split_size, excl_genes = excl_genes,
            genes = genes, verbose = verbose
        )
        result_length <- length(max_change_genes)
        attacks_length[[cell_type]] <- result_length
        if (verbose) {
            message(paste0(
                "At least ", result_length,
                " genes can be modified with the ", adv_method,
                " method, and the cluster will still be classified as ",
                cell_type
            ))
        }
    }
    df_result <- data.frame(at_least_gene_number = unlist(attacks_length))
    rownames(df_result) <- names(attacks_length)
    df_result
}

#' Run an approximation of advMaxChange on all clusters
#' for a given modification. The precision of the
#' approximation can be controlled using the `max_split_size`
#' parameter, with lower values resulting in greater precision
#' but longer processing time. The purpose of this
#' approximation is to determine which clusters are most
#' susceptible to max change adversarial attacks.
#'
#' @param exprs a matrix or dataframe of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param clusters a list of the clusters to which the cells belong
#' @param classifier a classifier in the suitable format
#' @param excl_genes a list of genes to exclude from the analysis
#' @param genes a list of genes in case you want to limit the
#' analysis on a subset of genes
#' @param modifications the list of the modifications to study
#' @param adv_method the name of the method to use
#' @param adv_fixed_value the numeric value to use in case of
#' adv_method=`fixed`
#' @param adv_fct the function to use in case adv_method
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param max_split_size max size of dichotomic slices.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a list of genes/new classification tuples
#' @examples
#' maxChangeOverview(rna_expression, clusters_id,
#' myClassifier, modifications = list(list("perc1"), list("perc99")))
#' @export
maxChangeOverview <- function(exprs, clusters, classifier, excl_genes = c(),
                            genes = c(),
                            modifications = list(list("perc1"), list("perc99")),
                            adv_method = "perc99", adv_fixed_value = 3,
                            adv_fct = NULL, max_split_size = 100,
                            verbose = FALSE) {
    if (length(modifications) == 0) {
        df_result <- .maxOverArgModifs(exprs, clusters, classifier,
            excl_genes, genes, adv_method, adv_fixed_value,
            adv_fct, max_split_size, verbose)
    } else {
        df_result <- .maxOverListModifs(exprs, clusters, classifier, excl_genes,
            genes, modifications, max_split_size, verbose)
    }
    df_result
}


.minOverListModifs <- function(exprs, clusters, classifier, excl_genes,
                            genes, modifications, first_dichot, max_split_size,
                            change_type, verbose){
    df_result <- data.frame(todel = unique(clusters))
    rownames(df_result) <- unique(clusters); df_names <- c()
    for (modif_ind in seq_len(length(modifications))) {
        mod1 <- modifications[[modif_ind]][[1]]
        attacks_length <- c()
        for (cell_type in unique(clusters)) {
            if (verbose) {
                message("Running minChange attack on ", cell_type,
                    ", with a max_split_size of: ", max_split_size)
                message("The smaller the max_split_size, the more precise",
                    " the result will be, but it will take longer.")
                message("Modification: ",
                    paste(modifications[[modif_ind]], collapse = " "))
            }
            if (length(modifications[[modif_ind]]) == 1) {
                min_change_genes <- advMinChange(exprs, clusters, cell_type,
                    classifier, excl_genes = excl_genes, genes = genes,
                    adv_method = mod1, max_split_size = max_split_size,
                    first_dichot = first_dichot, change_type = change_type,
                    verbose = verbose)
            } else {
                mod2 <- modifications[[modif_ind]][[2]]
                min_change_genes <- advMinChange(exprs, clusters, cell_type,
                    classifier, excl_genes = excl_genes, genes = genes,
                    adv_method = mod1, adv_fixed_value = mod2, adv_fct = mod2,
                    max_split_size = max_split_size,
                    first_dichot = first_dichot, change_type = change_type,
                    verbose = verbose)
            }
            result_length <- length(min_change_genes)
            attacks_length <- c(attacks_length, result_length)
            if (verbose) {
                message("An approximation gives about ", result_length,
                    " genes can cause a one gene min change attack on the ",
                    cell_type, " cell type for the modification ",
                    paste(modifications[[modif_ind]], collapse = " "))
            }
        }
        df_names <- c(df_names,
            paste(modifications[[modif_ind]], collapse = "_"))
        df_result <- cbind(df_result, attacks_length)
        df_result$todel <- NULL
        colnames(df_result) <- df_names
    }
    df_result
}


.minOverArgModifs <- function(exprs, clusters, classifier, excl_genes,
                    genes, adv_method, adv_fixed_value, adv_fct,
                    first_dichot, max_split_size, change_type, verbose) {
    attacks_length <- c()
    for (cell_type in unique(clusters)) {
        if (verbose) {
            message(paste0(
                "Running minChange attack on ", cell_type,
                ", with a max_split_size of: ", max_split_size
            ))
            message(
                "The smaller the max_split_size, the more",
                " precise the result will be, but it will take longer."
            )
        }
        min_change_genes <- advMinChange(exprs, clusters, cell_type,
            classifier,
            excl_genes = excl_genes,
            genes = genes, adv_method = adv_method,
            adv_fixed_value = adv_fixed_value,
            adv_fct = adv_fct,
            max_split_size = max_split_size,
            first_dichot = first_dichot, change_type = change_type,
            verbose = verbose
        )
        result_length <- length(min_change_genes)
        attacks_length[[cell_type]] <- result_length
        if (verbose) {
            message(paste0(
                "An approximation gives about ",
                result_length,
                " genes can cause a one gene min change attack on the ",
                cell_type, " cell type"
            ))
        }
    }
    df_result <- data.frame(about_gene_number = unlist(attacks_length))
    rownames(df_result) <- names(attacks_length)
    df_result
}

#' Run an approximation of advMinChange on all clusters
#' for a given modification. The precision of the
#' approximation can be controlled using the `max_split_size`
#' parameter, with lower values resulting in greater precision
#' but longer processing time. The purpose of this
#' approximation is to determine which clusters are most
#' susceptible to min change adversarial attacks.
#'
#' @param exprs a matrix or dataframe of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param clusters a list of the clusters to which the cells belong
#' @param classifier a classifier in the suitable format
#' @param excl_genes a list of genes to exclude from the analysis
#' @param genes a list of genes in case you want to limit the
#' analysis on a subset of genes
#' @param modifications the list of the modifications to study
#' @param adv_method the name of the method to use
#' @param adv_fixed_value the numeric value to use in case of
#' adv_method=`fixed`
#' @param adv_fct the function to use in case adv_method
#' belongs to the following list: `full_row_fct`, `target_row_fct`,
#' `target_matrix_fct`, `full_matrix_fct`
#' @param first_dichot the initial number of slices before
#' the dichotomic search
#' @param max_split_size max size of dichotomic slices.
#' @param change_type `any` consider each misclassification,
#'  `not_na` consider each misclassification but NA.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return a list of genes/new classification tuples
#' @examples
#' minChangeOverview(rna_expression, clusters_id,
#' myClassifier, modifications = list(list("perc1"), list("perc99")))
#' @export
minChangeOverview <- function(exprs, clusters, classifier, excl_genes = c(),
            genes = c(), modifications = list(list("perc1"), list("perc99")),
            adv_method = "perc99", adv_fixed_value = 3, adv_fct = NULL,
            first_dichot = 100, max_split_size = 100, change_type = "any",
            verbose = FALSE) {
    if (length(modifications) == 0) {
        df_result <- .minOverArgModifs(exprs, clusters, classifier, excl_genes,
                            genes, adv_method, adv_fixed_value, adv_fct,
                            first_dichot, max_split_size, change_type, verbose)
    } else {
        df_result <- .minOverListModifs(exprs, clusters, classifier, excl_genes,
                            genes, modifications, first_dichot, max_split_size,
                            change_type, verbose)
    }
    df_result
}

.gridWarning <- function(modifications, genes, iamsure){
    if ((length(modifications) + 1)^length(genes) > 100000 & !iamsure) {
        message("Exit because of too many combinations to test: ",
            (length(modifications) + 1)^length(genes))
        message("This will probably make your computer freeze.")
        message("You should lower the number of genes and/or",
            " of modifications. For example 5 genes and 2 modifications",
            " gives 243 combinations to test")
        message("You can use the iamsure=T option to run the function anyway.")
        return(TRUE)
    } else {
        return (FALSE)
    }
}

#' Grid search of min change adversarial attack. Tries each
#' combination on a cluster, given a list of genes and a list of modifications.
#'
#' @param exprs a matrix or dataframe of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param clusters a list of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format
#' @param genes the list of genes to study
#' @param modifications the list of the modifications to study
#' @param return_first_found set to TRUE to return result when a
#' the first misclassification is found
#' @param verbose logical, set to TRUE to activate verbose mode
#' @param iamsure logical, prevents from expansive calculations
#' when `genes` list is too long, set to `TRUE` to run anyway.
#' @return dataframe results of the classification of all the grid combinations
#' @examples
#' advGridMinChange(rna_expression, clusters_id, "Memory CD4+",
#' myClassifier, genes=c("CD4","CD8","IL2A","IL2B"),
#' modifications = list(list("perc1"), list("perc99")))
#' @export
advGridMinChange <- function(exprs, clusters, target, classifier,
                genes, modifications = list(list("perc1"), list("perc99")),
                return_first_found = FALSE, verbose = FALSE, iamsure = FALSE) {
    function_results <- data.frame(matrix(ncol = length(genes) + 5, nrow = 0))
    colnames(function_results) <- c("prediction", "odd", "genes_modified",
        "type_modified", "iteration", genes)
    if (.gridWarning(modifications, genes, iamsure))
        return(function_results)
    tests_grid <- gtools::permutations(n = length(modifications) + 1,
        r = length(genes), repeats.allowed = TRUE)
    results <- data.frame(matrix(ncol = length(genes) + 2, nrow = 0))
    colnames(results) <- c(genes, "prediction", "odd"); i <- 1
    while (i <= nrow(tests_grid)) {
        message("Running combination: ", i, " on ", nrow(tests_grid))
        exprs_temp <- exprs; row_results <- c()
        for (gene_ind in seq_len(length(genes))) {
            modif_ind <- tests_grid[i, gene_ind]
            if (modif_ind != length(modifications) + 1) {
                mod1 <- modifications[[modif_ind]][[1]]
                if (length(modifications[[modif_ind]]) == 1) {
                    exprs_temp <- advModifications(exprs_temp, genes[gene_ind],
                        clusters, target, adv_method = mod1, verbose = verbose)
                } else { mod2 <- modifications[[modif_ind]][[2]]
                    exprs_temp <- advModifications(exprs_temp, genes[gene_ind],
                        clusters, target, adv_method = mod1, adv_fct = mod2,
                        adv_fixed_value = mod2, verbose = verbose)}
                row_results <- c(row_results,
                    paste(modifications[[modif_ind]], collapse = " "))
            } else { row_results <- c(row_results, "NA")}
        }
        class_results <- classifier(exprs_temp, clusters, target)
        row_results <- c(row_results, class_results[1], class_results[2])
        results <- rbind(results, row_results)
        colnames(results) <- c(genes, "prediction", "odd")
        ifelse(class_results[1] == target, i <- i + 1, i <- nrow(tests_grid)+1)
    }
    results$genes_modified <-
        length(genes) - apply(results, 1, function(x) sum(x == "NA"))
    results$type_modified <- results$prediction != target
    results <- results[, c( "prediction", "odd",
            "genes_modified", "type_modified", genes)]
    function_results <- rbind(results, function_results)
    function_results <-
        function_results[order(as.numeric(function_results$genes_modified)), ]
    function_results <- function_results[order(function_results$type_modified,
            decreasing = TRUE), ]
    function_results
}

.randWalkTryNewVector <- function(exprs, genes, new_walk_params, modifications,
        clusters, target, classifier, best_attack_params, i,
        function_results, change_type, verbose){
    exprs_temp <- exprs
    row_results <- c()
    row_results_int <- c()
    # Try the new params
    for (gene_ind in seq_len(length(genes))) {
        modif_ind <- as.numeric(new_walk_params[gene_ind])
        if (modif_ind != length(modifications) + 1) {
            mod1 <- modifications[[modif_ind]][[1]]
            if (length(modifications[[modif_ind]]) == 1) {
                exprs_temp <- advModifications(exprs_temp,
                    genes[gene_ind], clusters, target,
                    adv_method = mod1, verbose = verbose)
            } else {
                mod2 <- modifications[[modif_ind]][[2]]
                exprs_temp <- advModifications(exprs_temp, genes[gene_ind],
                    clusters, target, adv_method = mod1,
                    adv_fixed_value = mod2, adv_fct = mod2,
                    verbose = verbose)
            }
            row_results <- c(row_results, paste(modifications[[modif_ind]],
                        collapse = " "))
            row_results_int <- c(row_results_int, modif_ind)
        } else {
            row_results <- c(row_results, "NA")
            row_results_int <- c(row_results_int, modif_ind)
        }
    }
    genes_modified <- length(genes) - sum(row_results == "NA")
    previous_genes_modified <- length(genes) -
        sum(best_attack_params == (length(modifications) + 1))
    class_results <- classifier(exprs_temp, clusters, target)
    type_modified <- class_results[1] != target
    if ( change_type == "not_na")
        type_modified <- (class_results[1] != target) &&
            class_results[1] != "NA"
    if (type_modified && genes_modified <= previous_genes_modified) {
        best_attack_params <- row_results_int
        message( "Better attack with only ", genes_modified,
            " genes modified")
    }
    row_results <- c(class_results[1], class_results[2],
        genes_modified, type_modified, (i + 1), row_results)
    function_results <- rbind(row_results, function_results)
    function_results
}

.randWalkNewVector <- function(previous_str_comb, modifications,
    best_attack_params, step_change_ratio, while_max_count, function_results){
    new_walk_params <- c()
    new_walk_step <- -1
    kill_count <- 0
    while (new_walk_step %in% previous_str_comb ||
        length(new_walk_params[
            new_walk_params != (length(modifications) + 1)
        ]) >=
            length(best_attack_params[best_attack_params !=
                (length(modifications) + 1)])) {
        new_walk_params <- c()
        for (param in best_attack_params) {
            if (param != length(modifications) + 1) {
                if (sample(1:round(1 / step_change_ratio), 1) == 1) {
                    new_walk_params <- c( new_walk_params,
                        sample(1:(length(modifications) + 1), 1))
                } else {
                    new_walk_params <- c(new_walk_params, param)
                }
            } else {
                if (sample(1:round(sum(new_walk_params ==
                    (length(modifications) + 1))^(1 / 2) /
                    step_change_ratio), 1) == 1) {
                    new_walk_params <- c( new_walk_params,
                        sample(1:(length(modifications) + 1), 1))
                } else {
                    new_walk_params <- c(new_walk_params, param)
                }
            }
        }
        new_walk_step <- paste(new_walk_params, collapse = " ")
        # It inifite loop kill function and return results
        if (kill_count > while_max_count) {
            function_results <- function_results[order(
                    as.numeric(function_results$genes_modified)), ]
            function_results <- function_results[order(
                function_results$type_modified, decreasing = TRUE), ]
            message("Inifite loop, kill function and return results")
            return(list(FALSE, function_results))
        }
        kill_count <- kill_count + 1
    }
    return(list(TRUE, new_walk_params, new_walk_step))
}

.randWalkBeforeWalk <- function(exprs, genes, modifications,
            clusters, target, classifier, first_batch, verbose){
    init_objs <- .initRandWalk(modifications, genes, first_batch)
    previous_str_comb <- init_objs[[1]]; tests_grid <- init_objs[[2]]
    function_results <- init_objs[[3]]; function_results_int <- init_objs[[4]]
    results <- init_objs[[5]]; results_int <- init_objs[[6]]
    rw_seed <- .randWalkGetSeed(tests_grid, exprs, genes, modifications,
        clusters, target, classifier, results, results_int,
        previous_str_comb, verbose)
    results <- rw_seed[[1]]; results_int <- rw_seed[[2]]
    previous_str_comb <- rw_seed[[3]]

    results$genes_modified <-
        length(genes) - apply(results, 1, function(x) sum(x == "NA"))
    results$type_modified <- results$prediction != target
    results$iteration <- 1
    results <- results[, c("prediction", "odd", "genes_modified",
            "type_modified", "iteration", genes)]
    results_int$genes_modified <- length(genes) - apply( results_int, 1,
                        function(x) sum(x == (length(modifications) + 1)))
    results_int$type_modified <- results_int$prediction != target
    results_int$iteration <- 1
    results_int <- results_int[, c("prediction", "odd", "genes_modified",
        "type_modified", "iteration", genes)]
    function_results <- rbind(results, function_results)
    function_results <-
        function_results[order(function_results$genes_modified), ]
    function_results <- function_results[order(function_results$type_modified,
            decreasing = TRUE), ]
    function_results_int <- rbind(results_int, function_results_int)
    function_results_int <-
        function_results_int[order(function_results_int$genes_modified), ]
    function_results_int <- function_results_int[
        order(function_results_int$type_modified, decreasing = TRUE), ]
    best_attack_params <- function_results_int[1, 6:ncol(function_results_int)]
    list(previous_str_comb, function_results, best_attack_params)
}

.randWalkGetSeed <- function(tests_grid, exprs, genes, modifications, clusters,
    target, classifier, results, results_int, previous_str_comb, verbose){
    i<- 1; new_previous_str_comb <- c()
    while (i <= nrow(tests_grid)) {
        message("Running first batch to determine walk seed: ",
            i, " on ", nrow(tests_grid))
        exprs_temp <- exprs
        new_previous_str_comb <- c(new_previous_str_comb, previous_str_comb[i])
        row_results <- c()
        row_results_int <- c()
        for (gene_ind in seq_len(length(genes))) {
            modif_ind <- tests_grid[i, gene_ind]
            if (modif_ind != length(modifications) + 1) {
                mod1 <- modifications[[modif_ind]][[1]]
                if (length(modifications[[modif_ind]]) == 1) {
                    exprs_temp <- advModifications(exprs_temp,
                        genes[gene_ind], clusters, target,
                        adv_method = mod1, verbose = verbose)
                } else {
                    mod2 <- modifications[[modif_ind]][[2]]
                    exprs_temp <- advModifications(exprs_temp,
                        genes[gene_ind], clusters, target,
                        adv_method = mod1,adv_fixed_value = mod2,
                        adv_fct = mod2, verbose = verbose)
                }
                row_results <- c( row_results,
                    paste(modifications[[modif_ind]], collapse = " "))
                row_results_int <- c(row_results_int, modif_ind)
            } else {
                row_results <- c(row_results, "NA")
                row_results_int <- c(row_results_int, modif_ind)
            }
        }
        class_results <- classifier(exprs_temp, clusters, target)
        row_results <- c(row_results, class_results[1], class_results[2])
        row_results_int <- c(unlist(row_results_int),
                    class_results[1], class_results[2])
        results <- rbind(results, row_results)
        results_int <- rbind(results_int, row_results_int)
        # colnames(results) <- c(genes, "prediction", "odd")
        colnames(results) <- c(genes)
        #colnames(results_int) <- c(genes, "prediction", "odd")
        colnames(results_int) <- c(genes)
        if (class_results[1] == target) {
            i <- i + 1
        } else {
            i <- nrow(tests_grid) + 1
        }
    }
    list(results, results_int, new_previous_str_comb)
}

.initRandWalk <- function(modifications, genes, first_batch){
    previous_str_comb <- c(-1)
    if ((length(modifications) + 1)^length(genes) > first_batch) {
        tests_grid <- data.frame(matrix(ncol = length(genes), nrow = 0))
        while (nrow(tests_grid) < (first_batch - length(modifications))) {
            rand_sample <- sample(1:(length(modifications) + 1),
                length(genes),
                replace = TRUE
            )
            str_rand_sample <- paste(rand_sample, collapse = " ")
            if (!str_rand_sample %in% previous_str_comb) {
                tests_grid <- rbind(tests_grid, rand_sample)
            }
            previous_str_comb <- c(previous_str_comb, str_rand_sample)
        }
        for (i in seq_len(length(modifications))) {
            tests_grid <- rbind(rep(i, length(genes)), tests_grid)
        }
        tests_grid <- as.matrix(tests_grid)
    } else {
        tests_grid <- gtools::permutations(n = length(modifications) + 1,
            r = length(genes), repeats.allowed = TRUE)
    }
    function_results <- data.frame(matrix(ncol = length(genes) + 5, nrow = 0))
    colnames(function_results) <- c("prediction", "odd", "genes_modified",
        "type_modified", "iteration", genes)
    function_results_int <- data.frame(
        matrix(ncol = length(genes) + 5, nrow = 0)
    )
    colnames(function_results_int) <- c( "prediction", "odd",
        "genes_modified", "type_modified", "iteration", genes)
    results <- data.frame(matrix(ncol = length(genes) + 2, nrow = 0))
    colnames(results) <- c(genes, "prediction", "odd")
    results_int <- data.frame(matrix(ncol = length(genes) + 2, nrow = 0))
    colnames(results_int) <- c(genes, "prediction", "odd")
    return (list(previous_str_comb, tests_grid, function_results,
        function_results_int, results, results_int))
}

#' Random walk search of min change adversarial attack.
#' Step 1 is to find a seed by trying random combinations of
#' genes and modifications on a cluster until the classification is altered.
#' Step 2 is to perform a random walk search to reduce the number of genes
#' needed to change the classification."
#' 
#' @param exprs a matrix or dataframe of numeric RNA expression,
#' cells are rows and genes are columns.
#' @param clusters a list of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format
#' @param genes the list of genes to study
#' @param modifications the list of the modifications to study
#' @param first_batch the maximum number of try in step 1
#' @param walk_length the maximum number of try in step 2
#' @param step_change_ratio ratio of parameters change in new walk step
#' @param while_max_count the maximum number of try when looking
#' for new combination of parameters
#' @param change_type `any` consider each misclassification,
#'  `not_na` consider each misclassification but NA.
#' @param verbose logical, set to TRUE to activate verbose mode
#' @return dataframe results of the classification of all the grid combinations
#' @examples
#' advRandWalkMinChange(rna_expression, clusters_id, "Memory CD4+",
#' myClassifier, genes=c("CD4","CD8","IL2A","IL2B"),
#' modifications = list(list("perc1"), list("perc99")))
#'
#' # Stop at first attack discovery, whitout going into the walk
#' # parameter search.
#' advRandWalkMinChange(rna_expression, clusters_id, "Memory CD4+",
#' myClassifier, genes=c("CD4","CD8","IL2A","IL2B"),
#' modifications = list(list("perc1"), list("perc99")), walk_length=0)
#' @export
advRandWalkMinChange <- function(exprs, clusters, target, classifier, genes,
        modifications = list(list("perc1"), list("perc99")), first_batch = 100,
        walk_length = 100, step_change_ratio = 0.2, while_max_count = 10000,
        change_type = "any", verbose = FALSE) {
    before_walk <- .randWalkBeforeWalk(exprs, genes, modifications,
        clusters, target, classifier, first_batch, verbose)
    previous_str_comb <- before_walk[[1]]
    function_results <- before_walk[[2]]
    best_attack_params <- before_walk[[3]]
    if (function_results[1, "type_modified"] == FALSE) {
        message("No modified type, try with a higher first_batch argument")
        return(function_results)
    }
    for (i in 1:walk_length) {
        message("Walk step ", i, " on ", walk_length)
        if (verbose) {
            message("Current best attack length: ",
                length(best_attack_params), " genes")
        }
        # Create a new vector of parameters to test
        rw_new_vector <- .randWalkNewVector(previous_str_comb, modifications,
            best_attack_params, step_change_ratio, while_max_count,
            function_results)
        # If infinite loop kill function and return result
        if ( !rw_new_vector[[1]] ) {return(rw_new_vector[[2]])}
        new_walk_params <- rw_new_vector[[2]]
        new_walk_step <- rw_new_vector[[3]]

        previous_str_comb <- c(previous_str_comb, new_walk_step)
        function_results <- .randWalkTryNewVector(exprs, genes,
            new_walk_params, modifications, clusters, target,
            classifier, best_attack_params, i, function_results,
            change_type, verbose)
    }
    function_results <- function_results[order(
            as.numeric(function_results$genes_modified)), ]
    function_results <- function_results[order(
            function_results$type_modified, decreasing = TRUE), ]
    function_results
}

#' @export
todel <- function(a, b) {
    a + b
}
