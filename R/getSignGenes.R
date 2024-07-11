#' The function getSignGenes orders the genes by maximizing the
#' significance of the gene to differentiate the clusters and ensures
#' that they represent at most the variations across all possible pairs
#' of clusters.
#' @details The function getSignGenes orders the genes by maximizing the
#' significance of the gene to differentiate the clusters and ensures
#' that they represent at most the variations across all possible pairs
#' of clusters.
#' 
#' The getDistantCouples function is used to generate all possible pairs
#' of clusters and order them in a way that the following pairs are as
#' distant as possible.
#' cell_types = c("B", "CD4 T", "NK", "CD8 T", "DC")
#' getDistantCouples(cell_types)
#'      'B___CD4 T', 'NK___CD8 T', 'B___DC', 'CD4 T___NK',
#'      'CD8 T___DC', 'B___NK', 'CD4 T___CD8 T', 'NK___DC',
#'      'B___CD8 T', 'CD4 T___DC'
#' 
#' @param expr A matrix of gene expression data. Rows are cells and columns
#' are genes.
#' @param clusters a character vector of the clusters to which the cells belong
#' @param method the statistical test to use. Either "wilcox" for the Wilcoxon
#' rank sum test or "ttest" for the t-test. Default is "wilcox".
#' @param verbose logical, set to TRUE to activate verbose mode
#' @examples
#' rna_expression <- data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'     CD8B=c(2,2,3,3))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","T cell","T cell","T cell")
#' 
#' getSignGenes(rna_expression, clusters_id, method="wilcox", verbose=TRUE)
#' 
#' @export 
getSignGenes <- function(expr, clusters, method="wilcox", verbose=FALSE){
    combinations <- getDistantCouples(unique(clusters))
    pvalList <- list()
    
    for(strComb in combinations){
        a_a <- unlist(strsplit(strComb, "___"))
        clust1 <- a_a[1]
        clust2 <- a_a[2]
        if (verbose) {message("Cluster ",clust1," vs ",clust2)}
        pvals <- apply(t(expr), 1, function(x){
                c1 <- x[clusters == clust1]
                c1 <- c1[!is.na(c1)]
                c2 <- x[clusters == clust2]
                c2 <- c2[!is.na(c2)]
                if ( length(c1) == 0 || length(c2) == 0){
                    return(1)
                }
                if (length(unique(c1))==1){
                    c1[1] = c1[1] + 0.01
                }
                if (length(unique(c2))==1){
                    c2[1] = c2[1] + 0.01
                }
                if (method=="wilcox"){
                    return(wilcox.test(c1, c2)$p.value)
                }
                if (method=="ttest"){
                    return(t.test(c1, c2)$p.value)
                }    
        })
        dfPvals <- data.frame(gene=colnames(expr), pval=unname(pvals))
        rownames(dfPvals) <- dfPvals$gene
        dfPvals <- dfPvals[order(dfPvals$pval),]
        pvalList[[strComb]] <- dfPvals
    }
    
    dfResults <- as.data.frame(matrix(nrow=0, ncol=2))
    colnames(dfResults) <- c("gene", "pval")
    
    for (i in 1:ncol(expr)){
        for(strCouple in names(pvalList)){
            strGene <- pvalList[[strCouple]][i,"gene"]
            if (!strGene %in% dfResults$gene){
                dfResults <- rbind(dfResults, pvalList[[strCouple]][i,])
            }
        }
    }

    results <- list(results=dfResults, pvalList=pvalList)
    return(results)
}

getDistantCouples <- function(cell_types, verbose=FALSE){
    combinations <- combn(cell_types, 2, simplify = FALSE)
    strCombs <- c()
    for (i in combinations){
        strCombs <- c(strCombs, paste0(i[1], "___", i[2]))
    }
    twoPrevs <- c()
    selCombs <- c()
    intInd <- 1
    killLoop <- FALSE
    intKill <- 0
    while(length(selCombs)<length(combinations)){
        if (intKill == 0 && killLoop){
            selCombs <- c(selCombs, setdiff(strCombs, selCombs))
        }
        intKill <- intKill + 1
        killLoop <- TRUE
        c1 <- combinations[[intInd]][1]
        c2 <- combinations[[intInd]][2]
        myKey <- paste0(c1,"___",c2)

        if ( !c1 %in% twoPrevs && !c2 %in% twoPrevs && !myKey %in% selCombs){
            if (verbose){message(myKey)}
            intKill <- 0
            killLoop <- FALSE
            selCombs <- c(selCombs, myKey)
            twoPrevs <- c(c1, c2)
        }
        intInd <- intInd + 1
        if (intKill > length(combinations) +1){
            intKill <- 0
        }
        if (intInd > length(combinations)){
            intInd <- 1
        }
    }
    return(selCombs)
}