#' Implementation of the IDG4C algorithm.
#' 
#' @details this function is an implementation of the IDG4C algorithm
#' which permits to generate adversarial attacks on a classifier on a given
#' cluster. The attack is done by modifying the expression of the genes by 
#' gradient descent after having inferred if by the finite difference method.
#' Two parameters alpha and epsilon are used to control the step size of the 
#' modifications.
#' 
#' @param expr a matrix, a data.fram or a Seurat object
#' @param clusters a character vector of the clusters to which the cells belong
#' @param target the name of the cluster to modify
#' @param classifier a classifier in the suitable format.
#' A classifier function should be formated as follow:
#' classifier = function(expr, clusters, target){
#'      # Making the classification
#'      c("cell type", score)
#' }
#' `score` should be numeric between 0 and 1, 1 being the highest confidance
#' into the cell type classification.
#' The matrix `expr` contains RNA expression values, the vector `clusters`
#' consists of the cluster IDs for each cell in `expr`, and `target` is the
#' ID of the cluster for which we want to have a classification.
#' The function returns a vector with the classification result, and a score.
#' @param genes the character vector of genes to study
#' @param exclNewTargets the character vector of cell types to exclude as new target
#' @param newTarget the name of the new cell type cell type target, this will be the
#' new prediction cell type after the attack
#' @param alpha the alpha parameter of the IDG4C algorithm
#' @param epsilon the epsilon parameter of the IDG4C algorithm
#' @param slot the slot to modify in case of Seurat object
#' @param stopAtSwitch logical, set to TRUE to stop the attack when the new target
#' set to FALSE to continue the attack until all genes are tested
#' @param verbose logical, set to TRUE to activate verbose mode
#' @param returnCTpredMat logical, set to TRUE to return the cell type predictions matrix
#' @return a list containing the modified expression matrix, the list of modified genes,
#' the summary of the attack by gene, the summary of the attack,
#' the new cell types predictions and the original cell types predictions
#' @examples
#' MyClassifier <- function(expr, clusters, target) {
#'    typePredictions <- as.data.frame(matrix(nrow=nrow(expr), ncol=length(unique(clusters))))
#'    colnames(typePredictions) <- unique(clusters)
#'    typePredictions[unique(clusters)[1]] <- c(1,0,0,0)
#'    typePredictions[unique(clusters)[2]] <- c(0,1,1,1)
#'    rownames(typePredictions) <- 1:4
#'
#'    list(prediction="T cell", odd=1,
#'                      typePredictions=t(typePredictions),
#'                      cellTypes=c("B cell","T cell","T cell","T cell"))
#' }
#' rna_expression <- data.frame(CD4=c(0,0,0,0), CD8A=c(1,1,1,1),
#'     CD8B=c(2,2,3,3))
#' genes <- c("CD4", "CD8A")
#' clusters_id <- c("B cell","T cell","T cell","T cell")
#' 
#' advCGD(rna_expression, clusters_id, "T cell",
#'   MyClassifier, genes=genes, verbose=TRUE)
#' 
#' 
#' @export
advCGD <- function(expr, clusters, target, classifier, genes=NULL,
                     exclNewTargets=NULL, newTarget=NULL,
                     alpha=0.1, epsilon=0, slot=NULL,
			         stopAtSwitch=TRUE, verbose=FALSE,
                     returnCTpredMat=FALSE){
    # orExpr <- expr
    predictions <- classifier(expr, clusters, target)
    orTarget <- predictions$prediction
    if (verbose){
      message(orTarget)
      message(dim(t(predictions$typePredictions)))
      message(length(clusters))
      message(sum(clusters == target))
    }
    if (class(expr) == "Seurat"){
        cluster_cell_ids <- rownames(expr@assays$RNA@cells)[clusters == target]
    } else {
        cluster_cell_ids <- names(clusters)[clusters == target]
    }
    
    dfSumUp <- as.data.frame(matrix(nrow=0, ncol=7))
    
    
    cellPredictions <- apply(t(predictions$typePredictions[, clusters == target]), 1, function(x){
            names(x[x == max(x)])[1]
        })
    oriCellPredictions <- cellPredictions
    #display(table(unlist(cellPredictions)))
#     display(names(cellPredictions[cellPredictions=='Myeloid Dendritic cells']))

    # Define the cell type after the attack
    if (is.null(newTarget)){
        orTable <- table(cellPredictions)
        orTable <- orTable[order(orTable, decreasing=T)]
        if (length(names(orTable)) == 1 ){
            # only on cell type predicted, so take the second cell type with the higher odd
            clustTypeOdds <- apply(predictions$typePredictions[, clusters == target], 1, function(x){
                mean(x)
            })
            clustTypeOdds <- clustTypeOdds[order(clustTypeOdds, decreasing = TRUE)]
            if (names(clustTypeOdds)[1] == orTarget){
                newTarget <- names(clustTypeOdds)[2]
            } else {
                newTarget <- names(clustTypeOdds)[1]
            }
        } else {
            # take second most present predicted cell type
            newTarget <- setdiff(names(orTable), c(exclNewTargets, orTarget))[1]
        }
    }
    rightCells <- names(cellPredictions[cellPredictions==newTarget])
    orCells <- names(cellPredictions[cellPredictions==orTarget])
    nOrTarg <- length(orCells)
    c4summaryStart <- c(nOrTarg/sum(clusters==target),
                    length(rightCells)/sum(clusters==target))
    prevMeanOrCell <- mean(unlist(predictions$typePredictions[orTarget,orCells]))
    if ( length(rightCells) == 0){
        prevMeanRightCell <- 0
    } else {
        prevMeanRightCell <- mean(unlist(predictions$typePredictions[newTarget,rightCells]))
    }
    deltMeanRightCell <- 0
    deltMeanOrCell <- 0
    
    if(verbose){message("New cluster target: ",newTarget)}

    
    if (is.null(genes) && class(expr) != "Seurat"){
        genes <- colnames(expr)[sample(1:ncol(expr), ncol(expr))]
    }
    if ( is.null(genes) && class(expr) == "Seurat"){
        genes <- expr@assays$RNA@features[sample(1:length(expr@assays$RNA@features), length(expr@assays$RNA@features))]
   }

    modExpr <- expr
    tempExpr <- expr
    modifiedGenes <- c()
    geneCount <- 1
    geneInDS <- TRUE
    dfSumUp <- rbind(dfSumUp, c(0, "#", 0, nOrTarg/sum(clusters==target),
                              length(rightCells)/sum(clusters==target),
                               prevMeanOrCell, prevMeanRightCell,
                               mean(unlist(predictions$typePredictions[orTarget,cluster_cell_ids])),
                              mean(unlist(predictions$typePredictions[newTarget,cluster_cell_ids]))))
    for(myGene in genes){
        if(verbose){
            message(myGene, " ", geneCount)
            message("Number of original annot ", orTarget," : ", nOrTarg)
            message(" mean ",  prevMeanOrCell, " delt ", deltMeanOrCell)
            message("Number of ", newTarget, " : ", length(rightCells))
            message(" mean ", prevMeanRightCell, " delt ", deltMeanRightCell)
            
        }
        
        if (class(expr) == "Seurat"){
            dfScaled <- as.data.frame(modExpr@assays$RNA@layers[[slot]])
            colnames(dfScaled) <- rownames(modExpr@assays$RNA@cells)
            rownames(dfScaled) <- rownames(modExpr@assays$RNA@features)
            dfScaled <- as.data.frame(t(dfScaled))
            # modify only cells out of target cells
            cells2modify <- rownames(dfScaled)[clusters == target]
            cells2modify <- setdiff(cells2modify, rightCells)
            if ( !myGene %in% colnames(dfScaled)){
                myGene <- sub( "\\.","-",myGene)
                if (!myGene %in% colnames(dfScaled)){
                    geneInDS <- FALSE
                }
            }
            if (geneInDS){
                dfScaled[cells2modify, myGene][dfScaled[cells2modify, myGene]>0] <- dfScaled[cells2modify, myGene][dfScaled[cells2modify, myGene]>0]*(1+alpha) + epsilon
                dfScaled[cells2modify, myGene][dfScaled[cells2modify, myGene]<= 0] <- dfScaled[cells2modify, myGene][dfScaled[cells2modify, myGene] <= 0]*(1-alpha) + epsilon
                modExpr@assays$RNA@layers[[slot]] <- t(dfScaled)
            }
        } else {
            # modify only cells out of target cells
            cells2modify <- rownames(modExpr)[clusters == target]
            cells2modify <- setdiff(cells2modify, rightCells)
            # modExpr[cells2modify, myGene] <- modExpr[cells2modify, myGene]*(1+alpha) + epsilon    
            modExpr[cells2modify, myGene][modExpr[cells2modify, myGene]>0] <- modExpr[cells2modify, myGene][modExpr[cells2modify, myGene]>0]*(1+alpha) + epsilon
            modExpr[cells2modify, myGene][modExpr[cells2modify, myGene]<= 0] <- modExpr[cells2modify, myGene][modExpr[cells2modify, myGene] <= 0]*(1-alpha) + epsilon
        }
        
        newPredictions <- classifier(modExpr, clusters, target)
        newCellPredictions <- apply(t(newPredictions$typePredictions[, clusters == target]), 1, function(x){
            names(x[x == max(x)])[1]
        })
        
        dfnp <- as.data.frame(t(newPredictions$typePredictions[,colnames(predictions$typePredictions[, cells2modify])]))
        if (nrow(dfnp) == 0){
            message("break")
            break
        }
        dfp <- t(predictions$typePredictions[, cells2modify])
        allCols <- unique(c(colnames(dfnp), colnames(dfp)))
        dfnp[,setdiff(allCols, colnames(dfnp))] <- 0
        dfp[,setdiff(allCols, colnames(dfp))] <- 0
        dfnp <- dfnp[, order(colnames(dfnp))]
        dfp <- dfp[,colnames(dfnp)]
        
#         display(head(dfnp))
#         display(head(dfp))
#         display(all(rownames(dfnp) == rownames(dfp)))
#         display(all(colnames(dfnp) == colnames(dfp)))
        
        # Make sure colnames are the same, same length
        colnames(dfnp) <- paste0("new_",colnames(dfnp))
        colnames(dfp) <-paste0("or_",colnames(dfp))
        df <- cbind(dfnp, dfp)
        scoreDelts <- t(apply(df, 1, function(x){
            newVals <- x[grep("^new_", names(x), value=TRUE)]
            orVals <- x[grep("^or_", names(x), value=TRUE)]
            newVals - orVals
        }))
        

#         display(head(df2disp[,c("newCellPredictions", "new_Myeloid Dendritic cells", colnames(scoreDelts))]))

#         cellNewPredictions <- apply(t(newPredictions$typePredictions[, clusters == target]), 1, function(x){
#             names(x[x == max(x)])[1]
#         })
        # If attack succeeded
            # try attack on tempExpr
                # return result with summary of what has changed
        # Which cells benefit from modif
        scoreDelts <- as.data.frame(scoreDelts)
        scoreDelts$cell_id <- rownames(scoreDelts)
        scoreDelts$prevCellPredictions <- cellPredictions[cells2modify]

        # Which cells benefit from modif
        cbfm <- unlist(apply(scoreDelts[scoreDelts$prevCellPredictions != newTarget,], 1, function(x){
            if (as.numeric(unname(x[paste0("new_",newTarget)])) > as.numeric(unname(x[paste0("new_",orTarget)]))) {
                return(x["cell_id"])
            }
        }))
        # Which cells benefit from opposit modif
        cbfom <- unlist(apply(scoreDelts[scoreDelts$prevCellPredictions != newTarget,], 1, function(x){
            if (as.numeric(unname(x[paste0("new_",newTarget)])) < as.numeric(unname(x[paste0("new_",orTarget)]))) {
                return(x["cell_id"])
            }
        }))
        
        if (verbose){
            message("Number of modified cells ", (length(cbfm) + length(cbfom)))
        }
        modExpr <- tempExpr
        
        if ((length(cbfm) + length(cbfom)) > 0){

            modifiedGenes <- c(modifiedGenes, myGene)
            if (class(expr) == "Seurat"){
                dfScaled <- as.data.frame(modExpr@assays$RNA@layers[[slot]])
                colnames(dfScaled) <- rownames(modExpr@assays$RNA@cells)
                rownames(dfScaled) <- rownames(modExpr@assays$RNA@features)
                dfScaled <- as.data.frame(t(dfScaled))
                dfScaled[cbfm, myGene][dfScaled[cbfm, myGene]>0] <- dfScaled[cbfm, myGene][dfScaled[cbfm, myGene]>0]*(1+alpha) + epsilon
                dfScaled[cbfm, myGene][dfScaled[cbfm, myGene]<= 0] <- dfScaled[cbfm, myGene][dfScaled[cbfm, myGene] <= 0]*(1-alpha) + epsilon
#                 dfScaled[cbfm, myGene] <- dfScaled[cbfm, myGene]*(1+alpha) + epsilon
                dfScaled[cbfom, myGene][dfScaled[cbfom, myGene]>0] <- dfScaled[cbfom, myGene][dfScaled[cbfom, myGene]>0]*(1-alpha) - epsilon
                dfScaled[cbfom, myGene][dfScaled[cbfom, myGene]<= 0] <- dfScaled[cbfom, myGene][dfScaled[cbfom, myGene] <= 0]*(1+alpha) - epsilon
#                 dfScaled[cbfom, myGene] <- dfScaled[cbfom, myGene]*(1-alpha) - epsilon
                modExpr@assays$RNA@layers[[slot]] <- t(dfScaled)
            } else {
                modExpr[cbfm, myGene] <- modExpr[cbfm, myGene]*(1+alpha) + epsilon
                modExpr[cbfom, myGene] <- modExpr[cbfom, myGene]*(1-alpha) - epsilon
            }
            
            tempExpr <- modExpr

            tempPredictions <- classifier(modExpr, clusters, target)
            predictions <- tempPredictions
            cellPredictions <- apply(t(tempPredictions$typePredictions[, clusters == target]),
                                     1, function(x){
                names(x[x == max(x)])[1]
            })

            orTable <- table(cellPredictions)
            orTable <- orTable[order(orTable, decreasing=T)]
            if (length(names(orTable)) == 1 ){
                # only on cell type predicted, so take the second cell type with the higher odd
                clustTypeOdds <- apply(predictions$typePredictions[, clusters == target], 1, function(x){
                    mean(x)
                })
                clustTypeOdds <- clustTypeOdds[order(clustTypeOdds, decreasing = TRUE)]
                if (names(clustTypeOdds)[1] == orTarget){
                    tempNewTarget <- names(clustTypeOdds)[2]
                } else {
                    tempNewTarget <- names(clustTypeOdds)[1]
                }
            } else {
                # take second most present predicted cell type
                tempNewTarget <- setdiff(names(orTable), c(exclNewTargets, orTarget))[1]
            }
            if(tempNewTarget != newTarget && verbose){
                message("New cluster target: ",tempNewTarget)
            }
            newTarget <- tempNewTarget
            rightCells <- names(cellPredictions[cellPredictions==newTarget])
            if (length(rightCells) == 0){
                tempPrevMeanRightCell <- 0
            } else {
                tempPrevMeanRightCell <- mean(unlist(tempPredictions$typePredictions[newTarget,rightCells]))
            }
            deltMeanRightCell <- tempPrevMeanRightCell - prevMeanRightCell
            prevMeanRightCell <- tempPrevMeanRightCell
            
            orCells <- names(cellPredictions[cellPredictions==orTarget])
            nOrTarg <- length(orCells)
            tempPrevMeanOrCell <- mean(unlist(tempPredictions$typePredictions[orTarget,orCells]))
            #display("if if if")
            #display(tempPrevMeanOrCell)
            #display(prevMeanOrCell)
            deltMeanOrCell <- tempPrevMeanOrCell - prevMeanOrCell
            prevMeanOrCell <- tempPrevMeanOrCell
            
            if (length(rightCells)>nOrTarg && stopAtSwitch){
                #display(1)
                #display(table(cellPredictions)[order(table(cellPredictions), decreasing = T)])
                
                dfSumUp <- rbind(dfSumUp, c(geneCount, myGene,
                    (length(cbfm) + length(cbfom)),
                    nOrTarg/sum(clusters==target),
                    length(rightCells)/sum(clusters==target),
                    prevMeanOrCell, prevMeanRightCell,
                    mean(unlist(tempPredictions$typePredictions[orTarget,cluster_cell_ids])),
                    mean(unlist(tempPredictions$typePredictions[newTarget,cluster_cell_ids]))))
                colnames(dfSumUp) <- c("order","gene", "nbModifiedCells",
                   "oriCell_pourc", "targCell_pourc",
                   "oriCell_odd", "targCell_odd",
                   "wholeClust_oriOdd", "wholeClust_targOdd")
                cSummary <- c(orTarget,sum(clusters==target), newTarget, alpha, epsilon, 0,
                    length(modifiedGenes), geneCount,
                    c4summaryStart,
                    nOrTarg/sum(clusters==target),
                    length(rightCells)/sum(clusters==target),
                    dfSumUp[1, "wholeClust_oriOdd"],
                    dfSumUp[1, "wholeClust_targOdd"],
                    dfSumUp[nrow(dfSumUp), "wholeClust_oriOdd"],
                    dfSumUp[nrow(dfSumUp), "wholeClust_targOdd"]
                    )
                return(list(expr=modExpr, modGenes=modifiedGenes,
                    byGeneSummary=dfSumUp, oneRowSummary=cSummary,
                    modCellTypes=cellPredictions,
                    oriCellTypes=oriCellPredictions))
            }
            dfSumUp <- rbind(dfSumUp, c(geneCount, myGene,
                (length(cbfm) + length(cbfom)),
                nOrTarg/sum(clusters==target),
                length(rightCells)/sum(clusters==target),
                prevMeanOrCell, prevMeanRightCell,
                mean(unlist(tempPredictions$typePredictions[orTarget,cluster_cell_ids])),
                mean(unlist(tempPredictions$typePredictions[newTarget,cluster_cell_ids]))))
        } else {
            deltMeanOrCell <- 0
            deltMeanRightCell <- 0
            colnames(dfSumUp) <- c("order", "gene", "nbModifiedCells",
                   "oriCell_pourc", "targCell_pourc",
                   "oriCell_odd", "targCell_odd",
                   "wholeClust_oriOdd", "wholeClust_targOdd")
            dfSumUp <- rbind(dfSumUp, c(geneCount, myGene, (length(cbfm) + length(cbfom)),
                    nOrTarg/sum(clusters==target),
                    length(rightCells)/sum(clusters==target),
                    prevMeanOrCell, prevMeanRightCell,
                    dfSumUp[nrow(dfSumUp), "wholeClust_oriOdd"],
                    dfSumUp[nrow(dfSumUp), "wholeClust_targOdd"])
                    )
        }

        geneCount <- geneCount + 1
    }
    
    #display(2)
    tempPredictions <- classifier(modExpr, clusters, target)
    cellPredictions <- apply(t(tempPredictions$typePredictions[, clusters == target]),
                                 1, function(x){
            names(x[x == max(x)])[1]
        })
    #display(table(cellPredictions)[order(table(cellPredictions), decreasing = T)])
    colnames(dfSumUp) <- c("order","gene", "nbModifiedCells",
                       "oriCell_pourc", "targCell_pourc",
                       "oriCell_odd", "targCell_odd",
                   "wholeClust_oriOdd", "wholeClust_targOdd")
    
    cSummary <- c(orTarget, sum(clusters==target), newTarget, alpha, epsilon, 0,
                length(modifiedGenes), length(genes),
                 c4summaryStart, nOrTarg/sum(clusters==target),
                 length(rightCells)/sum(clusters==target),
                 # mean odd of ori cell type in the cluster at the beginning and at the end
                 dfSumUp[1, "wholeClust_oriOdd"],
                 dfSumUp[1, "wholeClust_targOdd"],
                 dfSumUp[nrow(dfSumUp), "wholeClust_oriOdd"],
                 dfSumUp[nrow(dfSumUp), "wholeClust_targOdd"]
                 )
    
    if (!returnCTpredMat){
        return(list(expr=modExpr, modGenes=modifiedGenes,
                byGeneSummary=dfSumUp, oneRowSummary=cSummary,
                modCellTypes=cellPredictions,
                oriCellTypes=oriCellPredictions))
    } else {
        return(list(expr=modExpr, modGenes=modifiedGenes,
                byGeneSummary=dfSumUp, oneRowSummary=cSummary,
                modCellTypes=cellPredictions,
                oriCellTypes=oriCellPredictions,
                modCellTypesMat=tempPredictions$typePredictions,
                oriCellTypesMat=predictions$typePredictions))
    }
}
