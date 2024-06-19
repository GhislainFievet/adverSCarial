#' @export
advTraining <- function(expr, classifier, genes=NULL,
                    coef=0.1, add=0, slot=NULL, geneNumber=100,
                    modGeneNumber=10000000, verbose=FALSE){
    if (class(expr) == "Seurat"){
        cell_number <- nrow(expr@assays$RNA@cells)
    } else {
        cell_number <- nrow(expr)
    }
    predictions <- classifier(expr, rep("all", cell_number), "all")
    cellPredictions <- predictions$cellTypes

    if (is.null(genes) && class(expr) != "Seurat"){
        genes <- colnames(expr)[sample(1:ncol(expr), ncol(expr))]
    }
    if ( is.null(genes) && class(expr) == "Seurat"){
        genes <- expr@assays$RNA@features[sample(1:length(expr@assays$RNA@features), length(expr@assays$RNA@features))]
   }

    modExpr <- expr
    tempExpr <- expr
    geneCount <- 1
    modGeneCount <- 0
    geneInDS <- TRUE
    modifiedGenes <- c()
    while(geneCount <= geneNumber && modGeneCount <= modGeneNumber){
        
        myGene <- genes[geneCount]
        if(verbose){
            message(myGene, " ", geneCount)
        }
    
        if (class(expr) == "Seurat"){
            dfScaled <- as.data.frame(modExpr@assays$RNA@layers[[slot]])
            colnames(dfScaled) <- rownames(modExpr@assays$RNA@cells)
            rownames(dfScaled) <- rownames(modExpr@assays$RNA@features)
            dfScaled <- as.data.frame(t(dfScaled))
            # modify only cells out of target cells
            if ( !myGene %in% colnames(dfScaled)){
                myGene <- sub( "\\.","-",myGene)
                if (!myGene %in% colnames(dfScaled)){
                    geneInDS <- FALSE
                }
            }
            if (geneInDS){
                dfScaled[, myGene][dfScaled[, myGene]>0] <- dfScaled[, myGene][dfScaled[, myGene]>0]*(1+coef) + add
                dfScaled[, myGene][dfScaled[, myGene]<= 0] <- dfScaled[, myGene][dfScaled[, myGene] <= 0]*(1-coef) + add
                modExpr@assays$RNA@layers[[slot]] <- t(dfScaled)
            }
        } else {
            # modExpr[, myGene] <- modExpr[, myGene]*(1+coef) + add    
            modExpr[, myGene][modExpr[, myGene]>0] <- modExpr[, myGene][modExpr[, myGene]>0]*(1+coef) + add
            modExpr[, myGene][modExpr[, myGene]<= 0] <- modExpr[, myGene][modExpr[, myGene] <= 0]*(1-coef) + add
        }
        
        newPredictions <- classifier(modExpr, rep("all", cell_number), "all")
        newPredictions$cellTypes
        newPredictions$typePredictions

        dfnp <- t(newPredictions$typePredictions[,colnames(predictions$typePredictions)])
        dfp <- t(predictions$typePredictions)
        allCols <- unique(c(colnames(dfnp), colnames(dfp)))
        dfnp[,setdiff(allCols, colnames(dfnp))] <- 0
        dfp[,setdiff(allCols, colnames(dfp))] <- 0
        dfnp <- dfnp[, order(colnames(dfnp))]
        dfp <- dfp[,colnames(dfnp)]
        
        # Make sure colnames are the same, same length
        colnames(dfnp) <- paste0("new_",colnames(dfnp))
        colnames(dfp) <-paste0("or_",colnames(dfp))
        df <- cbind(dfnp, dfp)
        scoreDelts <- t(apply(df, 1, function(x){
            newVals <- x[grep("^new_", names(x), value=TRUE)]
            orVals <- x[grep("^or_", names(x), value=TRUE)]
            newVals - orVals
        }))
        
        # Which cells benefit from modif
        scoreDelts <- as.data.frame(scoreDelts)
        scoreDelts$cell_id <- rownames(scoreDelts)
        scoreDelts$prevCellPredictions <- cellPredictions

        # Which cells benefit from modif
        cbfm <- unlist(apply(scoreDelts, 1, function(x){
            y <- x[order(x, decreasing = FALSE)]
            cellOrTarget <- names(y)[1]
            cellNewTarget <- names(y)[2]
            if (as.numeric(unname(x[cellNewTarget])) > as.numeric(unname(x[cellOrTarget]))) {
                return(x["cell_id"])
            }
        }))
        # Which cells benefit from opposit modif
        cbfom <- unlist(apply(scoreDelts, 1, function(x){
            y <- x[order(x, decreasing = FALSE)]
            cellOrTarget <- names(y)[1]
            cellNewTarget <- names(y)[2]
            if (as.numeric(unname(x[cellNewTarget])) < as.numeric(unname(x[cellOrTarget]))) {
                return(x["cell_id"])
            }
        }))

        modExpr <- tempExpr

        if ((length(cbfm) + length(cbfom)) > 0){
            modifiedGenes <- c(modifiedGenes, myGene)
            if (class(expr) == "Seurat"){
                dfScaled <- as.data.frame(modExpr@assays$RNA@layers[[slot]])
                colnames(dfScaled) <- rownames(modExpr@assays$RNA@cells)
                rownames(dfScaled) <- rownames(modExpr@assays$RNA@features)
                dfScaled <- as.data.frame(t(dfScaled))
                dfScaled[cbfm, myGene][dfScaled[cbfm, myGene]>0] <- dfScaled[cbfm, myGene][dfScaled[cbfm, myGene]>0]*(1+coef) + add
                dfScaled[cbfm, myGene][dfScaled[cbfm, myGene]<= 0] <- dfScaled[cbfm, myGene][dfScaled[cbfm, myGene] <= 0]*(1-coef) + add
#                 dfScaled[cbfm, myGene] <- dfScaled[cbfm, myGene]*(1+coef) + add
                dfScaled[cbfom, myGene][dfScaled[cbfom, myGene]>0] <- dfScaled[cbfom, myGene][dfScaled[cbfom, myGene]>0]*(1-coef) - add
                dfScaled[cbfom, myGene][dfScaled[cbfom, myGene]<= 0] <- dfScaled[cbfom, myGene][dfScaled[cbfom, myGene] <= 0]*(1+coef) - add
#                 dfScaled[cbfom, myGene] <- dfScaled[cbfom, myGene]*(1-coef) - add
                modExpr@assays$RNA@layers[[slot]] <- t(dfScaled)
                modGeneCount <- modGeneCount + 1
            } else {
                modExpr[cbfm, myGene] <- modExpr[cbfm, myGene]*(1+coef) + add
                modExpr[cbfom, myGene] <- modExpr[cbfom, myGene]*(1-coef) - add
            }
            tempExpr <- modExpr
        }
        geneCount <- geneCount + 1
    }
    return(modExpr)
}