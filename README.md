
# Introduction
**adverSCarial is a package for generating and evaluating vulnerability to adversarial attacks on single-cell RNA sequencing classifiers.**
Single cell analysis are on the way to be used in clinical routine. As a critical use, algorithms must address the challenge of ensuring reliability, one concern being the susceptibility to adversarial attacks.


Zhang, J., Wang, W., Huang, J., Wang, X., & Zeng, Y. (2020). How far is single-cell sequencing from clinical application?. Clinical and translational medicine, 10(3), e117. https://doi.org/10.1002/ctm2.117

de Hond, A.A.H., Leeuwenberg, A.M., Hooft, L. et al. Guidelines and quality criteria for artificial intelligence-based prediction models in healthcare: a scoping review. npj Digit. Med. 5, 2 (2022). https://doi.org/10.1038/s41746-021-00549-7


The package is designed to generate and analyze the vulnerability of scRNA-seq classifiers to adversarial attacks. The package is versatile and provides a format for integrating
any type of classifier. It offers functions for studying and generating two types of attacks, single gene attack and max change attack.
The single gene attack involves making a small modification to the input to alter the classification. This is an issue when the change of classification is made on a gene that is not known to be biologicaly involved in the change of this cell type.
The max change attack involves making a large modification to the input without changing its classification. This is an issue when a significant percentage of genes, sometimes as high as 99%, can be modified with the classifier still predicting the same cell type.

# Jupyter Notebook examples
- [Train, format and CGD attack on a scAnnotatR classifier](https://github.com/GhislainFievet/adverSCarial_article/blob/master/04_train_and_attack_scAnnotatR.ipynb)
- [Train, format and max-change attack on a random forest classifier](https://github.com/GhislainFievet/adverSCarial_article/blob/master/05_train_format_random_forest_classifier_scRF.ipynb)
- [Train, format and single-gene attack on a multi layer perceptron classifier](https://github.com/GhislainFievet/adverSCarial_article/blob/master/06_train_format_classifier_scMLP.ipynb)


# Installation
```{r installation, eval=FALSE}
## Install BiocManager is necessary
if (!require("BiocManager")) {
    install.packages("BiocManager")
}
BiocManager::install(version='devel')
BiocManager::install('adverSCarial')
```

# Generate an adversarial attack

There are three types of adversarial attacks: `single-gene attack` which modifies the expression of a gene in order to alter the classification, `max-change attack` which introduce the largest possible perturbations to the input while still keeping the same classification, and `Cluster Gradient Descent (CGD)` which modifies the input matrix gene by gene by following an approximated gradient until the cluster has switch its classification.

Load libraries

```{r load libraries, warning = FALSE, message=FALSE}
library(adverSCarial)
library(LoomExperiment)
```

Load a pbmc raw data

```{r min change attack, message=FALSE, warning = FALSE}
pbmcPath <- system.file("extdata", "pbmc_short.loom", package="adverSCarial")
lfile <- import(pbmcPath, type="SingleCellLoomExperiment")
```

```{r extract matrix data}
matPbmc <- counts(lfile)
```

```{r visualize1, message=FALSE, warning = FALSE}
matPbmc[1:5, 1:5]
```


and its cell type classification based on [Seurat documentation](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
```{r load annnots, message=FALSE, warning = FALSE}
cellTypes <- rowData(lfile)$cell_type
```

```{r visualize2, message=FALSE, warning = FALSE}
head(cellTypes)
```

The package contains a marker based classifier working on the pbmc3k dataset: `MClassifier`.
We verify it is classifying properly according to the `cellTypes`.

```{r first classif, message=FALSE}
for ( cell_type in unique(cellTypes)){
    resClassif <- MClassifier(matPbmc, cellTypes, cell_type)
    print(paste0("Expected: ", cell_type, ", predicted: ", resClassif[1]))
}
```

We want to run attacks on the "DC" cluster but without modifying the genes used by human to make manual classification. We can find this list on [Seurat documentation](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)

```{r markers}
# Known markers for each cell type
markers <- c("IL7R", "CCR7", "CD14", "LYZ", "S100A4", "MS4A1", "CD8A", "FCGR3A", "MS4A7",
              "GNLY", "NKG7", "FCER1A", "CST3", "PPBP")
```

## Single-gene attack

Modify the value of a gene inside a cell cluster in order to change its classification.

The function `advSingleGene` runs a dichotomic search of one gene attacks, given the cluster to attack and a type of modification. We suggest two main modifications: `perc1` replacing the value of the gene by its first percentile, and `perc99` replacing the value of the gene by its 99th percentile. The `adv_fct` argument allows users to choose custom modifications of the gene.

The `excl_genes` argument allows users to exclude certain genes from being modified. Here we exclude the genes usually used as markers of specific cell types, as defined previously.

Computation times can be lengthy, we use the `return_first_found` argument to return the result as soon as it is found.

We can specify the `change_type = "not_na"` argument to indicate that we want the attack to misclassify the cluster as an existing cell type, rather than as NA.

```{r search min change attack, message=FALSE, warning = FALSE}
genesMinChange <- advSingleGene(matPbmc, cellTypes, "DC",
                        MClassifier, exclGenes = markers, advMethod = "perc99",
                        returnFirstFound = TRUE, changeType = "not_na",
                        firstDichot = 10)
```
```{r search min change attack bis, warning = FALSE}
genesMinChange
```

The function found that modifying the single gene `r names(genesMinChange@values)[1]` with the `perc99` modification on the cluster leads to a new classification, `r genesMinChange@values[[1]][1]`.

Let's run this attack and verify if it is successful.

First we modify the `matPbmc` matrix on the target cluster.

```{r run min change attack, message=FALSE, warning = FALSE}
matAdver <- advModifications(matPbmc, names(genesMinChange@values)[1], cellTypes, "DC")
```

Then classify the "DC" cluster with `MClassifier`.
```{r verify the min change attack, warning = FALSE, message = FALSE}
resClassif <- MClassifier(matAdver, cellTypes, "DC")
```
```{r verify the min change attack bis, warning = FALSE}
resClassif
```

## Max-change attack

Modify the maximum number of genes inside a cell cluster without altering its classification.

The function `advMaxChange` runs a dichotomic search of gene subsets, given the cluster to attack and a type of modification, such that the classification does not change.

The `maxSplitSize` argument is the maximum size of dichotomic slices. Set to 1 to have better results, but it will take longer to compute.

```{r search max change attack, warning = FALSE, message=FALSE}
genesMaxChange <- advMaxChange(matPbmc, cellTypes, "Memory CD4 T", MClassifier,
                    exclGenes = markers, advMethod = "perc99")
```
```{r search max change attack bis, warning = FALSE}
length(genesMaxChange@values)
```


The function found `r length(genesMaxChange@values)` genes that you can modify with the `perc99` modification, and the cluster is still classified as `Memory CD4 T`.

Let's run this attack and verify if it is successful.

First we modify the `matPbmc` matrix on the target cluster, on the genes previously determined.

```{r run max change attack, message=FALSE, warning=FALSE}
matMaxAdver <- advModifications(matPbmc, genesMaxChange@values, cellTypes, "Memory CD4 T")
```

Then we verify that classification is still `Memory CD4 T`.
```{r verify max change attack, warning = FALSE, message=FALSE}
resClassif <- MClassifier(matMaxAdver, cellTypes, "Memory CD4 T")
```
```{r verify max change attack bis, warning = FALSE}
resClassif
```


## CGD attack

Apply the Cluster Gradient Descent (CGD) attack on the "Memory CD4 T" cells cluster with parameters alpha=epsilon=1

```{r CGD1, warning = FALSE}
resCGD <- advCGD(as.data.frame(matPbmc), cellTypes,
                "Memory CD4 T", MClassifier, alpha=1, epsilon=1,
                genes=colnames(matPbmc)[ncol(matPbmc):1])
```

The algorithm test each gene and stops when the cluster prediction has switched
to another cell type prediction. Here it stops after having modified the CCR5 gene which make the classifier switch its prediction to Naive CD4 T.

```{r CGD2, warning = FALSE}
tail(resCGD$byGeneSummary)
```
```{r CGD2_bis, warning = FALSE}
resCGD$oneRowSummary
```
You can retrieve the modified matrix.

```{r CGD3, warning = FALSE}
modifiedMatrix <- resCGD$expr
```

And visualize the list of modified genes.

```{r CGD4, warning = FALSE}
resCGD$modGenes
```

Check the new classification of the modified matrix.

```{r CGD5, warning = FALSE}
MClassifier(modifiedMatrix, cellTypes, "Memory CD4 T")$prediction
```

```{r session info}
sessionInfo()
```
