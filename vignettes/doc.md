---
title: "adverSCarial"
shorttitle: "Short title for headers"
author: Ghislain FIEVET <ghislain.fievet@gmail.com>
package: adverSCarial
abstract: >
  Document summary
output:
    BiocStyle::html_document:
        toc: true
        toc_depth: 2
vignette: >
  %\VignetteIndexEntry{adverSCarial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Introduction
**adverSCarial is a package for generating and evaluating vulnerability to adversarial attacks on single-cell RNA sequencing classifiers.**

# Installation

```r
## Install BiocManager is necessary
if (!require("BiocManager")) {
    install.packages("BiocManager")
}
BiocManager::install('adverSCarial')
```

# Generate an adversarial attack

There are two types of adversarial attacks: `min change attack` which modifies slightly the input in order to alter the classification, and `max change attack` which introduce the largest possible perturbations to the input while still keeping the same classification.

Load libraries



Load a pbmc SingleCellExperiment and its cell type classification.



The package contains a random forest based classifier working on the pbmc3k dataset: `RFClassifier`. We verify it is classifying properly.


```r
RFClassifier(mat_pbmc, cell_types$cell_type, "DC")
```

```
## Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('randomForest.formula', 'randomForest')"
```

We want to run attacks on the "DC" cluster but without modifying the genes used by human to make manual classification. We can find this list on [Seurat documentation](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)


```r
# Known markers for each cell type
markers = c("IL7R", "CCR7", "CD14", "LYZ", "S100A4", "MS4A1", "CD8A", "FCGR3A", "MS4A7",
              "GNLY", "NKG7", "FCER1A", "CST3", "PPBP")
```

## Min change attack

Modify the value of a gene inside a cell cluster in order to change its classification.

The function `advMinChange` runs a dichotomic search of one gene attacks, given the cluster to attack and a type of modification. We suggest two main modifications: `perc1` replacing the value of the gene by its first percentile, and `perc99` replacing the value of the gene by its 99th percentile. The `adv_fct` argument allows users to choose custom modifications of the gene.

The `excl_genes` argument allows users to exclude certain genes from being modified. Here we exclude the genes usually used as markers of specific cell types, as defined previously.

Computation times can be lengthy, we use the `return_first_found` argument to return the result as soon as it is found.

We can specify the `change_type = "not_na"` argument to indicate that we want the attack to misclassify the cluster as an existing cell type, rather than as NA.


```r
adv_min_change = advMinChange(mat_pbmc, cell_types, "DC",
                        RFClassifier, excl_genes = markers, adv_method = "perc99",
                        return_first_found = T, change_type = "not_na",
                        first_dichot = 10)
```

```
## Warning in split.default(genes_index, 1:first_dichot): data length is not a
## multiple of split variable
```

```
## Warning in split.default(unlist(genes), 1:2): data length is not a multiple of
## split variable
```

```
## Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('randomForest.formula', 'randomForest')"
```

```r
print(adv_min_change)
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'print': object 'adv_min_change' not found
```

The function found that modifying the single gene `names(adv_min_change)[1]` with the `perc99` modification on the cluster leads to a new classification, `adv_min_change[[1]][1]`.

Let's run this attack and verify if it is successful.

First we modify the `mat_pbmc` matrix on the target cluster.


```
## Error in unique(genes): object 'adv_min_change' not found
```

Then classify the "DC" cluster with `RFClassifier`.

```r
RFClassifier(mat_adver, cell_types, "DC")
```

```
## Error in is.data.frame(x): object 'mat_adver' not found
```


## Max change attack

Modify the maximum number of genes inside a cell cluster without altering its classification.

The function `advMaxChange` runs a dichotomic search of gene subsets, given the cluster to attack and a type of modification, such that the classification does not change.

The `max_split_size` argument is the maximum size of dichotomic slices. Set to 1 to have better results, but it will take longer to compute.


```r
adv_max_change = advMaxChange(mat_pbmc, cell_types, "DC", RFClassifier,
                    excl_genes = markers, max_split_size = 1000, adv_method = "perc99")
```

```
## new version 3
```

```
## genes size: 22028
```

```
## current gene results length: 0
```

```
## Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('randomForest.formula', 'randomForest')"
```

```r
print(dim(adv_max_change))
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'print': object 'adv_max_change' not found
```

The function found `length(adv_max_change)` genes that you can modify with the `perc99` modification, and the cluster is still classified as `DC`.

Let's run this attack and verify if it is successful.

First we modify the `mat_pbmc` matrix on the target cluster, on the genes previously determined.


```
## Error in unique(genes): object 'adv_max_change' not found
```

Then we verify that classification is still `DC`.

```r
RFClassifier(mat_adver, cell_types, "DC")
```

```
## Error in is.data.frame(x): object 'mat_adver' not found
```


# Prepare a classifier with `CHETAH`

# Advanced attacks with `advGridMinChange` and `advRandWalkMinChange`

# Run vulnerability analysis with `minChangeOverview` and `maxChangeOverview`
