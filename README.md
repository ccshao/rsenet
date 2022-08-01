# Introduction
`rsenet` aims to resolve the spatial location of single cells (profiled in scRNA-seq experiments) with previous known information (bins). Based on the same genes detected in both query and reference datasets, `rsenet` perform regularized regression (Elastic Net/ Lasso) with the expression from single cell experiments as responsible variables, and values from bins as predictor variables. The coefficients of bins considered as the assignment probabilities.Additionally, an adaptive technique that uses ridge regression to find the penalties is employed to address the problem of strongly correlated bins. Consensus assignment are obtained via running the analysis multiple times.

Generally rsenet sloves the following equations, while c1, c2, ..., are cells to be aligned, and B1, ..B1 are known bins. The coefficients b1, ..., are foced to be non-negative. 
$$E_{c1g1} = B0 + b_{1} * B_{1} + b_{2} * B_{2} + b_{3} * B_{3} +  ...$$
$$E_{c1g2} = B0 + b_{1} * B_{1} + b_{2} * B_{2} + b_{3} * B_{3} +  ...$$
$$...$$

# Analysis
First install the pacakge form github.
```r
remotes::install_github("ccshao/rsenet")
```

It is highly recommeded to install parallelization backends.
```r
install.packages("doRNG")
install.packages("doFuture")
```

A subset of scRNA-seq and reference bins data from fly embyro [1]  are included in the package. Please refer to the publication for more information, including steps of normalization. The inference could be run easily as the following codes.
```r
library(doFuture)
library(progressr)

registerDoFuture()
plan(multicore, workers = 4)

#- Reference and scRNA-seq data, respectively.
data(example_bin)
data(example_expr)

#- Returns bins by cells matrix.
res <- enet(example_bin, example_expr)

#- With progress bar.
with_progress(res <- enet(example_bin, example_expr))
```

# Comparison
Nikos Karaiskos et al. [1] developed an method, `DistMap`, to assign cells to bins. As shown bellow, `DistMap` generate a continuous mcc score (correlation between cells and bins) which have little difference among bins, thus could not give an accurate estimation of cells locations.

![mcc_dist](./inst/image/highest_mcc_score.png)
![mcc_diff](./inst/image/diff_in_perc.png)

`rsent` is able to sparsely allocate cells to bins because of regularized regressions. In addtions, cells that are unlikely to locate in the same spatial regions are identified. The results from `rsent` are comparable with `DistMap`, with about 60% of bins predicted are among the top bins within 1sd.

![enet_prob](./inst/image/enet_pred_dis.png)
![enet_comp](./inst/image/mcc_sore_1sd.png)

# Reference
[1] Karaiskos N, Wahle P, Alles J, Boltengagen A, Ayoub S, Kipar C, et al. The Drosophila embryo at single-cell transcriptome resolution. Science. 2017;358:194–9.
