[![DOI](https://zenodo.org/badge/614626969.svg)](https://zenodo.org/badge/latestdoi/614626969)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/guidedPLS)](
https://cran.r-project.org/package=guidedPLS)
[![Downloads](https://cranlogs.r-pkg.org/badges/guidedPLS)](https://CRAN.R-project.org/package=guidedPLS)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/guidedPLS?color=orange)](https://CRAN.R-project.org/package=guidedPLS)
[![:name status badge](https://rikenbit.r-universe.dev/badges/:name)](https://rikenbit.r-universe.dev)
[![:registry status badge](https://rikenbit.r-universe.dev/badges/:registry)](https://rikenbit.r-universe.dev)
[![:total status badge](https://rikenbit.r-universe.dev/badges/:total)](https://rikenbit.r-universe.dev)
[![guidedPLS status badge](https://rikenbit.r-universe.dev/badges/guidedPLS)](https://rikenbit.r-universe.dev)
![GitHub Actions](https://github.com/rikenbit/guidedPLS/actions/workflows/build_test_push.yml/badge.svg)

# guidedPLS
R package for supervised dimensional reduction by guided partial least squares

## Installation

======
~~~~
git clone https://github.com/rikenbit/guidedPLS/
R CMD INSTALL guidedPLS
~~~~
or type the code below in the R console window
~~~~
library(devtools)
devtools::install_github("rikenbit/guidedPLS")
~~~~

## Sparse matrix support (v1.2.0)

`guidedPLS()` accepts `Matrix::dgCMatrix` (sparse matrices) for `X1` and `X2`.
This is particularly useful for single-cell data (scRNA-seq, scATAC-seq) where feature matrices are large and sparse.

Internally, the row-centered cross-product is computed without materializing the dense centered matrix, using the identity:

```
X_centered %*% Y = X %*% Y - mu * (1^T Y)
```

where `mu` is the vector of row means and `1^T Y = colSums(Y)`.
Only the small `k x p` result matrix is dense, keeping memory usage proportional to the number of non-zero entries.

## References
======
- **Sparse Partial Least Squares (sPLS)/Sparse Partial Least Squares Discriminant Analysis (sPLS-DA)**
  - Lê Cao, et al. A Sparse PLS for Variable Selection when Integrating Omics Data. Statistical Applications in Genetics and Molecular Biology, 7(1), 2008
- **Guided Principal Component Analysis (guided-PCA)**
  - Reese S E, et al. A new statistic for identifying batch effects in high-throughput genomic data that uses guided principal component analysis. Bioinformatics, 29(22), 2877-2883, 2013
- **Multiset Canonical Correlations Analysis**
  - Nielsen A A. Multiset Canonical Correlations Analysis and Multispectral, Truly Multitemporal Remote Sensing Data. IEEE Transactions on Image Processing, 11(3), 293-305, 2002

## Contributing

If you have suggestions for how `guidedPLS` could be improved, or want to report a bug, open an issue! We'd love all and any contributions.

For more, check out the [Contributing Guide](CONTRIBUTING.md).

## Authors
- Koki Tsuyuzaki
