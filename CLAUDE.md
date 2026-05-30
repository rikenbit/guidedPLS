# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
# Install package locally
R CMD INSTALL .

# Run all tests via devtools
Rscript -e 'devtools::test()'

# Run a single test file
Rscript -e 'devtools::test(filter="sparse")'

# Full R CMD check
R CMD build . && R CMD check guidedPLS_*.tar.gz

# Load package in development mode (no install needed)
Rscript -e 'devtools::load_all(".")'
```

## Architecture

An R package for supervised dimensionality reduction via guided partial least squares.

### Core algorithm (guidedPLS)

The main function `guidedPLS(X1, X2, Y1, Y2, k)` computes metadata-guided cross-covariance between two datasets:

1. **Input**: X1 (n1 x p1), X2 (n2 x p2) data matrices; Y1 (n1 x k), Y2 (n2 x k) metadata matrices
2. **Init** (`.initguidedPLS`): Row-standardize X and Y, compute cross-product `YX = scale(t(Y) %*% X)`, then `M = t(YX1) %*% YX2`
3. **SVD**: Decompose M to get loadings; project X onto loadings for scores

Two modes: default SVD-based and SUMCOR-based CCA (`sumcor=TRUE`, requires geigen).

### Sparse matrix support (v1.2.0)

X1/X2 accept `Matrix::dgCMatrix`. The sparse path avoids materializing the dense row-centered matrix using the identity `X_c %*% Y = X %*% Y - mu * (1^T Y)`. Helper functions in `guidedPLS-internal.R`:

- `.crossprod_row_standardized(X, Y_std)` — computes `t(Y_std) %*% row_standardize(X)` without densifying X
- `.tcrossprod_row_centered(X, Y_c)` — computes `t(row_center(X)) %*% Y_c` without densifying X
- `.row_mean_sd(X)` — row means/sds via `rowMeans(X)` and `rowSums(X^2)` (both sparse-safe)

SUMCOR mode still densifies X because it needs p x p covariance matrices.

### Function organization

Each algorithm follows the pattern: `.check*()` for input validation, `.init*()` for preprocessing, then the main function body. Internal functions use dot-prefix naming.

### Test structure

Tests are in `tests/testthat/` and registered in `tests/testthat.R`. Test files must be added to `testthat.R` manually. Tests validate output dimensions and numerical equivalence (e.g., sparse vs dense paths at tolerance 1e-8). SUMCOR tests are conditional on geigen being installed.
