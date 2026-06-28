# Test sparse/dense numerical equivalence for guidedPLS

library(Matrix)

# --- Toy data ---
set.seed(42)
n1 <- 30
n2 <- 25
p1 <- 50
p2 <- 40
k <- 4

# Dense matrices (samples x features)
X1_dense <- matrix(rpois(n1 * p1, lambda = 0.5), nrow = n1, ncol = p1)
X2_dense <- matrix(rpois(n2 * p2, lambda = 0.3), nrow = n2, ncol = p2)

# Convert to sparse
X1_sparse <- Matrix(X1_dense, sparse = TRUE)
X2_sparse <- Matrix(X2_dense, sparse = TRUE)

# Metadata matrices (one-hot, samples x groups)
Y1 <- dummyMatrix(rep(1:k, length.out = n1), center = FALSE)
Y2 <- dummyMatrix(rep(1:k, length.out = n2), center = FALSE)

# --- Test 1: .crossprod_row_standardized equivalence ---
test_that("sparse crossprod_row_standardized matches dense path", {
    Y1_std <- t(scale(t(Y1), center = TRUE, scale = TRUE))

    # Dense reference
    X1_std_dense <- t(scale(t(X1_dense), center = TRUE, scale = TRUE))
    ref <- t(Y1_std) %*% X1_std_dense

    # Sparse path
    result <- guidedPLS:::.crossprod_row_standardized(X1_sparse, Y1_std)

    expect_equal(as.vector(result), as.vector(ref), tolerance = 1e-10)
})

# --- Test 2: .tcrossprod_row_centered equivalence ---
test_that("sparse tcrossprod_row_centered matches dense path", {
    Y1_c <- t(scale(t(Y1), center = TRUE, scale = FALSE))

    # Dense reference
    X1_c_dense <- t(scale(t(X1_dense), center = TRUE, scale = FALSE))
    ref <- t(X1_c_dense) %*% Y1_c

    # Sparse path
    result <- guidedPLS:::.tcrossprod_row_centered(X1_sparse, Y1_c)

    expect_equal(as.vector(result), as.vector(ref), tolerance = 1e-10)
})

# --- Test 3: Full guidedPLS (SVD mode) equivalence ---
test_that("guidedPLS with sparse X matches dense X results", {
    out_dense <- guidedPLS(X1_dense, X2_dense, Y1, Y2, k = 3)
    out_sparse <- guidedPLS(X1_sparse, X2_sparse, Y1, Y2, k = 3)

    # Singular values must match
    expect_equal(out_dense$res$d, out_sparse$res$d, tolerance = 1e-8)

    # Loadings may differ by sign flip; compare absolute values
    expect_equal(abs(out_dense$loadingYX1), abs(out_sparse$loadingYX1), tolerance = 1e-8)
    expect_equal(abs(out_dense$loadingYX2), abs(out_sparse$loadingYX2), tolerance = 1e-8)

    # Scores may differ by sign flip
    expect_equal(abs(out_dense$scoreX1), abs(out_sparse$scoreX1), tolerance = 1e-8)
    expect_equal(abs(out_dense$scoreX2), abs(out_sparse$scoreX2), tolerance = 1e-8)
    expect_equal(abs(out_dense$scoreYX1), abs(out_sparse$scoreYX1), tolerance = 1e-8)
    expect_equal(abs(out_dense$scoreYX2), abs(out_sparse$scoreYX2), tolerance = 1e-8)
})

# --- Test 4: guidedPLS with irlba (fullrank=FALSE) ---
test_that("guidedPLS with sparse X and irlba matches dense", {
    out_dense <- guidedPLS(X1_dense, X2_dense, Y1, Y2, k = 2, fullrank = FALSE)
    out_sparse <- guidedPLS(X1_sparse, X2_sparse, Y1, Y2, k = 2, fullrank = FALSE)

    expect_equal(out_dense$res$d, out_sparse$res$d, tolerance = 1e-8)
    expect_equal(abs(out_dense$loadingYX1), abs(out_sparse$loadingYX1), tolerance = 1e-8)
    expect_equal(abs(out_dense$loadingYX2), abs(out_sparse$loadingYX2), tolerance = 1e-8)
})

# --- Test 5: guidedPLS with cortest ---
test_that("guidedPLS with sparse X and cortest matches dense", {
    out_dense <- guidedPLS(X1_dense, X2_dense, Y1, Y2, k = 2, cortest = TRUE)
    out_sparse <- guidedPLS(X1_sparse, X2_sparse, Y1, Y2, k = 2, cortest = TRUE)

    expect_equal(out_dense$res$d, out_sparse$res$d, tolerance = 1e-8)
    expect_equal(abs(out_dense$corYX1), abs(out_sparse$corYX1), tolerance = 1e-8)
    expect_equal(abs(out_dense$corYX2), abs(out_sparse$corYX2), tolerance = 1e-8)
    expect_equal(out_dense$pvalYX1, out_sparse$pvalYX1, tolerance = 1e-8)
    expect_equal(out_dense$pvalYX2, out_sparse$pvalYX2, tolerance = 1e-8)
})

# --- Test 6: SUMCOR mode sparse equivalence ---
if(requireNamespace("iTensor", quietly = TRUE)){
    test_that("guidedPLS SUMCOR with sparse X matches dense", {
        out_dense <- guidedPLS(X1_dense, X2_dense, Y1, Y2, k = 2,
                               sumcor = TRUE, lambda = 1e-4)
        out_sparse <- guidedPLS(X1_sparse, X2_sparse, Y1, Y2, k = 2,
                                sumcor = TRUE, lambda = 1e-4)

        expect_equal(out_dense$res$d, out_sparse$res$d, tolerance = 1e-6)
        expect_equal(abs(out_dense$loadingYX1), abs(out_sparse$loadingYX1), tolerance = 1e-6)
        expect_equal(abs(out_dense$loadingYX2), abs(out_sparse$loadingYX2), tolerance = 1e-6)
    })
}

# --- Test 7: Mixed sparse/dense (X1 sparse, X2 dense) ---
test_that("guidedPLS works with mixed sparse/dense inputs", {
    out_ref <- guidedPLS(X1_dense, X2_dense, Y1, Y2, k = 2)
    out_mixed <- guidedPLS(X1_sparse, X2_dense, Y1, Y2, k = 2)

    expect_equal(out_ref$res$d, out_mixed$res$d, tolerance = 1e-8)
    expect_equal(abs(out_ref$loadingYX1), abs(out_mixed$loadingYX1), tolerance = 1e-8)
    expect_equal(abs(out_ref$loadingYX2), abs(out_mixed$loadingYX2), tolerance = 1e-8)
})

# --- Test 8: Dimension checks preserved ---
test_that("output dimensions are correct for sparse input", {
    out <- guidedPLS(X1_sparse, X2_sparse, Y1, Y2, k = 3)

    expect_equal(dim(out$loadingYX1), c(p1, 3))
    expect_equal(dim(out$loadingYX2), c(p2, 3))
    expect_equal(dim(out$scoreX1), c(n1, 3))
    expect_equal(dim(out$scoreX2), c(n2, 3))
    expect_equal(dim(out$scoreYX1), c(ncol(Y1), 3))
    expect_equal(dim(out$scoreYX2), c(ncol(Y2), 3))
    expect_equal(length(out$res$d), 3)
})

# --- Test 9: Zero-density edge case (very sparse) ---
test_that("guidedPLS handles very sparse matrices", {
    X1_very_sparse <- Matrix(0, nrow = n1, ncol = p1, sparse = TRUE)
    X1_very_sparse[1, 1] <- 1
    X1_very_sparse[2, 5] <- 2
    X1_very_sparse[10, 20] <- 3

    X2_very_sparse <- Matrix(0, nrow = n2, ncol = p2, sparse = TRUE)
    X2_very_sparse[1, 1] <- 1
    X2_very_sparse[3, 10] <- 2

    # Should run without error
    out <- guidedPLS(X1_very_sparse, X2_very_sparse, Y1, Y2, k = 2)
    expect_equal(length(out$res$d), 2)
})
