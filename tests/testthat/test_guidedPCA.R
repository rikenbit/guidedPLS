context("guidedPCA")

test_that("guidedPCA works with basic inputs", {
  # Create test data
  set.seed(123)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- data.frame(
    group = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    score = rnorm(n),
    binary = sample(c(TRUE, FALSE), n, replace = TRUE)
  )
  
  # Test basic functionality
  result <- guidedPCA(X, Y, k = 3)
  
  expect_equal(ncol(result$loadingX), 3)
  expect_equal(ncol(result$loadingY), 3)
  expect_equal(ncol(result$scoreX), 3)
  expect_equal(ncol(result$scoreY), 3)
  expect_equal(length(result$d), 3)
  expect_equal(nrow(result$scoreX), n)
  expect_equal(nrow(result$loadingX), p)
  
  # Check variance explained sums to <= 1
  expect_true(sum(result$variance_explained) <= 1)
  
  # Check contributions sum to 1
  expect_true(all(abs(colSums(result$contrib_features) - 1) < 1e-10))
  expect_true(all(abs(colSums(result$contrib_groups) - 1) < 1e-10))
})

test_that("guidedPCA handles different data types correctly", {
  set.seed(456)
  n <- 50
  p <- 30
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  # Test with various data types
  Y <- data.frame(
    factor_var = factor(sample(letters[1:4], n, replace = TRUE)),
    numeric_var = runif(n),
    integer_var = sample(1:10, n, replace = TRUE),
    logical_var = sample(c(TRUE, FALSE), n, replace = TRUE),
    character_var = sample(c("yes", "no"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  result <- guidedPCA(X, Y, k = 2)
  
  expect_equal(ncol(result$scoreX), 2)
  expect_true(ncol(result$Y_dummy) >= ncol(Y))  # Dummy encoding expands factors
  expect_equal(length(unique(result$Y_groups)), ncol(Y))
})

test_that("guidedPCA handles missing values", {
  set.seed(789)
  n <- 60
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  # Introduce NAs
  Y <- data.frame(
    group = factor(sample(c("A", "B", NA), n, replace = TRUE)),
    score = rnorm(n)
  )
  Y$score[sample(n, 5)] <- NA
  
  expect_error(result <- guidedPCA(X, Y, k = 2), NA)  # Should not error
  expect_equal(ncol(result$scoreX), 2)
})

test_that("guidedPCA deflation mode works", {
  set.seed(111)
  n <- 80
  p <- 40
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  
  result_deflation <- guidedPCA(X, Y, k = 3, deflation = TRUE)
  result_no_deflation <- guidedPCA(X, Y, k = 3, deflation = FALSE)
  
  # Both should return same dimensions
  expect_equal(dim(result_deflation$loadingX), dim(result_no_deflation$loadingX))
  expect_equal(dim(result_deflation$scoreX), dim(result_no_deflation$scoreX))
  
  # But values should be different
  expect_false(all(abs(result_deflation$loadingX - result_no_deflation$loadingX) < 1e-10))
})

test_that("guidedPCA centering and scaling options work", {
  set.seed(222)
  n <- 70
  p <- 25
  X <- matrix(rnorm(n * p, mean = 10, sd = 5), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  
  result_default <- guidedPCA(X, Y, k = 2)
  result_no_center <- guidedPCA(X, Y, k = 2, center_X = FALSE)
  result_no_scale <- guidedPCA(X, Y, k = 2, scale_X = FALSE)
  result_no_normalize <- guidedPCA(X, Y, k = 2, normalize_Y = FALSE)
  
  # All should return valid results
  expect_equal(ncol(result_default$scoreX), 2)
  expect_equal(ncol(result_no_center$scoreX), 2)
  expect_equal(ncol(result_no_scale$scoreX), 2)
  expect_equal(ncol(result_no_normalize$scoreX), 2)
  
  # Results should differ based on preprocessing
  expect_false(all(abs(result_default$d - result_no_center$d) < 1e-10))
  expect_false(all(abs(result_default$d - result_no_scale$d) < 1e-10))
})

test_that("guidedPCA print and summary methods work", {
  set.seed(333)
  n <- 50
  p <- 30
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("Gene", 1:p)
  Y <- data.frame(
    celltype = factor(sample(c("TypeA", "TypeB", "TypeC"), n, replace = TRUE)),
    treatment = factor(sample(c("Control", "Treated"), n, replace = TRUE))
  )
  
  result <- guidedPCA(X, Y, k = 2)
  
  # These should not error
  expect_output(print(result))
  expect_output(summary(result))
})

test_that("guidedPCA handles single-level factors correctly", {
  set.seed(444)
  n <- 40
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  # Include a factor with only one level
  Y <- data.frame(
    constant = factor(rep("A", n)),
    variable = factor(sample(c("B", "C"), n, replace = TRUE)),
    numeric = rnorm(n)
  )
  
  result <- guidedPCA(X, Y, k = 2)
  
  expect_equal(ncol(result$scoreX), 2)
  expect_true("constant" %in% result$Y_groups)
})

test_that("guidedPCA works with matrix Y input", {
  set.seed(555)
  n <- 60
  p <- 35
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  colnames(Y) <- paste0("Meta", 1:4)
  
  result <- guidedPCA(X, Y, k = 3)
  
  expect_equal(ncol(result$scoreX), 3)
  expect_equal(ncol(result$Y_dummy), 4)  # Numeric columns remain as-is
  expect_equal(unique(result$Y_groups), colnames(Y))
})

test_that("guidedPCA error handling", {
  set.seed(666)
  n <- 30
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  
  # Test invalid k
  expect_error(guidedPCA(X, Y, k = 0))
  expect_error(guidedPCA(X, Y, k = 100))  # k too large
  
  # Test dimension mismatch
  Y_wrong <- matrix(rnorm((n + 5) * 3), nrow = n + 5, ncol = 3)
  expect_error(guidedPCA(X, Y_wrong, k = 2))
  
  # Test empty inputs
  expect_error(guidedPCA(matrix(nrow = 0, ncol = 0), Y, k = 2))
  expect_error(guidedPCA(X, data.frame(), k = 2))
})

test_that("guidedPCA contribution calculation can be disabled", {
  set.seed(777)
  n <- 50
  p <- 30
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  
  result <- guidedPCA(X, Y, k = 2, contribution = FALSE)
  
  expect_null(result$contrib_features)
  expect_null(result$contrib_groups)
  expect_equal(ncol(result$scoreX), 2)  # Other results should still be computed
})