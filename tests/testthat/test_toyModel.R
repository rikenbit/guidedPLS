out <- toyModel()

expect_equal(length(out), 6)
expect_equal(dim(out$X1), c(100, 300))
expect_equal(dim(out$X2), c(200, 150))
expect_equal(dim(out$Y1), c(100, 5))
expect_equal(dim(out$Y1_dummy), c(100, 5))
expect_equal(dim(out$Y2), c(200, 5))
expect_equal(dim(out$Y2_dummy), c(200, 5))
