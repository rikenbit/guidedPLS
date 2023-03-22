X <- matrix(runif(10*20), nrow=10, ncol=20)
Y <- matrix(runif(10*15), nrow=10, ncol=15)

out1 <- PLSSVD(X, Y, k=2)
out2 <- PLSSVD(X, Y, deflation=TRUE, k=2)
out3 <- PLSSVD(X, Y, fullrank=TRUE, k=2)
out4 <- PLSSVD(X, Y, deflation=TRUE, fullrank=TRUE, k=2)

expect_equal(dim(out1$u), c(ncol(X), 2))
expect_equal(length(out1$d), 2)
expect_equal(dim(out1$v), c(ncol(Y), 2))

expect_equal(dim(out2$u), c(ncol(X), 2))
expect_equal(length(out2$d), 2)
expect_equal(dim(out2$v), c(ncol(Y), 2))

expect_equal(dim(out3$u), c(ncol(X), 2))
expect_equal(length(out3$d), 2)
expect_equal(dim(out3$v), c(ncol(Y), 2))

expect_equal(dim(out4$u), c(ncol(X), 2))
expect_equal(length(out4$d), 2)
expect_equal(dim(out4$v), c(ncol(Y), 2))
