set.seed(1234)
X <- matrix(rbinom(10*20, 10, 0.1), nrow=10, ncol=20)
Y <- matrix(rbinom(10*15, 10, 0.1), nrow=10, ncol=15)

out1 <- sPLSDA(X, Y, k=2, thr=1, verbose=TRUE)
out2 <- sPLSDA(X, Y, k=2, thr=1, fullrank=TRUE, verbose=TRUE)

expect_equal(dim(out1$u), c(ncol(X), 2))
expect_equal(length(out1$d), 2)
expect_equal(dim(out1$v), c(ncol(Y), 2))

expect_equal(dim(out2$u), c(ncol(X), 2))
expect_equal(length(out2$d), 2)
expect_equal(dim(out2$v), c(ncol(Y), 2))
