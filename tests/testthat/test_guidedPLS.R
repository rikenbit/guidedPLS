X1 <- matrix(runif(20*10), nrow=20, ncol=10)
X2 <- matrix(runif(15*25), nrow=15, ncol=25)
Y1 <- matrix(runif(20*5), nrow=20, ncol=5)
Y2 <- matrix(runif(15*5), nrow=15, ncol=5)

out1 <- guidedPLS(X1, X2, Y1, Y2, k=2)

expect_equal(length(out1$res$d), 2)
expect_equal(dim(out1$res$u), c(ncol(X1), 2))
expect_equal(dim(out1$res$v), c(ncol(X2), 2))
expect_equal(dim(out1$loading1), c(ncol(X1), 2))
expect_equal(dim(out1$loading2), c(ncol(X2), 2))
expect_equal(dim(out1$score1), c(nrow(X1), 2))
expect_equal(dim(out1$score2), c(nrow(X2), 2))
expect_equal(dim(out1$score_cor1), c(ncol(Y1), 2))
expect_equal(dim(out1$score_cor2), c(ncol(Y2), 2))

out2 <- guidedPLS(X1, X2, Y1, Y2, k=3, fullrank=TRUE)

expect_equal(length(out2$res$d), 3)
expect_equal(dim(out2$res$u), c(ncol(X1), 3))
expect_equal(dim(out2$res$v), c(ncol(X2), 3))
expect_equal(dim(out2$loading1), c(ncol(X1), 3))
expect_equal(dim(out2$loading2), c(ncol(X2), 3))
expect_equal(dim(out2$score1), c(nrow(X1), 3))
expect_equal(dim(out2$score2), c(nrow(X2), 3))
expect_equal(dim(out2$score_cor1), c(ncol(Y1), 3))
expect_equal(dim(out2$score_cor2), c(ncol(Y2), 3))

out3 <- guidedPLS(X1, X2, Y1, Y2, k=2, cortest=TRUE)

expect_equal(dim(out3$cor1), c(ncol(X1), 2))
expect_equal(dim(out3$pval1), c(ncol(X1), 2))
expect_equal(dim(out3$qval1), c(ncol(X1), 2))
expect_equal(dim(out3$cor2), c(ncol(X2), 2))
expect_equal(dim(out3$pval2), c(ncol(X2), 2))
expect_equal(dim(out3$qval2), c(ncol(X2), 2))
