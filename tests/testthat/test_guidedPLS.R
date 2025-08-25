X1 <- matrix(runif(20*10), nrow=20, ncol=10)
X2 <- matrix(runif(15*25), nrow=15, ncol=25)
Y1 <- matrix(runif(20*5), nrow=20, ncol=5)
Y2 <- matrix(runif(15*5), nrow=15, ncol=5)

out1 <- guidedPLS(X1, X2, Y1, Y2, k=2)

expect_equal(length(out1$res$d), 2)
expect_equal(dim(out1$res$u), c(ncol(X1), 2))
expect_equal(dim(out1$res$v), c(ncol(X2), 2))
expect_equal(dim(out1$loadingYX1), c(ncol(X1), 2))
expect_equal(dim(out1$loadingYX2), c(ncol(X2), 2))
expect_equal(dim(out1$scoreX1), c(nrow(X1), 2))
expect_equal(dim(out1$scoreX2), c(nrow(X2), 2))
expect_equal(dim(out1$scoreYX1), c(ncol(Y1), 2))
expect_equal(dim(out1$scoreYX2), c(ncol(Y2), 2))

out2 <- guidedPLS(X1, X2, Y1, Y2, k=3, fullrank=TRUE)

expect_equal(length(out2$res$d), 3)
expect_equal(dim(out2$res$u), c(ncol(X1), 3))
expect_equal(dim(out2$res$v), c(ncol(X2), 3))
expect_equal(dim(out2$loadingYX1), c(ncol(X1), 3))
expect_equal(dim(out2$loadingYX2), c(ncol(X2), 3))
expect_equal(dim(out2$scoreX1), c(nrow(X1), 3))
expect_equal(dim(out2$scoreX2), c(nrow(X2), 3))
expect_equal(dim(out2$scoreYX1), c(ncol(Y1), 3))
expect_equal(dim(out2$scoreYX2), c(ncol(Y2), 3))

out3 <- guidedPLS(X1, X2, Y1, Y2, k=2, cortest=TRUE)

expect_equal(length(out3$res$d), 2)
expect_equal(dim(out3$res$u), c(ncol(X1), 2))
expect_equal(dim(out3$res$v), c(ncol(X2), 2))
expect_equal(dim(out3$loadingYX1), c(ncol(X1), 2))
expect_equal(dim(out3$loadingYX2), c(ncol(X2), 2))
expect_equal(dim(out3$scoreX1), c(nrow(X1), 2))
expect_equal(dim(out3$scoreX2), c(nrow(X2), 2))
expect_equal(dim(out3$scoreYX1), c(ncol(Y1), 2))
expect_equal(dim(out3$scoreYX2), c(ncol(Y2), 2))

expect_equal(dim(out3$corYX1), c(ncol(X1), 2))
expect_equal(dim(out3$corYX2), c(ncol(X2), 2))
expect_equal(dim(out3$pvalYX1), c(ncol(X1), 2))
expect_equal(dim(out3$pvalYX2), c(ncol(X2), 2))
expect_equal(dim(out3$qvalYX1), c(ncol(X1), 2))
expect_equal(dim(out3$qvalYX2), c(ncol(X2), 2))

# Test SUMCOR-based CCA implementation
if(requireNamespace("geigen", quietly = TRUE)){
    out4 <- guidedPLS(X1, X2, Y1, Y2, k=2, sumcor=TRUE, lambda=1e-6)
    
    expect_equal(length(out4$res$d), 2)
    expect_equal(dim(out4$res$u), c(ncol(X1), 2))
    expect_equal(dim(out4$res$v), c(ncol(X2), 2))
    expect_equal(dim(out4$loadingYX1), c(ncol(X1), 2))
    expect_equal(dim(out4$loadingYX2), c(ncol(X2), 2))
    expect_equal(dim(out4$scoreX1), c(nrow(X1), 2))
    expect_equal(dim(out4$scoreX2), c(nrow(X2), 2))
    expect_equal(dim(out4$scoreYX1), c(ncol(Y1), 2))
    expect_equal(dim(out4$scoreYX2), c(ncol(Y2), 2))
    
    # Test that loadings are normalized
    for(i in 1:2){
        expect_equal(sum(out4$loadingYX1[,i]^2), 1, tolerance=1e-10)
        expect_equal(sum(out4$loadingYX2[,i]^2), 1, tolerance=1e-10)
    }
    
    # Test SUMCOR with correlation test
    out5 <- guidedPLS(X1, X2, Y1, Y2, k=3, sumcor=TRUE, cortest=TRUE, lambda=1e-4)
    
    expect_equal(length(out5$res$d), 3)
    expect_equal(dim(out5$res$u), c(ncol(X1), 3))
    expect_equal(dim(out5$res$v), c(ncol(X2), 3))
    expect_equal(dim(out5$loadingYX1), c(ncol(X1), 3))
    expect_equal(dim(out5$loadingYX2), c(ncol(X2), 3))
    expect_equal(dim(out5$scoreX1), c(nrow(X1), 3))
    expect_equal(dim(out5$scoreX2), c(nrow(X2), 3))
    expect_equal(dim(out5$scoreYX1), c(ncol(Y1), 3))
    expect_equal(dim(out5$scoreYX2), c(ncol(Y2), 3))
    
    expect_equal(dim(out5$corYX1), c(ncol(X1), 3))
    expect_equal(dim(out5$corYX2), c(ncol(X2), 3))
    expect_equal(dim(out5$pvalYX1), c(ncol(X1), 3))
    expect_equal(dim(out5$pvalYX2), c(ncol(X2), 3))
    expect_equal(dim(out5$qvalYX1), c(ncol(X1), 3))
    expect_equal(dim(out5$qvalYX2), c(ncol(X2), 3))
}
