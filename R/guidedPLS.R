guidedPLS <- function(X1, X2, Y1, Y2, k=.minDim(X1, X2, Y1, Y2),
    cortest=FALSE, fullrank=TRUE, verbose=FALSE){
    # Argument Check
    .checkguidedPLS(X1, X2, Y1, Y2, k, cortest, verbose)
    # Initialization
    int <- .initguidedPLS(X1, X2, Y1, Y2, verbose)
    YX1 <- int$YX1
    YX2 <- int$YX2
    M <- int$M
    # guided PLS
    if(verbose){
        cat("# SVD Step...\n")
    }
    if(fullrank){
        res <- svd(M)
        res$u <- res$u[, seq_len(k)]
        res$d <- res$d[seq_len(k)]
        res$v <- res$v[, seq_len(k)]
    }else{
        res <- irlba(M, k)
    }
    # Loading
    loadingYX1 <- as.matrix(res$u[, seq_len(k)])
    loadingYX2 <- as.matrix(res$v[, seq_len(k)])
    # Score
    scoreX1 <- X1 %*% loadingYX1
    scoreX2 <- X2 %*% loadingYX2
    # Smaller Score
    scoreYX1 <- YX1 %*% loadingYX1
    scoreYX2 <- YX2 %*% loadingYX2
    if(cortest){
        if(verbose){
            cat("# Correlation Test Step...\n")
        }
        # Correlation Coefficient
        corYX1 <- matrix(0, nrow=ncol(YX1), ncol=k)
        corYX2 <- matrix(0, nrow=ncol(YX2), ncol=k)
        for(i in seq_len(k)){
            corYX1[,i] <- apply(YX1, 2, function(x){
                cor(x, scoreYX2[, i])
            })
            corYX2[,i] <- apply(YX2, 2, function(x){
                cor(x, scoreYX1[, i])
            })
        }
        # P-value / Q-value
        pvalYX1 <- matrix(0, nrow=ncol(YX1), ncol=k)
        pvalYX2 <- matrix(0, nrow=ncol(YX2), ncol=k)
        qvalYX1 <- matrix(0, nrow=ncol(YX1), ncol=k)
        qvalYX2 <- matrix(0, nrow=ncol(YX2), ncol=k)
        for(i in seq_len(k)){
            pvalYX1[, i] <- apply(YX1, 2, function(x){
                cor.test(x, scoreYX2[, i])$p.value
            })
            pvalYX2[, i] <- apply(YX2, 2, function(x){
                cor.test(x, scoreYX1[, i])$p.value
            })
            qvalYX1[, i] <- p.adjust(pvalYX1[, i], "BH")
            qvalYX2[, i] <- p.adjust(pvalYX2[, i], "BH")
        }
    }else{
        corYX1 = NULL
        corYX2 = NULL
        pvalYX1 = NULL
        pvalYX2 = NULL
        qvalYX1 = NULL
        qvalYX2 = NULL
    }
    # Output
    list(res=res, loadingYX1=loadingYX1, loadingYX2=loadingYX2,
        scoreX1=scoreX1, scoreX2=scoreX2,
        scoreYX1=scoreYX1, scoreYX2=scoreYX2,
        corYX1=corYX1, corYX2=corYX2,
        pvalYX1=pvalYX1, pvalYX2=pvalYX2, qvalYX1=qvalYX1, qvalYX2=qvalYX2)
}

# Check Function
.checkguidedPLS <- function(X1, X2, Y1, Y2, k, cortest, verbose){
    if(verbose){
        cat("# Input Check Step...\n")
    }
stopifnot(is.matrix(X1))
    stopifnot(is.matrix(X2))
    stopifnot(is.matrix(Y1))
    stopifnot(is.matrix(Y2))
    stopifnot(nrow(X1) == nrow(Y1))
    stopifnot(nrow(X2) == nrow(Y2))
    stopifnot(ncol(Y1) == ncol(Y2))
    stopifnot(is.numeric(k))
    stopifnot(k <= .minDim(X1, X2, Y1, Y2))
    stopifnot(is.logical(cortest))
    stopifnot(is.logical(verbose))
}

# Initialization Function
.initguidedPLS <- function(X1, X2, Y1, Y2, verbose){
    if(verbose){
        cat("# Initialization Step...\n")
    }
    # Auto scaling
    X1 <- t(scale(t(X1), center=TRUE, scale=TRUE))
    X2 <- t(scale(t(X2), center=TRUE, scale=TRUE))
    Y1 <- t(scale(t(Y1), center=TRUE, scale=TRUE))
    Y2 <- t(scale(t(Y2), center=TRUE, scale=TRUE))
    # Auto scaling of Cross-product Matrix
    YX1 <- scale(t(Y1) %*% X1, center=TRUE, scale=TRUE)
    YX2 <- scale(t(Y2) %*% X2, center=TRUE, scale=TRUE)
    # Cross-Cross-product Matrix
    M <- t(YX1) %*% YX2
    list(YX1=YX1, YX2=YX2, M=M)
}
