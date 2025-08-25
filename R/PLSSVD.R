PLSSVD <- function(X, Y, k=.minDim(X, Y), cortest=FALSE, deflation=FALSE, fullrank=TRUE,
    verbose=FALSE){
    # Argument Check
    .checkPLSSVD(X, Y, k, verbose)
	# Initialization
	int <- .initPLSSVD(X, Y, verbose)
	X <- int$X
	Y <- int$Y
    M <- int$M
    if(deflation){
        if(verbose){
            cat("# SVD w deflation Step...\n")
        }
        loadingX <- matrix(0, nrow=ncol(X), ncol=k)
        d <- rep(0, length=k)
        loadingY <- matrix(0, nrow=ncol(Y), ncol=k)
        for(i in seq_len(k)){
			if(fullrank){
	            res <- svd(M)
			}else{
				res <- irlba(M, 1)
			}
            loadingX[, i] <-res$u[, 1]
            d[i] <- res$d[1]
            loadingY[, i] <- res$v[, 1]
            M <- M - d[i] * loadingX[, i] %*% t(loadingY[, i])
        }
    }else{
        if(verbose){
            cat("# SVD w/o deflation Step...\n")
        }
		if(fullrank){
			res <- svd(M)
		}else{
	        res <- irlba(M, k)
		}
		loadingX <- res$u[, seq_len(k)]
		d <- res$d[seq_len(k)]
		loadingY <- res$v[, seq_len(k)]
    }
    # Projection
    scoreX <- X %*% loadingX
    scoreY <- Y %*% loadingY
    if(cortest){
        if(verbose){
            cat("# Correlation Test Step...\n")
        }
        # Correlation Coefficient
        corX <- matrix(0, nrow=ncol(X), ncol=k)
        corY <- matrix(0, nrow=ncol(Y), ncol=k)
        for(i in seq_len(k)){
            corX[,i] <- apply(X, 2, function(x){
                cor(x, scoreY[, i])
            })
            corY[,i] <- apply(Y, 2, function(x){
                cor(x, scoreX[, i])
            })
        }
        # P-value / Q-value
        pvalX <- matrix(0, nrow=ncol(X), ncol=k)
        qvalX <- matrix(0, nrow=ncol(X), ncol=k)
        pvalY <- matrix(0, nrow=ncol(Y), ncol=k)
        qvalY <- matrix(0, nrow=ncol(Y), ncol=k)
        for(i in seq_len(k)){
            pvalX[, i] <- apply(X, 2, function(x){
                cor.test(x, scoreY[, i])$p.value
            })
            pvalY[, i] <- apply(Y, 2, function(x){
                cor.test(x, scoreX[, i])$p.value
            })
            qvalX[, i] <- p.adjust(pvalX[, i], "BH")
            qvalY[, i] <- p.adjust(pvalY[, i], "BH")
        }
    }else{
        corX = NULL
        corY = NULL
        pvalX = NULL
        pvalY = NULL
        qvalX = NULL
        qvalY = NULL
    }
    # Output
    list(scoreX=scoreX, loadingX=loadingX, scoreY=scoreY, loadingY=loadingY, d=d, corX=corX, corY=corY, pvalX=pvalX, pvalY=pvalY, qvalX=qvalX, qvalY=qvalY)
}

# Check Function
.checkPLSSVD <- function(X, Y, k, verbose){
    if(verbose){
        cat("# Input Check Step...\n")
    }
    stopifnot(is.matrix(X))
    stopifnot(is.matrix(Y))
    stopifnot(nrow(X) == nrow(Y))
    stopifnot(is.numeric(k))
    stopifnot(k <= .minDim(X, Y))
    stopifnot(is.logical(verbose))
}

# Initialization Function
.initPLSSVD <- function(X, Y, verbose){
    if(verbose){
        cat("# Initialization Step...\n")
    }
	# Auto scaling
    X <- scale(X, center=TRUE, scale=TRUE)
    Y <- scale(Y, center=TRUE, scale=TRUE)
    # Cross-product Matrix
    M <- t(X) %*% Y
    list(X=X, Y=Y, M=M)
}
