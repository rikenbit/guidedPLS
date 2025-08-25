sPLSDA <- function(X, Y, k=.minDim(X, Y), cortest=FALSE, lambda=1, thr=1e-10, fullrank=TRUE,
    num.iter=10, verbose=FALSE){
    # Argument Check
	.checksPLSDA(X, Y, k, lambda, thr, fullrank, num.iter, verbose)
	# Initialization
	int <- .initsPLSDA(X, Y, k, verbose)
	X <- int$X
	Y <- int$Y
	M <- int$M
	scoreX <- int$scoreX
	loadingX <- int$loadingX
	scoreY <- int$scoreY
	loadingY <- int$loadingY
	d <- int$d
	for(i in seq_len(k)){
        if(verbose){
            cat(paste0("Outer iter = ", i, "\n"))
			cat("# SVD Step...\n")
        }
		if(fullrank){
			res <- svd(M)
		}else{
			res <- irlba(M, 1)
		}
		loadingX[, i] <-res$u[, 1]
		d[i] <- res$d[1]
		loadingY[, i] <- res$v[, 1]
		r <- 10
		iter <- 1
        if(verbose){
			cat("# Soft-thresholding Step...\n")
        }
		while((r > thr) && (iter <= num.iter)){
			# Soft-thresholding
			loadingX_old <- softThr(M %*% loadingY[, i], lambda)
			loadingY_old <- softThr(t(M) %*% loadingX[, i], lambda)
			# Normalization
			loadingX[, i] <- loadingX_old / sqrt(sum(loadingX_old^2))
			loadingY[, i] <- loadingY_old / sqrt(sum(loadingY_old^2))
			# Convergence
			r <- max(sqrt(sum(loadingX[,i] - loadingX_old)^2), sqrt(sum(loadingY[,i] - loadingY_old)^2))
			if(verbose){
				cat(paste0("Inner iter = ", iter, ": r = ", r, "\n"))
			}
			iter <- iter + 1
		}
		# Score
		scoreX[, i] <- X %*% loadingX[, i]
		scoreY[, i] <- Y %*% loadingY[, i]
		# Deflation
		c_i <- t(X) %*% (scoreX[, i] / as.numeric(t(scoreX[, i]) %*% scoreX[, i]))
		d_i <- t(Y) %*% (scoreX[, i] / as.numeric(t(scoreX[, i]) %*% scoreY[, i]))
		X <- X - scoreX[, i] %*% t(c_i)
		Y <- Y - scoreX[, i] %*% t(d_i)
		M <- t(X) %*% Y
	}
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
.checksPLSDA <- function(X, Y, k, lambda, thr, fullrank, num.iter, verbose){
    if(verbose){
        cat("# Input Check Step...\n")
    }
    stopifnot(is.matrix(X))
    stopifnot(is.matrix(Y))
    stopifnot(nrow(X) == nrow(Y))
    stopifnot(is.numeric(k))
    stopifnot(k <= .minDim(X, Y))
    stopifnot(is.numeric(lambda))
    stopifnot(lambda >= 0)
    stopifnot(is.numeric(thr))
    stopifnot(thr >= 0)
    stopifnot(is.logical(fullrank))
    stopifnot(is.numeric(num.iter))
    stopifnot(num.iter >= 0)
    stopifnot(is.logical(verbose))
}

# Initialization Function
.initsPLSDA <- function(X, Y, k, verbose){
    if(verbose){
        cat("# Initialization Step...\n")
    }
	# Auto scaling
    X <- scale(X, center=TRUE, scale=TRUE)
    Y <- scale(Y, center=TRUE, scale=TRUE)
    # Cross-product Matrix
    M <- t(X) %*% Y
	# Output object for SVD
	scoreX <- matrix(0, nrow=nrow(X), ncol=k)
	loadingX <- matrix(0, nrow=ncol(X), ncol=k)
	d <- rep(0, length=k)
	scoreY <- matrix(0, nrow=nrow(Y), ncol=k)
	loadingY <- matrix(0, nrow=ncol(Y), ncol=k)
    list(X=X, Y=Y, M=M, scoreX=scoreX, loadingX=loadingX, scoreY=scoreY, loadingY=loadingY, d=d)
}
