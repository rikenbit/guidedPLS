guidedPLS <- function(X1, X2, Y1, Y2, k=.minDim(X1, X2, Y1, Y2),
    cortest=FALSE, fullrank=TRUE, sumcor=FALSE, lambda=1e-6, verbose=FALSE){
    # Argument Check
    .checkguidedPLS(X1, X2, Y1, Y2, k, cortest, sumcor, lambda, verbose)
    # Initialization
    int <- .initguidedPLS(X1, X2, Y1, Y2, verbose)
    YX1 <- int$YX1
    YX2 <- int$YX2
    M <- int$M
    
    if(sumcor){
        # SUMCOR-based CCA implementation
        if(verbose){
            cat("# SUMCOR-based CCA Step...\n")
        }
        # Check if geigen is available
        if(!requireNamespace("geigen", quietly = TRUE)){
            stop("Package 'geigen' is required for sumcor=TRUE. Please install it.")
        }
        
        # Compute generalized eigenvalue decomposition
        # Objective: max Tr(W1^T X1^T Y1 Y2^T X2 W2)
        # Constraints: W1^T X1^T X1 W1 = I, W2^T X2^T X2 W2 = I
        
        # Center the data (row-wise centering: subtract row means)
        Y1_centered <- t(scale(t(Y1), center=TRUE, scale=FALSE))
        Y2_centered <- t(scale(t(Y2), center=TRUE, scale=FALSE))

        # Dimensions
        p1 <- ncol(X1)
        p2 <- ncol(X2)
        n1 <- nrow(X1)
        n2 <- nrow(X2)

        # SUMCOR requires p x p covariance matrices, so densify sparse X
        if(.is_sparse(X1)){
            if(verbose) cat("# SUMCOR: converting sparse X1 to dense (p x p covariance needed)\n")
            X1 <- as.matrix(X1)
        }
        if(.is_sparse(X2)){
            if(verbose) cat("# SUMCOR: converting sparse X2 to dense (p x p covariance needed)\n")
            X2 <- as.matrix(X2)
        }
        X1_centered <- t(scale(t(X1), center=TRUE, scale=FALSE))
        X2_centered <- t(scale(t(X2), center=TRUE, scale=FALSE))

        # Build the concatenated cross-covariance matrix C
        # C = [0, C12; C21, 0] where C12 = X1^T Y1 Y2^T X2 / n
        C12 <- (t(X1_centered) %*% Y1_centered %*% t(Y2_centered) %*% X2_centered) / (n1 * n2)

        # Build the concatenated covariance matrix D (block diagonal)
        # D = [D11, 0; 0, D22] where D11 = X1^T X1 / n, D22 = X2^T X2 / n
        D11 <- (t(X1_centered) %*% X1_centered) / n1
        D22 <- (t(X2_centered) %*% X2_centered) / n2

        C21 <- t(C12)

        # Create block matrix C
        C <- matrix(0, nrow=p1+p2, ncol=p1+p2)
        C[1:p1, (p1+1):(p1+p2)] <- C12
        C[(p1+1):(p1+p2), 1:p1] <- C21

        # Add regularization for numerical stability
        D11_reg <- D11 + lambda * diag(p1)
        D22_reg <- D22 + lambda * diag(p2)
        
        # Create block diagonal matrix D
        D <- matrix(0, nrow=p1+p2, ncol=p1+p2)
        D[1:p1, 1:p1] <- D11_reg
        D[(p1+1):(p1+p2), (p1+1):(p1+p2)] <- D22_reg
        
        # Solve generalized eigenvalue problem: C*W = D*W*Lambda
        eigen_result <- geigen::geigen(C, D, symmetric=FALSE)
        
        # Select top k components (by absolute value of eigenvalues)
        idx <- order(abs(eigen_result$values), decreasing=TRUE)[seq_len(k)]
        
        # Extract the concatenated eigenvectors
        W <- eigen_result$vectors[, idx, drop=FALSE]
        
        # Split W into W1 and W2
        loadingYX1 <- as.matrix(W[1:p1, , drop=FALSE])
        loadingYX2 <- as.matrix(W[(p1+1):(p1+p2), , drop=FALSE])
        
        # Normalize each loading vector to have unit norm
        for(i in seq_len(k)){
            norm1 <- sqrt(sum(loadingYX1[,i]^2))
            norm2 <- sqrt(sum(loadingYX2[,i]^2))
            if(norm1 > 0) loadingYX1[,i] <- loadingYX1[,i] / norm1
            if(norm2 > 0) loadingYX2[,i] <- loadingYX2[,i] / norm2
        }
        
        # Create res object for compatibility
        res <- list(u=loadingYX1, v=loadingYX2, d=abs(eigen_result$values[idx]))
        
    }else{
        # Original SVD-based implementation
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
    }
    
    # Score
    scoreX1 <- as.matrix(X1 %*% loadingYX1)
    scoreX2 <- as.matrix(X2 %*% loadingYX2)
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
.checkguidedPLS <- function(X1, X2, Y1, Y2, k, cortest, sumcor, lambda, verbose){
    if(verbose){
        cat("# Input Check Step...\n")
    }
    stopifnot(is.matrix(X1) || .is_sparse(X1))
    stopifnot(is.matrix(X2) || .is_sparse(X2))
    stopifnot(is.matrix(Y1))
    stopifnot(is.matrix(Y2))
    stopifnot(nrow(X1) == nrow(Y1))
    stopifnot(nrow(X2) == nrow(Y2))
    stopifnot(ncol(Y1) == ncol(Y2))
    stopifnot(is.numeric(k))
    stopifnot(k <= .minDim(X1, X2, Y1, Y2))
    stopifnot(is.logical(cortest))
    stopifnot(is.logical(sumcor))
    stopifnot(is.numeric(lambda))
    stopifnot(lambda >= 0)
    stopifnot(is.logical(verbose))
}

# Initialization Function
.initguidedPLS <- function(X1, X2, Y1, Y2, verbose){
    if(verbose){
        cat("# Initialization Step...\n")
    }
    # Y matrices are always small (samples x k), so row-standardize normally
    Y1 <- t(scale(t(Y1), center=TRUE, scale=TRUE))
    Y2 <- t(scale(t(Y2), center=TRUE, scale=TRUE))
    # Cross-product with row-standardized X
    use_sparse <- .is_sparse(X1) || .is_sparse(X2)
    if(use_sparse){
        if(verbose){
            cat("# Sparse-aware path for X matrices...\n")
            if(.is_sparse(X1)){
                cat(sprintf("#   X1: %d x %d, density=%.4f\n",
                    nrow(X1), ncol(X1), .matrix_density(X1)))
            }
            if(.is_sparse(X2)){
                cat(sprintf("#   X2: %d x %d, density=%.4f\n",
                    nrow(X2), ncol(X2), .matrix_density(X2)))
            }
        }
        # Compute t(Y_std) %*% X_std without materializing dense X_std
        YX1 <- .crossprod_row_standardized(X1, Y1)
        YX2 <- .crossprod_row_standardized(X2, Y2)
    }else{
        # Original dense path
        X1 <- t(scale(t(X1), center=TRUE, scale=TRUE))
        X2 <- t(scale(t(X2), center=TRUE, scale=TRUE))
        YX1 <- t(Y1) %*% X1
        YX2 <- t(Y2) %*% X2
    }
    # Auto scaling of Cross-product Matrix (always small: k x p)
    YX1 <- scale(YX1, center=TRUE, scale=TRUE)
    YX2 <- scale(YX2, center=TRUE, scale=TRUE)
    # Cross-Cross-product Matrix
    M <- t(YX1) %*% YX2
    list(YX1=YX1, YX2=YX2, M=M)
}
