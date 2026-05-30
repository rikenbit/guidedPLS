.minDim <- function(...){
    min(unlist(lapply(list(...), dim)))
}

softThr <- function(y, lambda=1){
	if(max(abs(y)) < lambda){
		msg <- paste0(
			"All the elements are converted to 0.\n",
			"Please specify smaller lambda value.")
		stop(msg)
	}
	sign(y) * sapply(abs(y)-lambda, function(x){max(x, 0)})
}

# Check if a matrix is a sparse Matrix (dgCMatrix, dgRMatrix, etc.)
.is_sparse <- function(X) {
    inherits(X, "sparseMatrix")
}

# Compute row means and row standard deviations without densifying
# X: n x p matrix (dense or sparse)
# Returns list(mu, s) where mu = rowMeans, s = rowSds
.row_mean_sd <- function(X) {
    p <- ncol(X)
    if (.is_sparse(X)) {
        mu <- Matrix::rowMeans(X)
        row_sq_sum <- Matrix::rowSums(X^2)
    } else {
        mu <- rowMeans(X)
        row_sq_sum <- rowSums(X^2)
    }
    s <- sqrt((row_sq_sum - p * mu^2) / (p - 1))
    s[s == 0 | is.na(s)] <- 1
    list(mu = mu, s = s)
}

# Compute t(Y_std) %*% X_std without materializing X_std
# X: n x p matrix (sparse or dense)
# Y_std: n x k dense matrix (already row-standardized)
# Returns: k x p dense matrix = t(Y_std) %*% row_standardize(X)
.crossprod_row_standardized <- function(X, Y_std) {
    stats <- .row_mean_sd(X)
    mu <- stats$mu
    s <- stats$s
    # Weight Y by inverse row-sd of X: Y_w[i,a] = Y_std[i,a] / s[i]
    Y_w <- Y_std / s  # n x k, s recycled along columns
    # t(Y_w) %*% X - outer(t(Y_w) %*% mu, rep(1, ncol(X)))
    term1 <- crossprod(Y_w, X)  # k x p
    correction <- as.vector(crossprod(Y_w, mu))  # k-vector
    # Subtract correction from each column (R recycles along columns)
    result <- as.matrix(term1) - correction
    result
}

# Compute t(X_c) %*% Y_c where X_c is row-centered X (centering only, no scaling)
# X: n x p matrix (sparse or dense)
# Y_c: n x k dense matrix (already row-centered)
# Returns: p x k dense matrix
.tcrossprod_row_centered <- function(X, Y_c) {
    if (.is_sparse(X)) {
        mu <- Matrix::rowMeans(X)
    } else {
        mu <- rowMeans(X)
    }
    # t(X_c) %*% Y_c = t(X) %*% Y_c - 1_p %*% t(mu) %*% Y_c
    # where 1_p is a p-vector of ones, and t(mu) %*% Y_c is 1 x k
    term1 <- as.matrix(crossprod(X, Y_c))  # p x k
    mu_Y <- as.vector(crossprod(mu, Y_c))  # k-vector
    # Subtract mu_Y from each row of term1
    p <- nrow(term1)
    correction <- matrix(rep(mu_Y, each = p), nrow = p)  # p x k
    term1 - correction
}

# Compute density of a sparse matrix
.matrix_density <- function(X) {
    if (.is_sparse(X)) {
        Matrix::nnzero(X) / (as.double(nrow(X)) * ncol(X))
    } else {
        1.0
    }
}

dummyMatrix <- function(y, center=TRUE){
    stopifnot(is.numeric(y))
    group <- unique(y)
    l <- length(group)
	out <- matrix(0, nrow=length(y), ncol=l)
	colnames(out) <- group
	for(i in seq_len(l)){        
		out[which(y == group[i]), i] <- 1
	}
    if(center){
    	out <- scale(out, center=TRUE, scale=FALSE)
    }
    out
}
