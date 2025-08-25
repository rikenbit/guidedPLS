#' Guided PCA (Principal Component Analysis with Label Guidance)
#' 
#' @description
#' Performs guided PCA by finding principal components that maximize covariance 
#' between data matrix X and label/metadata matrix Y. This method extends PLSSVD
#' to automatically handle mixed data types and provide detailed contribution analysis.
#' 
#' @details
#' The algorithm works as follows:
#' 1. **Y preprocessing**: Mixed data types in Y are handled automatically:
#'    - Categorical variables (factor/character) are converted to dummy variables
#'    - Continuous variables (numeric) are used as-is  
#'    - Logical variables are converted to 0/1
#'    - Missing values are handled (NA in factors become a separate category, NA in numerics become 0)
#' 
#' 2. **Normalization**: When normalize_Y=TRUE (default), each Y column is normalized 
#'    to unit L2 norm. This ensures equal weight across different metadata types,
#'    preventing continuous variables with large scales from dominating categorical ones.
#' 
#' 3. **Core computation**: Computes SVD of the cross-product matrix M = X^T Y,
#'    where X is the centered/scaled data matrix and Y is the normalized metadata matrix.
#'    This finds linear combinations that maximize covariance between X and Y.
#' 
#' @param X A numeric matrix (samples x features)
#' @param Y A matrix or data.frame with label/metadata (samples x variables). 
#'          Can contain any mix of numeric (continuous), factor, character (categorical),
#'          or logical columns. Each column type is handled appropriately.
#' @param k Number of components to compute (default: min dimensions)
#' @param center_X Logical, whether to center X columns (default: TRUE)
#' @param scale_X Logical, whether to scale X columns to unit variance (default: TRUE)
#' @param normalize_Y Logical, whether to normalize Y columns to unit L2 norm (default: TRUE).
#'                    This is recommended to balance contributions from different metadata types.
#' @param contribution Logical, whether to calculate feature contributions (default: TRUE)
#' @param deflation Logical, whether to use deflation for sequential component extraction (default: FALSE)
#' @param fullrank Logical, whether to use full SVD or truncated SVD (default: TRUE)
#' @param verbose Logical, whether to print progress messages (default: FALSE)
#' 
#' @return A list containing:
#' \itemize{
#'   \item loadingX: Loading matrix for X (features x components)
#'   \item loadingY: Loading matrix for Y (dummy variables x components)
#'   \item scoreX: Score matrix for X (samples x components)
#'   \item scoreY: Score matrix for Y (samples x components)
#'   \item d: Singular values
#'   \item Y_dummy: The dummy-encoded Y matrix used internally
#'   \item Y_groups: Group labels for dummy variables
#'   \item contrib_features: Feature contributions to each component (if contribution=TRUE)
#'   \item contrib_groups: Grouped contributions by original Y variables (if contribution=TRUE)
#'   \item variance_explained: Variance explained by each component
#' }
#' 
#' @export
#' @examples
#' # Example with mixed data types
#' X <- matrix(rnorm(100*50), 100, 50)
#' Y <- data.frame(
#'   celltype = factor(sample(c("A", "B", "C"), 100, replace=TRUE)),
#'   treatment = factor(sample(c("ctrl", "treated"), 100, replace=TRUE)),
#'   score = rnorm(100)
#' )
#' result <- guidedPCA(X, Y, k=3)
guidedPCA <- function(X, Y, k=NULL, center_X=TRUE, scale_X=TRUE, 
                      normalize_Y=TRUE, contribution=TRUE,
                      deflation=FALSE, fullrank=TRUE, verbose=FALSE) {
  
  # Input validation
  if(verbose) cat("# Input validation...\n")
  if(!is.matrix(X)) X <- as.matrix(X)
  stopifnot(nrow(X) > 0 && ncol(X) > 0)
  
  # Convert Y to dummy matrix
  if(verbose) cat("# Converting Y to dummy variables...\n")
  dummy_result <- .dummify_with_groups(Y)
  Y_dummy <- dummy_result$Y
  Y_groups <- dummy_result$groups
  
  # Remove zero-variance columns
  Y_var <- apply(Y_dummy, 2, stats::var, na.rm=TRUE)
  keep_cols <- Y_var > .Machine$double.eps
  if(sum(keep_cols) == 0) {
    stop("All Y columns have zero variance after dummy encoding")
  }
  Y_dummy <- Y_dummy[, keep_cols, drop=FALSE]
  Y_groups <- Y_groups[keep_cols]
  
  # Set default k
  if(is.null(k)) {
    k <- min(nrow(X)-1, ncol(X), ncol(Y_dummy))
  }
  stopifnot(k > 0 && k <= min(nrow(X)-1, ncol(X), ncol(Y_dummy)))
  
  # Center and scale X
  if(verbose) cat("# Preprocessing X...\n")
  X_processed <- X
  if(center_X) {
    X_processed <- scale(X_processed, center=TRUE, scale=FALSE)
  }
  if(scale_X) {
    X_scales <- apply(X_processed, 2, stats::sd, na.rm=TRUE)
    X_scales[X_scales == 0] <- 1
    X_processed <- sweep(X_processed, 2, X_scales, "/")
  }
  
  # Normalize Y
  if(verbose) cat("# Preprocessing Y...\n")
  Y_processed <- Y_dummy
  if(normalize_Y) {
    Y_norms <- sqrt(colSums(Y_processed^2))
    Y_norms[Y_norms == 0] <- 1
    Y_processed <- sweep(Y_processed, 2, Y_norms, "/")
  }
  
  # Compute cross-product matrix
  if(verbose) cat("# Computing cross-product matrix...\n")
  M <- t(X_processed) %*% Y_processed
  
  # Perform SVD
  if(deflation) {
    if(verbose) cat("# Performing SVD with deflation...\n")
    loadingX <- matrix(0, nrow=ncol(X_processed), ncol=k)
    d <- rep(0, length=k)
    loadingY <- matrix(0, nrow=ncol(Y_processed), ncol=k)
    
    M_deflated <- M
    for(i in seq_len(k)) {
      if(fullrank) {
        res <- svd(M_deflated)
      } else {
        if(!requireNamespace("irlba", quietly=TRUE)) {
          stop("Package 'irlba' is required for fullrank=FALSE")
        }
        res <- irlba::irlba(M_deflated, 1)
      }
      loadingX[, i] <- res$u[, 1]
      d[i] <- res$d[1]
      loadingY[, i] <- res$v[, 1]
      M_deflated <- M_deflated - d[i] * loadingX[, i] %*% t(loadingY[, i])
    }
  } else {
    if(verbose) cat("# Performing SVD...\n")
    if(fullrank) {
      res <- svd(M)
      loadingX <- res$u[, seq_len(k), drop=FALSE]
      d <- res$d[seq_len(k)]
      loadingY <- res$v[, seq_len(k), drop=FALSE]
    } else {
      if(!requireNamespace("irlba", quietly=TRUE)) {
        stop("Package 'irlba' is required for fullrank=FALSE")
      }
      res <- irlba::irlba(M, k)
      loadingX <- res$u
      d <- res$d
      loadingY <- res$v
    }
  }
  
  # Compute scores
  if(verbose) cat("# Computing scores...\n")
  scoreX <- X_processed %*% loadingX
  scoreY <- Y_processed %*% loadingY
  
  # Add column names
  comp_names <- paste0("PC", seq_len(k))
  colnames(loadingX) <- comp_names
  colnames(loadingY) <- comp_names
  colnames(scoreX) <- comp_names
  colnames(scoreY) <- comp_names
  if(!is.null(colnames(X))) rownames(loadingX) <- colnames(X)
  if(!is.null(colnames(Y_dummy))) rownames(loadingY) <- colnames(Y_dummy)
  
  # Calculate variance explained
  total_var <- sum(X_processed^2)
  var_explained <- d^2 / total_var
  names(var_explained) <- comp_names
  
  # Feature contribution analysis
  contrib_features <- NULL
  contrib_groups <- NULL
  
  if(contribution) {
    if(verbose) cat("# Calculating feature contributions...\n")
    
    # Feature contributions to each component
    contrib_features <- matrix(0, nrow=ncol(X), ncol=k)
    rownames(contrib_features) <- colnames(X)
    colnames(contrib_features) <- comp_names
    
    for(i in seq_len(k)) {
      # Contribution of each feature to component i
      # This is based on the loading magnitude and the data variance
      contrib_features[, i] <- abs(loadingX[, i]) * apply(X_processed, 2, stats::var)
      contrib_features[, i] <- contrib_features[, i] / sum(contrib_features[, i])
    }
    
    # Group contributions (from Y side)
    contrib_groups <- matrix(0, nrow=length(unique(Y_groups)), ncol=k)
    rownames(contrib_groups) <- unique(Y_groups)
    colnames(contrib_groups) <- comp_names
    
    for(i in seq_len(k)) {
      # Contribution of each Y variable group to component i
      dummy_contrib <- abs(loadingY[, i])
      for(grp in unique(Y_groups)) {
        contrib_groups[grp, i] <- sum(dummy_contrib[Y_groups == grp])
      }
      contrib_groups[, i] <- contrib_groups[, i] / sum(contrib_groups[, i])
    }
  }
  
  # Return results
  result <- list(
    loadingX = loadingX,
    loadingY = loadingY,
    scoreX = scoreX,
    scoreY = scoreY,
    d = d,
    Y_dummy = Y_dummy,
    Y_groups = Y_groups,
    contrib_features = contrib_features,
    contrib_groups = contrib_groups,
    variance_explained = var_explained
  )
  
  class(result) <- c("guidedPCA", "list")
  return(result)
}

#' Convert data frame to dummy variables with group tracking
#' 
#' @param df A data frame or matrix
#' @return A list with Y (dummy matrix) and groups (group labels)
#' @keywords internal
.dummify_with_groups <- function(df) {
  if (!is.data.frame(df)) df <- as.data.frame(df)
  n <- nrow(df)
  
  mats <- list()
  grps <- list()
  
  for (nm in names(df)) {
    x <- df[[nm]]
    
    if (is.factor(x) || is.character(x)) {
      # Categorical variables
      f <- as.factor(x)
      levs <- levels(f)
      
      # Handle NA as a separate category
      if (any(is.na(f))) {
        f <- addNA(f)
        levs <- levels(f)
      }
      
      if (length(levs) == 1) {
        # Single level - create constant column
        mat <- matrix(1, nrow = n, ncol = 1)
        colnames(mat) <- paste0(nm, "_", levs[1])
      } else {
        # Multiple levels - one-hot encoding (drop first for identifiability)
        mat <- stats::model.matrix(~ f - 1)
        # Remove intercept column if it exists
        if (colnames(mat)[1] == "(Intercept)") {
          mat <- mat[, -1, drop = FALSE]
        }
        colnames(mat) <- paste0(nm, "_", levs)
      }
      
      mats[[length(mats) + 1]] <- mat
      grps[[length(grps) + 1]] <- rep(nm, ncol(mat))
      
    } else if (is.logical(x)) {
      # Logical variables
      mat <- matrix(as.numeric(x), nrow = n, ncol = 1)
      mat[is.na(mat)] <- 0
      colnames(mat) <- nm
      
      mats[[length(mats) + 1]] <- mat
      grps[[length(grps) + 1]] <- nm
      
    } else {
      # Numeric variables
      mat <- matrix(as.numeric(x), nrow = n, ncol = 1)
      mat[is.na(mat)] <- 0
      colnames(mat) <- nm
      
      mats[[length(mats) + 1]] <- mat
      grps[[length(grps) + 1]] <- nm
    }
  }
  
  if (length(mats) == 0) {
    stop("No valid columns in input data")
  }
  
  Y <- do.call(cbind, mats)
  groups <- unlist(grps, use.names = FALSE)
  
  list(Y = Y, groups = groups)
}

#' Print method for guidedPCA
#' 
#' @param x A guidedPCA object
#' @param ... Additional arguments (ignored)
#' @export
print.guidedPCA <- function(x, ...) {
  cat("Guided PCA Results\n")
  cat("==================\n")
  cat("Number of components:", length(x$d), "\n")
  cat("Data dimensions: ", nrow(x$scoreX), "samples x", nrow(x$loadingX), "features\n")
  cat("Label dimensions:", ncol(x$Y_dummy), "dummy variables from", 
      length(unique(x$Y_groups)), "original variables\n")
  cat("\nSingular values:\n")
  print(round(x$d, 4))
  cat("\nVariance explained:\n")
  print(round(x$variance_explained, 4))
  
  if (!is.null(x$contrib_groups)) {
    cat("\nTop contributing label groups per component:\n")
    for (i in seq_len(ncol(x$contrib_groups))) {
      top3 <- head(sort(x$contrib_groups[, i], decreasing = TRUE), 3)
      cat(sprintf("  PC%d: %s\n", i, 
                  paste(sprintf("%s (%.1f%%)", names(top3), top3*100), 
                        collapse = ", ")))
    }
  }
  
  invisible(x)
}

#' Summary method for guidedPCA
#' 
#' @param object A guidedPCA object
#' @param ... Additional arguments (ignored)
#' @export
summary.guidedPCA <- function(object, ...) {
  cat("Guided PCA Summary\n")
  cat("==================\n")
  print(object)
  
  if (!is.null(object$contrib_features)) {
    cat("\nTop contributing features per component:\n")
    for (i in seq_len(ncol(object$contrib_features))) {
      top5 <- head(sort(object$contrib_features[, i], decreasing = TRUE), 5)
      cat(sprintf("\nPC%d:\n", i))
      for (j in seq_along(top5)) {
        cat(sprintf("  %2d. %s (%.1f%%)\n", j, names(top5)[j], top5[j]*100))
      }
    }
  }
  
  invisible(object)
}