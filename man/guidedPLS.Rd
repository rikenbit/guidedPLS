\name{guidedPLS}
\docType{methods}
\alias{guidedPLS}

\title{
Guided Partial Least Squares (guied-PLS)
}
\description{
Four matrices X1, X2, Y1, and Y2 are required.
X1 and Y1 are supposed to share the rows,
X2 and Y2 are supposed to share the rows,
and Y1 and Y2 are supposed to share the columns.
}
\usage{
guidedPLS(X1, X2, Y1, Y2, k=.minDim(X1, X2, Y1, Y2),
    cortest=FALSE, fullrank=TRUE, sumcor=FALSE, lambda=1e-6, verbose=FALSE)
}

\arguments{
\item{X1}{
The input matrix which has N-rows and M-columns.
}
\item{Y1}{
The input matrix which has N-rows and L-columns.
}
\item{X2}{
The input matrix which has O-rows and P-columns.
}
\item{Y2}{
The input matrix which has O-rows and L-columns.
}
\item{k}{
The number of low-dimension (k < \{N, M, L, O\}, Default: .minDim(X1, X2, Y1, Y2))
}
\item{cortest}{
If cortest is set as TRUE, t-test of correlation coefficient is performed (Default: FALSE)
}
\item{fullrank}{
If fullrank is set as TRUE, irlba is used, otherwise fullrank SVD is used (Default: TRUE)
}
\item{sumcor}{
If sumcor is set as TRUE, SUMCOR-based CCA using generalized eigenvalue decomposition is performed instead of SVD-based approach. This maximizes Tr(W1^T X1^T Y1 Y2^T X2 W2) subject to W1^T X1^T X1 W1 = I and W2^T X2^T X2 W2 = I (Default: FALSE). Requires the geigen package.
}
\item{lambda}{
Regularization parameter for numerical stability in SUMCOR-based CCA. Only used when sumcor=TRUE (Default: 1e-6). Larger values provide more regularization.
}
\item{verbose}{
Verbose option (Default: FALSE)
}
}

\value{
  res: object of svd()
  loadingYX1: Loading vector to project X1 to lower dimension via Y1 (M times k).
  loadingYX2: Loading vector to project X2 to lower dimension via Y2 (P times k).
  scoreX1: Projected X1 (N times k)
  scoreX2: Projected X2 (O times k)
  scoreYX1: Projected YX1 (L times k)
  scoreYX2: Projected YX2 (L times k)
  corYX1: Correlation Coefficient (Default: NULL)
  corYX2: Correlation Coefficient (Default: NULL)
  pvalYX1: P-value vector of corYX1 (Default: NULL)
  pvalYX2: P-value vector of corYX2 (Default: NULL)
  qvalYX1: Q-value vector of BH method against pvalYX1 (Default: NULL)
  qvalYX2: Q-value vector of BH method against pvalYX2 (Default: NULL)
}

\references{
Le Cao, et al. (2008). A Sparse PLS for Variable Selection when Integrating Omics Data. \emph{Statistical Applications in Genetics and Molecular Biology}, 7(1)

Reese S E, et al. (2013). A new statistic for identifying batch effects in high-throughput genomic data that uses guided principal component analysis. \emph{Bioinformatics}, 29(22), 2877-2883
}
\author{Koki Tsuyuzaki}

\examples{
# Test data
data <- toyModel()

# Simple usage
out <- guidedPLS(X1=data$X1, X2=data$X2, Y1=data$Y1, Y2=data$Y2, k=4)
}
\keyword{methods}