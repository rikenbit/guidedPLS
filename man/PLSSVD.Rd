\name{PLSSVD}
\docType{methods}
\alias{PLSSVD}

\title{
Partial Least Squares by Singular Value Decomposition (PLS-SVD)
}
\description{
Two matrices X and Y sharing a row are required
}
\usage{
PLSSVD(X, Y, k=.minDim(X, Y), cortest=FALSE,
  deflation=FALSE, fullrank=TRUE, verbose=FALSE)
}

\arguments{
\item{X}{
The input matrix which has N-rows and M-columns.
}
\item{Y}{
The input matrix which has N-rows and L-columns.
}
\item{k}{
The number of low-dimension (k < {N, M, L}, Default: .minDim(X, Y))
}
\item{cortest}{
If cortest is set as TRUE, t-test of correlation coefficient is performed (Default: FALSE)
}
\item{deflation}{
If deflation is set as TRUE, the score vectors are made orthogonal, otherwise the loading vectors are made orthogonal (Default: FALSE)
}
\item{fullrank}{
If fullrank is set as TRUE, irlba is used, otherwise fullrank SVD is used (Default: TRUE)
}
\item{verbose}{
Verbose option (Default: FALSE)
}
}

\value{
  scoreX : Score matrix which has M-rows and K-columns.
  loadingX : Loading matrix which has N-rows and K-columns.
  scoreY : Score matrix which has L-rows and K-columns.
  loadingY : Loading matrix which has N-rows and K-columns.
  d : K-length singular value vector of the cross-product matrix X'Y.
  corX: Correlation Coefficient (Default: NULL)
  corY: Correlation Coefficient (Default: NULL)
  pvalX: P-value vector of corX (Default: NULL)
  pvalY: P-value vector of corY (Default: NULL)
  qvalX: Q-value vector of BH method against pvalX (Default: NULL)
  qvalY: Q-value vector of BH method against pvalY (Default: NULL)
}

\references{
Le Cao, et al. (2008). A Sparse PLS for Variable Selection when Integrating Omics Data. \emph{Statistical Applications in Genetics and Molecular Biology}, 7(1)
}
\author{Koki Tsuyuzaki}

\examples{
# Test data
data <- toyModel()

# Simple usage
out <- PLSSVD(X=data$X1, Y=data$Y1, k=4)
}
\keyword{methods}