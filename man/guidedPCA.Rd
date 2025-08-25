\name{guidedPCA}
\alias{guidedPCA}
\title{Guided PCA (Principal Component Analysis with Label Guidance)}
\description{
Performs guided PCA by finding principal components that maximize covariance 
between data matrix X and label/metadata matrix Y. This method extends PLSSVD
to automatically handle mixed data types and provide detailed contribution analysis.
}
\usage{
guidedPCA(X, Y, k = NULL, center_X = TRUE, scale_X = TRUE, 
          normalize_Y = TRUE, contribution = TRUE, 
          deflation = FALSE, fullrank = TRUE, verbose = FALSE)
}
\arguments{
\item{X}{A numeric matrix (samples x features)}
\item{Y}{A matrix or data.frame with label/metadata (samples x variables). 
         Can contain any mix of numeric (continuous), factor, character (categorical),
         or logical columns. Each column type is handled appropriately.}
\item{k}{Number of components to compute (default: min dimensions)}
\item{center_X}{Logical, whether to center X columns (default: TRUE)}
\item{scale_X}{Logical, whether to scale X columns to unit variance (default: TRUE)}
\item{normalize_Y}{Logical, whether to normalize Y columns to unit L2 norm (default: TRUE).
                   This is recommended to balance contributions from different metadata types.}
\item{contribution}{Logical, whether to calculate feature contributions (default: TRUE)}
\item{deflation}{Logical, whether to use deflation for sequential component extraction (default: FALSE)}
\item{fullrank}{Logical, whether to use full SVD or truncated SVD (default: TRUE)}
\item{verbose}{Logical, whether to print progress messages (default: FALSE)}
}
\details{
The algorithm works as follows:

1. Y preprocessing: Mixed data types in Y are handled automatically:
   - Categorical variables (factor/character) are converted to dummy variables
   - Continuous variables (numeric) are used as-is  
   - Logical variables are converted to 0/1
   - Missing values are handled (NA in factors become a separate category, NA in numerics become 0)

2. Normalization: When normalize_Y=TRUE (default), each Y column is normalized 
   to unit L2 norm. This ensures equal weight across different metadata types,
   preventing continuous variables with large scales from dominating categorical ones.

3. Core computation: Computes SVD of the cross-product matrix M = X^T Y,
   where X is the centered/scaled data matrix and Y is the normalized metadata matrix.
   This finds linear combinations that maximize covariance between X and Y.
}
\value{
A list of class "guidedPCA" containing:
\itemize{
  \item loadingX: Loading matrix for X (features x components)
  \item loadingY: Loading matrix for Y (dummy variables x components)
  \item scoreX: Score matrix for X (samples x components)
  \item scoreY: Score matrix for Y (samples x components)
  \item d: Singular values
  \item Y_dummy: The dummy-encoded Y matrix used internally
  \item Y_groups: Group labels for dummy variables
  \item contrib_features: Feature contributions to each component (if contribution=TRUE)
  \item contrib_groups: Grouped contributions by original Y variables (if contribution=TRUE)
  \item variance_explained: Variance explained by each component
}
}
\examples{
# Example with mixed data types
X <- matrix(rnorm(100*50), 100, 50)
Y <- data.frame(
  celltype = factor(sample(c("A", "B", "C"), 100, replace=TRUE)),
  treatment = factor(sample(c("ctrl", "treated"), 100, replace=TRUE)),
  score = rnorm(100)
)
result <- guidedPCA(X, Y, k=3)
print(result)
summary(result)
}
\author{
Koki Tsuyuzaki
}
\references{
Reese S E, et al. A new statistic for identifying batch effects in high-throughput 
genomic data that uses guided principal component analysis. Bioinformatics, 
29(22), 2877-2883, 2013
}