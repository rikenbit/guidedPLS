\name{toyModel}
\docType{methods}
\alias{toyModel}

\title{
Toy model data for using PLSSVD, sPLSDA, and guidedPLS
}
\description{
The data is used for confirming the algorithm are properly working.
}
\usage{
toyModel(model="Easy", seeds=123)
}

\arguments{
  \item{model}{
"Easy" and "Hard" are available (Default: "Easy").
}
  \item{seeds}{
Random number for setting set.seeds in the function (Default: 123).
}
}

\value{
A list object containing a set of matrices X1, X2, Y1, Y1_dummy, Y2, Y1_dummy.
}
\author{Koki Tsuyuzaki}

\seealso{\code{\link{PLSSVD}},\code{\link{sPLSDA}},\code{\link{guidedPLS}}}
\examples{
data <- toyModel(seeds=123)
}

\keyword{methods}