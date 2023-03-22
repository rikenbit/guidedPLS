\name{softThr}
\docType{methods}
\alias{softThr}

\title{
Soft-thresholding to make a sparse vector sparse
}
\description{
The degree of the sparseness of vector is controlled by the lambda parameter.
}
\usage{
softThr(y, lambda=1)
}

\arguments{
  \item{y}{
A numerical vector.
}
  \item{lambda}{
Threshold value to convert a value 0. If the absolute value of an element of vector is less than lambda, the value is converted to 0 (Default: 1).
}
}

\value{
A numerical vector, whose length is the same as that of y.
}
\author{Koki Tsuyuzaki}

\examples{
y <- seq(-2, 2, 0.1)
softThr(y)
}

\keyword{methods}