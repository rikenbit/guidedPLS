\name{dummyMatrix}
\docType{methods}
\alias{dummyMatrix}

\title{
Toy model data for using dNMF, dSVD, dsiNMF, djNMF, dPLS, dNTF, and dNTD
}
\description{
A label vector is converted to a dummy matrix.
}
\usage{
dummyMatrix(y, center=TRUE)
}

\arguments{
  \item{y}{
A label vector to specify the group of data.
}
  \item{center}{
An option to center the rows of matrix (Default: TRUE).
}
}

\value{
A matrix is generated. The number of row is equal to the length of y and the number of columns is the number of unique elements of y.
}
\author{Koki Tsuyuzaki}

\examples{
y <- c(1, 3, 2, 1, 4, 2)
dummyMatrix(y)
}

\keyword{methods}