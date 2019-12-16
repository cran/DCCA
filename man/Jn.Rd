\name{Jn}
\alias{Jn}
\title{Matrix J}
\usage{
Jn(n = 2)
}
\arguments{
  \item{n}{number of rows and columns in the J matrix.}
}
\value{
  an \eqn{n} by \eqn{n} lower triangular matrix with all non-zero entries equal to one. This is an auxiliary function.
}
\description{
Creates a \eqn{n} by \eqn{n} lower triangular matrix with all non-zero entries equal to one.
}
\examples{
J = Jn(n = 3)
J
}

