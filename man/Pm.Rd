\name{Pm}
\alias{Pm}
\title{Projection Matrix P}
\usage{
Pm(m = 2, nu = 0)
}
\arguments{
  \item{nu}{the degree of the polinomial fit.}
  \item{m}{a positive integer satisfying \eqn{m >= nu} indicating the size of the window for the polinomial fit.}

}
\value{
  an \eqn{m+1} by \eqn{m+1} matrix.
}
\description{
Creates the \eqn{m+1} by \eqn{m+1} projection matrix defined by \eqn{P = D(D'D)^{-1}D'} where \eqn{D} is the design matrix associated to a polynomial regression of degree nu + 1.
}

\details{
To perform matrix inversion, the code makes use of the routine DGETRI in LAPACK, which applies an LU decomposition approach to obtain the inverse matrix. See the LAPACK documentation available at \url{http://www.netlib.org/lapack}.
}
\examples{
P = Pm(m = 5, nu = 0)
P

n = 10
t = 1:n
D = cbind(rep(1,n),t,t^2)

# Calculating in R
PR = D\%*\%solve(t(D)\%*\%D)\%*\%t(D)
# Using the provided function
P = Pm(m = n-1, nu = 1)

# Difference:
sum(abs(P-PR))

}
\author{Taiane Schaedler Prass}
