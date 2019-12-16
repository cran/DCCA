\name{EFdcca}
\alias{EFdcca}
\title{Expected value of the detrended cross-covariance}
\usage{
EFdcca(m = 3, nu = 0, G, K = NULL)
}
\arguments{
  \item{m}{an integer or integer valued vector indicating the size of the window for the polinomial fit. \eqn{min(m)} must be greater or equal than \eqn{nu} or else it will result in an error.}
\item{nu}{a non-negative integer denoting the degree of the polinomial fit applied on the integrated series.}
  \item{G}{the cross-covariance matrix for the original time series. The dimension of \eqn{G} must be \eqn{(max(m)+1)} by \eqn{(max(m)+1)}.}
  \item{K}{optional: the matrix \eqn{K}. If this matrix and \eqn{m} are provided, then \eqn{nu} is ignored.}
}
\value{
  a size \eqn{length(m)} vector containing the expected values of the detrended cross-covariance corresponding to the values of \eqn{m} provided. This is expression (23) in Prass and Pumi (2019).
}
\description{
Calculates the expected value of the detrended cross-covariance given a cross-covariance matrix.
}
\examples{
m = 3
K = Km(m = m, nu = 0)
G = diag(m+1)
EFdcca(G = G, K = K)
# same as
EFdcca(m = 3, nu = 0, G = G)


}
\references{Prass, T.S. and Pumi, G. (2019). On the behavior of the DFA and DCCA in trend-stationary
processes <arXiv:1910.10589>.}
\author{Taiane Schaedler Prass}
