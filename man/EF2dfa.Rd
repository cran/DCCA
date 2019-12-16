\name{EF2dfa}
\alias{EF2dfa}
\title{Expected value of the detrended variance}
\usage{
EF2dfa(m = 3, nu = 0, G, K = NULL)
}
\arguments{
  \item{m}{an integer or integer valued vector indicating the size of the window for the polinomial fit. \eqn{min(m)} must be greater or equal than \eqn{nu} or else it will return an error.}
  \item{nu}{a non-negative integer denoting the degree of the polinomial fit applied on the integrated series.}
  \item{G}{the autocovariance matrix for the original time series. The dimension of \eqn{G} must be \eqn{(max(m)+1)} by \eqn{(max(m)+1)}.}
  \item{K}{optional: the matrix \eqn{K}. If this matrix is provided and \eqn{m} is an integer, then \eqn{nu} is ignored.}
}
\value{
A vector of size \eqn{length(m)} containing the expected values of the detrended variance corresponding to the values of \eqn{m} provided. This is expression (20) in Prass and Pumi (2019).
}
\description{
Calculates the expected value of the detrended variance.
}
\examples{
m = 3
K = Km(m = m, nu = 0)
G = diag(m+1)
EF2dfa(G = G, K = K)
# same as
EF2dfa(m = 3, nu = 0, G = G)

# An AR(1) example
phi = 0.4
n = 500
burn.in = 50
eps = rnorm(n + burn.in)
z.temp = numeric(n + burn.in)
z.temp[1] = eps[1]
for(i in 2:(n + burn.in)){
  z.temp[i] = phi*z.temp[i-1] + eps[i]
}
z = z.temp[(burn.in + 1):(n + burn.in)]

F2.dfa = F2dfa(z, m = 3:100, nu = 0, overlap = TRUE)
plot(3:100, F2.dfa, type="o", xlab = "m")

}
\references{Prass, T.S. and Pumi, G. (2019). On the behavior of the DFA and DCCA in trend-stationary
processes <arXiv:1910.10589>.}
\author{Taiane Schaedler Prass}
