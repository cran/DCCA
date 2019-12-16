\name{covFdcca}
\alias{covFdcca}
\title{Autocovariance function of the detrended cross-covariance}
\usage{
covFdcca(m = 3, nu = 0, h = 0, overlap = TRUE, G1, G2, G12, Cumulants = NULL)
}
\arguments{
  \item{m}{an integer or integer valued vector indicating the size of the window for the polinomial fit. \eqn{min(m)} must be greater or equal than \eqn{nu} or else it will return an error.}
    \item{nu}{a non-negative integer denoting the degree of the polinomial fit applied on the integrated series.}
  \item{h}{an integer or integer valued vector indicating the lags for which the autocovariance function is to be calculated. Negative values are not allowed.}
  \item{overlap}{logical: if true (the default), overlapping boxes are used for calculations. Otherwise, non-overlapping boxes are applied.}
  \item{G1, G2}{the autocovariance matrices for the original time series. The dimension of \eqn{G1} and \eqn{G2} must be compatible with the highest values in  vectors \eqn{m} and \eqn{h}. More specifically, the dimension of \eqn{G1} and \eqn{G2} is \eqn{(max(m)+max(h)+1)} by \eqn{(max(m)+max(h)+1)} if overlap = TRUE and \eqn{dim(G1) = dim(G2) = (max(m)+max(h))(max(h)+1)} by \eqn{(max(m)+max(h))(max(h)+1)} otherwise.  }
  \item{G12}{the cross-covariance matrix for the original time series. The dimension of \eqn{G12}  must be compatible with the highest values in  vectors \eqn{m} and \eqn{h}. If overlap = TRUE, \eqn{dim(G12) = [(max(m)+1)*(max(h)+1) - max(m)*max(h)]} by \eqn{[(max(m)+1)*(max(h)+1) - max(m)*max(h)]} and \eqn{dim(G12) = [(max(m)+1)*(max(h)+1)]} by \eqn{[max(m)+1)*(max(h)+1)]}, otherwise}
  \item{Cumulants}{The matrix of cumulants. If not provided, it is assumed that the cumulants are all zero.}
}
\value{
  A matrix of dimension \eqn{lenght(h)} by \eqn{length(m)} with the autocovariance of lag \eqn{h} (rows), for each value of \eqn{m} (columns) provided. This matrix is obtained from expressions (24) for \eqn{h = 0} and (25) for \eqn{h > 0} in Prass and Pumi (2019).
}
\description{
Calculates the autocovariance of the detrended cross-covariance.
}
\examples{
\dontrun{
ms = seq(3,100,1)
hs = seq(0,50,1)
overlap = TRUE
nu = 0
m_max = (max(ms)+1)*(max(hs)+1) - max(ms)*max(hs)*as.integer(overlap)

theta = c(c(1,(20:1)/10), rep(0, m_max - 20))
Gamma1 = diag(m_max+1)
Gamma2 = matrix(0, ncol = m_max+1, nrow = m_max+1)
Gamma12 = matrix(0, ncol = m_max+1, nrow = m_max+1)
for(t in 1:(m_max+1)){
    for(h in 0:(m_max+1-t)){
        Gamma2[t,t+h] = sum(theta[1:(length(theta)-h)]*theta[(1+h):length(theta)])
        Gamma2[t+h,t] = Gamma2[t,t+h]
        Gamma12[t,t+h] = theta[h+1]
    }
}

covdcca = covFdcca(m = ms, nu = 0, h = hs,
                   G1 = Gamma1, G2 = Gamma2, G12 = Gamma12)
                   }
}

\references{Prass, T.S. and Pumi, G. (2019). On the behavior of the DFA and DCCA in trend-stationary
processes <arXiv:1910.10589>.}
\author{Taiane Schaedler Prass}
