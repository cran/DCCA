\name{covF2dfa}
\alias{covF2dfa}
\title{Autocovariance function of the detrended variance}
\usage{
covF2dfa(m = 3, nu = 0, h = 0, overlap = TRUE, G, Cumulants = NULL)
}
\arguments{
  \item{m}{an integer or integer valued vector indicating the size of the window for the polinomial fit. \eqn{min(m)} must be greater or equal than \eqn{nu} or else it will return an error.}
  \item{nu}{a non-negative integer denoting the degree of the polinomial fit applied on the integrated series.}
  \item{h}{an integer or integer valued vector indicating the lags for which the autocovariance function is to be calculated.}
  \item{overlap}{logical: if true (the default), overlapping boxes are used for calculations. Otherwise, non-overlapping boxes are applied.}
  \item{G}{the autocovariance matrix for the original time series. The dimension of \eqn{G} must be \eqn{(max(m)+max(h)+1)} by \eqn{(max(m)+max(h)+1)} if overlap = TRUE and \eqn{(max(m)+max(h))(max(h)+1)} by \eqn{(max(m)+max(h))(max(h)+1)} otherwise.}
  \item{Cumulants}{The matrix containing the joint cumulants for lags. Dimension must be \eqn{(max(m)+1)*nrow(G)}. If not provided, it is assumed that the cumulants are all zero.}
}
\value{
  A matrix with the autocovariance of lag \eqn{h}, for each value of \eqn{m} provided. This matrix is obtained from expressions (21) for \eqn{h = 0} and (22) for \eqn{h > 0} in Prass and Pumi (2019).
}
\description{
Calculates the autocovariance of the detrended variance.
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

covdfa1 = covF2dfa(m = ms, nu = 0, h = hs,
                   overlap = TRUE, G = Gamma1, Cumulants = NULL)

covdfa2 = covF2dfa(m = ms, nu = 0, h = hs,
                  overlap = TRUE, G = Gamma2, Cumulants = NULL)

cr = rainbow(100)
plot(ms, covdfa1[,1], type = "l", ylim = c(0,20),
    xlab = "m", ylab = expression(gamma[DFA](h)), col = cr[1])
for(i in 2:ncol(covdfa1)){
  points(ms, covdfa1[,i],  type = "l", col = cr[i])
}

lattice::wireframe(covdfa1, drape = TRUE,
    col.regions = rev(rainbow(150))[50:150],
    zlab = expression(gamma[DFA]), xlab = "m", ylab = "h")
}
}

\references{Prass, T.S. and Pumi, G. (2019). On the behavior of the DFA and DCCA in trend-stationary
processes <arXiv:1910.10589>.}
\author{Taiane Schaedler Prass}
