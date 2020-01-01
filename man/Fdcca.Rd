\name{Fdcca}
\alias{Fdcca}
\title{Detrended Cross-covariance}
\usage{
Fdcca(y1, y2, m = 3, nu = 0, overlap = TRUE)
}
\arguments{
  \item{y1, y2}{vectors corresponding to the time series data. If \eqn{length(y1)} and \eqn{length(y2)} differ, the longer time series is coerced to match the lenght of the shorter.}
  \item{m}{an integer or integer valued vector indicating the size (or sizes) of the window for the polinomial fit. \eqn{min(m)} must be greater or equal than \eqn{nu} or else it will return an error.}
  \item{nu}{a non-negative integer denoting the degree of the polinomial fit applied on the integrated series.}
  \item{overlap}{logical: if true (the default), uses overlapping windows. Otherwise, non-overlapping boxes are applied.}
}
\value{
  A vector of size \eqn{length(m)} containing the detrended cross-covariance considering windows of size \eqn{m+1}, for each \eqn{m} supplied.
}
\description{
Calculates the detrended cross-covariance between two time series \eqn{y1} and \eqn{y2}.
}
\examples{
# Simple usage
y1 = rnorm(100)
y2 = rnorm(100)
F.dcca = Fdcca(y1, y2, m = 3, nu = 0, overlap = TRUE)
F.dcca

# A simple example where y1 and y2 are independent.

ms = 3:50
F.dcca1 = Fdcca(y1, y2, m = ms, nu = 0, overlap = TRUE)
F.dcca2 = Fdcca(y1, y2, m = ms, nu = 0, overlap = FALSE)

plot(ms, F.dcca1, type = "o", xlab = "m", col = "blue",
     ylim = c(min(F.dcca1,F.dcca2),max(F.dcca1,F.dcca2)),
     ylab = expression(F[DCCA]))
points(ms, F.dcca2, type = "o", col = "red")
legend("bottomright", legend = c("overlapping","non-overlapping"),
       col = c("blue", "red"), lty= 1, bty = "n", pch=1)


# A more elaborated example where y1 and y2 display cross-correlation for non-null lags.
# This example also showcases why overlapping windows are usually advantageous.
# The data generating process is the following:
# y1 is i.i.d. Gaussian while y2 is an MA(2) generated from y1.

n = 500
ms = 3:50
theta = c(0.4,0.5)

# Calculating the expected value of the DCCA in this scenario
m_max = max(ms)
vtheta = c(1,theta, rep(0, m_max - length(theta)))
G12 = matrix(0, ncol = m_max+1, nrow = m_max+1)
for(t in 1:(m_max+1)){
  for(h in 0:(m_max+1-t)){
    G12[t,t+h] = vtheta[h+1]
  }
}

EF.dcca = EFdcca(m = ms, nu = 0, G = G12)

# generating the series and calculating the DCCA
burn.in = 100
eps = rnorm(burn.in)

y1 = rnorm(n)
y2 = arima.sim(model = list(ma = theta), n, n.start = burn.in, innov = y1, start.innov = eps)

ms = 3:50
OF.dcca = Fdcca(y1, y2, m = ms, nu = 0, overlap = TRUE)
NOF.dcca = Fdcca(y1, y2, m = ms, nu = 0, overlap = FALSE)

plot(ms, OF.dcca, type = "o", xlab = "m", col = "blue",
     ylim = c(min(NOF.dcca,OF.dcca,EF.dcca),max(NOF.dcca,OF.dcca,EF.dcca)),
     ylab = expression(F[DCCA]))
points(ms, NOF.dcca, type = "o", col = "darkgreen")
points(ms, EF.dcca, type = "o", col = "red")
legend("bottomright", legend = c("overlapping","non-overlapping","expected"),
       col = c("blue", "darkgreen","red"), lty= 1, bty = "n", pch=1)

}

\references{Prass, T.S. and Pumi, G. (2019). On the behavior of the DFA and DCCA in trend-stationary
processes <arXiv:1910.10589>.}
\author{Taiane Schaedler Prass}
