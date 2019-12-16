\name{F2dfa}
\alias{F2dfa}
\title{Detrended Variance}
\usage{
F2dfa(y, m = 3, nu = 0, overlap = TRUE)
}
\arguments{
  \item{y}{vector corresponding to the time series data.}
  \item{m}{an integer or integer valued vector indicating the size (or sizes) of the window for the polinomial fit. \eqn{min(m)} must be greater or equal than \eqn{nu} or else it will return an error.}
\item{nu}{a non-negative integer denoting the degree of the polinomial fit applied on the integrated series.}
  \item{overlap}{logical: if true (the default), uses overlapping windows. Otherwise, non-overlapping boxes are applied.}
}
\value{
A vector of size \eqn{length(m)} containing the detrended variance considering windows of size \eqn{m+1}, for each \eqn{m} supplied.
}
\description{
Calculates the detrended variance based on a given time series.
}
\examples{
# Simple usage
y = rnorm(100)
F2.dfa = F2dfa(y, m = 3, nu = 0, overlap = TRUE)
F2.dfa

vF2.dfa = F2dfa(y, m = 3:5, nu = 0, overlap = TRUE)
vF2.dfa


###################################################
# AR(1) example showing how the DFA varies with phi

phi = (1:8)/10
n = 300
z = matrix(nrow = n, ncol = length(phi))
for(i in 1:length(phi)){
  z[,i] = arima.sim(model = list(ar = phi[i]), n)
}

ms = 3:50
F2.dfa = matrix(ncol = length(phi), nrow = length(ms))

for(j in 1:length(phi)){
  F2.dfa[,j] = F2dfa(z[,j], m = ms  , nu = 0, overlap = TRUE)
}

cr = rainbow(length(phi))
plot(ms, F2.dfa[,1], type = "o", xlab = "m", col = cr[1],
    ylim = c(0,max(F2.dfa)), ylab = "F2.dfa")
for(j in 2:length(phi)){
  points(ms, F2.dfa[,j], type = "o", col = cr[j])
}
legend("topleft", lty = 1, legend = phi, col = cr, bty = "n", title = expression(phi), pch=1)


##############################################################################
# An MA(2) example showcasing why overlapping windows are usually advantageous
n = 300
ms = 3:50
theta = c(0.4,0.5)

# Calculating the expected value of the DFA in this scenario
m_max = max(ms)
vtheta = c(c(1,theta, rep(0, m_max - length(theta))))
G = matrix(0, ncol = m_max+1, nrow = m_max+1)
for(t in 1:(m_max+1)){
  for(h in 0:(m_max+1-t)){
    G[t,t+h] = sum(vtheta[1:(length(vtheta)-h)]*vtheta[(1+h):length(vtheta)])
    G[t+h,t] = G[t,t+h]
  }
}

EF2.dfa = EF2dfa(m = ms, nu = 0, G = G)

z = arima.sim(model = list(ma = theta), n)

ms = 3:50
OF2.dfa = F2dfa(z, m = ms, nu = 0, overlap = TRUE)
NOF2.dfa = F2dfa(z, m = ms, nu = 0, overlap = FALSE)

plot(ms, OF2.dfa, type = "o", xlab = "m", col = "blue",
    ylim = c(0,max(OF2.dfa,NOF2.dfa,EF2.dfa)), ylab = "F2.dfa")
points(ms, NOF2.dfa, type = "o", col = "darkgreen")
points(ms, EF2.dfa, type = "o", col = "red")
legend("bottomright", legend = c("overlapping","non-overlapping","expected"),
            col = c("blue", "darkgreen","red"), lty= 1, bty = "n", pch=1)

}

\references{Prass, T.S. and Pumi, G. (2019). On the behavior of the DFA and DCCA in trend-stationary
processes <arXiv:1910.10589>.}
\author{Taiane Schaedler Prass}
