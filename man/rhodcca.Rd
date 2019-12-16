\name{rhodcca}
\alias{rhodcca}
\title{Detrended Cross-correlation coefficient}
\usage{
rhodcca(y1, y2, m = 3, nu = 0, overlap = TRUE)
}
\arguments{
  \item{y1, y2}{vectors corresponding to the time series data. If \eqn{length(y1)} and \eqn{length(y2)} differ, the longer time series is coerced to match the lenght of the shorter.}
  \item{m}{an integer value or a vector of integer values indicating the size of the window for the polinomial fit. \eqn{min(m)} must be greater or equal than \eqn{nu} or else it will return an error.}
  \item{nu}{the degree of the polynomial fit}
  \item{overlap}{logical: if true (the default), uses overlapping windows. Otherwise, non-overlapping boxes are applied.}
}
\value{A list containing the following elements, calculated considering windows of size \eqn{m+1}, for each \eqn{m} supplied:
\item{F2dfa1, F2dfa2}{The detrended variances for \eqn{y1} and \eqn{y2}, respectively.}
\item{Fdcca}{The detrended cross-covariance.}
\item{rhodcca}{The detrended cross-correlation coefficient.}
}
\description{
Calculates the detrended cross-correlation coefficient for two time series \eqn{y1} and \eqn{y2}.
}
\note{The time series \eqn{y1} and \eqn{y2} must have the same sample size.}
\examples{
y1 = rnorm(100)
y2 = rnorm(100)
rho.dccam1 = rhodcca(y1, y2, m = 3, nu = 0, overlap = TRUE)
rho.dccam1

rho.dccam2 = rhodcca(y1, y2, m = c(3,6,8), nu = 0, overlap = TRUE)
rho.dccam2
}
\references{Prass, T.S. and Pumi, G. (2019). On the behavior of the DFA and DCCA in trend-stationary
processes <arXiv:1910.10589>.}
\author{Taiane Schaedler Prass}
\seealso{\code{\link{F2dfa}} which calculated the DFA and \code{\link{Fdcca}} which calculated the DCCA of two given time series.}
