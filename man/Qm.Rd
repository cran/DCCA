\name{Qm}
\alias{Qm}
\title{Projection Matrix Q}
\usage{
Qm(m = 2, nu = 0, P = NULL)
}
\arguments{
  \item{nu}{the degree of the polinomial fit.}
  \item{m}{a positive integer satisfying \eqn{m >= nu} indicating the size of the window for the polinomial fit.}
  \item{P}{optional: the projection matrix such that \eqn{Q = I - P} (see function \code{\link{Pm}}). If this matrix is provided \eqn{m} and \eqn{nu} are ignored.}
}
\value{
  an \eqn{m+1} by \eqn{m+1} matrix.
}
\description{
Creates the \eqn{m+1} by \eqn{m+1} projection matrix defined by \eqn{Q = I-P} where \eqn{I} is the the \eqn{m+1} by \eqn{m+1} identity matrix and \eqn{P} is the \eqn{m+1} by \eqn{m+1} projection matrix into the space generated by polynomials of degree \eqn{nu + 1}.
}
\seealso{\code{\link{Pm}} which generates the projection matrix \eqn{P}.}
\examples{
Q = Qm(m = 3, nu = 0)
Q
# same as
P = Pm(m = 3, nu = 0)
Q = Qm(P = P)
Q
}
