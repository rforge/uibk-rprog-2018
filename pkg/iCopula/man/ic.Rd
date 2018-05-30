\name{ic_fit}
\alias{ic_fit}
\title{
Fitter for Implicit Copula Regression
}
\description{
This is a the fitter function for implicit copula regression.
}
\usage{
ic_fit(x, y, control = ic_control(), ...)
}
\arguments{
  \item{x}{A vector with a covariate.}
  \item{y}{A response vector.}
  \item{control}{A list with control parameter.}
  \item{...}{Not used yet.}
}
\value{A list as returned by \code{optim}.}
\examples{
data(mcycle, package = "MASS")
fit <- ic_fit(mcycle$times, mcycle$accel)
print(fit)
}
\keyword{models}
\keyword{regression}
