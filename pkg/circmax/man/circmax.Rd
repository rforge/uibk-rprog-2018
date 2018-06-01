\name{circmax}
\alias{circmax}
\alias{circmax_fit}
\alias{circmax_control}

\title{Maximum Likelihood Estimation of a Circular Response}

\description{
  Estimate the parameters of the von Mises distribution by maximum likelihood estimation.
}
\usage{
circmax(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE,
  control = circmax_control(\dots), \dots)

circmax_fit(x, y, z = NULL, control)

circmax_control(maxit = 5000, start = NULL, ...)
}
\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x | z} where
    \code{y} is the response and \code{x} and \code{z} are regressor variables
    for the location and the scale (density) of the von Mises distribution.}
  \item{data}{an optional data frame containing the variables occurring in the
    formulas.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y}{for \code{circmax}: logical. If \code{TRUE} the model matrix and
    response vector used for fitting are returned as components of the returned value.
    For \code{circmax.fit}: \code{x} is a design matrix with regressors for the
    location and \code{y} is a vector of observations.}
  \item{z}{a design matrix with regressors for the scale.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
  \item{control, maxit, start}{a list of control parameters passed to \code{\link{optim}} .}
}

\details{
\code{circmax} estimates the parameters of the the von Mises distribution by maximum likelihood estimation.

\code{circmax_fit} is the lower level function where the actual
fitting takes place.
}

\value{
  An object of class \code{"circmax"}.
}

\examples{
library(circular)
library(stats)

# small helper function
plot_circ <- function(x){
  plot(as.circular(x), stack = TRUE, col = colorspace::rainbow_hcl(length(x)))
}

## data-generating process
dgp <- function(n = 100, coef = c(2, 5, 5, 0)) {
  d <- data.frame(
    x1 = rnorm(n, 0, 0.2),
    x2 = rnorm(n, 0, 0.2)
  )
  mean <- coef[1] + 2 * atan(coef[2] * d$x1)
  kappa <- exp(coef[3] + coef[4] * d$x2)
  d$y <- NULL
  for(i in 1:n) {
    d[i, "y"] <- circular::rvonmises(1, mu = circular(mean[i]), kappa = kappa[i])
  }
  return(d)
}

set.seed(2018-05-29)
d <- dgp()
if(require("colorspace")) {
  plot_circ(d)
}

## Fitting with lm.circular
y <- d$y
x <- d$x1
(lm <- lm.circular(y = y, x = x, init = c(0), type='c-l', verbose=TRUE))

## Maximum Likelihood estimation (inside package)
(m1 <- circmax(y ~ x1 | x2, data = d))
}
\keyword{circular, regression}
