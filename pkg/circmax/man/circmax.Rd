\name{circmax}
\alias{circmax}
\alias{circmax_fit}
\alias{circmax_control}

\title{Circular Regression with Maximum Likelihood Estimation}

\description{
  Fitting a regression models for a circular response by maximum likelihood estimation.
}
\usage{
circmax(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE,
  control = circmax_control(\dots), which.method = "default", \dots)

circmax_fit(x, y, z = NULL, control, which.method = "default")

circmax_control(maxit = 5000, start = NULL, ...)
}
\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x | z} where
    \code{y} is the response and \code{x} and \code{z} are regressor variables
    for the location and the scale (density) of the von Mises distribution.}
  \item{data}{an optional data frame containing the variables occurring in the
    formulas; y has to be given in radians.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y}{for \code{circmax}: logical. If \code{TRUE} the model matrix and
    response vector used for fitting are returned as components of the returned value.
    For \code{circmax_fit}: \code{x} is a design matrix with regressors for the
    location and \code{y} is a vector of observations given in radians.}
  \item{z}{a design matrix with regressors for the scale.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
  \item{control, maxit, start}{a list of control parameters passed to \code{\link{optim}} .}
  \item{which.method}{Which nll is tested.}
}

\details{
\code{circmax} fits a regression model for a circular response assuming a von Mises distribution.


\code{circmax_fit} is the lower level function where the parameters of the the von Mises distribution 
are fitted by maximum likelihood estimation.
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

## Example 1: Simulated Data:

## data-generating process
dgp <- function(n = 1000, beta = c(3, 5, 2), gamma = c(3, 3)) {
  m <- length(beta) - 1  # here: number of betas minus intercept
  p <- length(gamma) -1  # here: number of gammas minus intercept

  #d <- sapply(1:(m + p), function(x) rnorm(n, runif(1, 0, 2 * pi), 0.2))
  #d <- sapply(1:(m + p), function(x) runif(n, 0, 1))
  d <- sapply(1:(m + p), function(x) rnorm(n, 0, 0.2))
  colnames(d) <- paste0("x", 1:(m + p))

  mu <- 2 * atan(beta[1]) + 2 * atan(crossprod(t(d[, 1:m, drop = FALSE]), beta[-1]))
  kappa <- exp(gamma[1] + crossprod(t(d[, m + 1:p, drop = FALSE]), gamma[-1]))
  d <- data.frame(d)
  d$y <- NULL
  for(i in 1:n) {
    d[i, "y"] <- circular::rvonmises(1, mu = circular::circular(mu[i]), kappa = kappa[i])
  }
  return(d)
}

set.seed(2018-05-29)
d <- dgp()
if(require("colorspace")) {
  plot(d[, colnames(d)[(colnames(d) != c("y"))]])
  plot_circ(d$y)
}

## Fitting with lm.circular
#y <- d$y
#x <- d$x1
#(lm <- lm.circular(y = y, x = x, init = c(0), type='c-l', verbose=TRUE))

## Maximum Likelihood estimation (inside package)
#(m1 <- circmax(y ~ x1 + x2 | x3, data = d, which.method = "version1", 
#  control = circmax_control(trace = 3)))
#(predict(m1))  ## predict currently just for default working
(m2 <- circmax(y ~ x1 + x2 | x3, data = d, which.method = "default", 
  control = circmax_control(trace = 6)))
#(m3 <- circmax(y ~ x1 + x2 | x3, data = d, which.method = "L-BFGS-B", 
#  control = circmax_control(trace = 6)))

## Example 2: Periwinkle Dataset of Fisher and Lee, 1992:
distance <- c(107, 46, 33, 67, 122, 69, 43, 30, 12, 25, 37, 69, 5, 83, 
  68, 38, 21, 1, 71, 60, 71, 71, 57, 53, 38, 70, 7, 48, 7, 21, 27)
directdeg <- c(67, 66, 74, 61, 58, 60, 100, 89, 171, 166, 98, 60, 197, 
  98, 86, 123, 165, 133, 101, 105, 71, 84, 75, 98, 83, 71, 74, 91, 38, 200, 56)
cdirect <- circular(directdeg * 2 * pi/360)
plot(as.numeric(cdirect) ~ distance, ylim = c(0, 4*pi), pch = 20)
points(as.numeric(cdirect) + 2*pi ~ distance, pch = 20)

pm1 <- lm.circular(type = "c-l", y = cdirect, x = distance, init = 0.0)
pm2 <- circmax(cdirect ~ distance, data = data.frame(cbind(distance, cdirect)), 
  which.method = "version1", control = circmax_control(trace = 3))
pm3 <- circmax(cdirect ~ distance, data = data.frame(cbind(distance, cdirect)), 
  which.method = "default", control = circmax_control(trace = 3))
#pm4 <- circmax(cdirect ~ distance, data = data.frame(cbind(distance, cdirect)), 
#  which.method = "L-BFGS-B", control = circmax_control(trace = 3))

ll_fun <- function(eta.mu, eta.kappa, y, log = FALSE) {
  mu <- 2 * atan(eta.mu)
  kappa <- exp(eta.kappa)
  ll <- dvonmises(y, mu = mu, kappa = kappa, log = TRUE)
  #ll <- NULL
  #for(i in 1:length(y)) {
  #  ll <- c(ll,circular::dvonmises(y[i], as.circular(mu), kappa, log = TRUE))
  #}
  if(log == TRUE) {
    return(sum(ll)) 
  } else { 
    return(exp(sum(ll)))
  }
}

par(mfrow = c(1, 2))
eta.mu <- seq(-100, 100, 0.1)
plot(eta.mu, sapply(eta.mu, ll_fun, y = cdirect, eta.kappa = coef(pm2, model = "scale")[1], 
  log = TRUE), xlab = expression(mu), ylab = "logLik", type = "l")

i.max <- which.max(sapply(eta.mu, ll_fun, y = cdirect, eta.kappa = coef(pm2, model = "scale")[1], 
  log = TRUE))
eta.kappa <- seq(0, 200, 0.1)
plot(eta.kappa, sapply(eta.kappa, ll_fun, y = cdirect, eta.mu = eta.mu[i.max], 
  log = TRUE), xlab = expression(kappa), ylab = "logLik", type = "l")


#eta.mu <- seq(-4, 4, 0.1)
#eta.kappa <- seq(0, 2, 0.1)
#tmp.ll <- tmp.mu <- tmp.kappa <- matrix(0, nrow = length(eta.mu), ncol = length(eta.kappa))
#for(i in 1:length(eta.mu)){
#  for(j in 1:length(eta.kappa)){
#    tmp.ll[i,j] <- ll_fun(eta.mu = eta.mu[i], eta.kappa = eta.kappa[j], y = cdirect, log = TRUE)
#    tmp.mu[i,j] <- eta.mu[i]
#    tmp.kappa[i,j] <- eta.kappa[j]
#  }
#}
#
#ll3d <- data.frame(eta.mu = as.numeric(tmp.mu), eta.kappa = as.numeric(tmp.kappa), 
#  ll = as.numeric(tmp.ll))
#bamlss::plot3d(ll3d, ticktype = "detailed", type = "mba", theta = 140, phi = 40, d.persp = 0.5)

}
\keyword{circular, regression}
