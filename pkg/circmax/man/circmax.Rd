\name{circmax}
\alias{circmax}

\title{Estimate the parameters of the von Mises distribution by maximum likelihood estimation.}
\usage{
circmax(y, x, z = x, start = NULL, ...)
}
\arguments{
\item{y}{Response vector.}

\item{x}{Design matrix for location.}

\item{z}{Design matrix for scale (density parameter kappa).}

\item{start}{Starting values.}

\item{...}{Additional parameters for optim.}
}
\description{
Estimate the parameters of the von Mises distribution by maximum likelihood estimation.
}

\examples{
library(circular)
library(stats)

# small helper function
plot_circ <- function(x) plot(as.circular(x), stack = TRUE, col = colorspace::rainbow_hcl(length(x)))

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
plot_circ(d)

## Fitting with lm.circular
y <- d$y
x <- d$x1
(lm <- lm.circular(y = y, x = x, init = c(0), type='c-l', verbose=TRUE))

## Maximum Likelihood estimation (inside package)
x <- cbind(1, d$x1)
z <- cbind(1, d$x2)
y <- d$y

(m1 <- circmax(y = y, x = x, z = z))
}

