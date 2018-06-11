\name{dgev.gevreg}
\alias{dgev.gevreg}
\alias{pgev.gevreg}
\alias{qgev.gevreg}


\title{GEV density, distribution and quantile function}

\description{Vectorized sensity, distribution and quantile function for the GEV distribution.}

\usage{
dgev.gevreg(x, loc = 0, scale = 1, shape = 0, log = FALSE)
pgev.gevreg(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
qgev.gevreg(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
}

\arguments{
 \item{x, q}{vector of quantiles for \code{dgev.gevreg} and \code{pgev.gevreg}.}
  \item{p}{vector of probabilitis for \code{qgev.gevreg}.}
  \item{loc, scale, shape}{vector of location, scale and shape parameters.}
  \item{log}{locigal; if \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are Pr[ X <= x], otherwise, Pr[X > x].}
}

\details{
\code{dgev.gevreg,pgev.gevreg,qgev.gevreg}{ GEV density, distribution and quantile function which can take location, scale and shape parameters as vector. Simply modifed versions of \code{\link[SpatialExtremes]{dgev}}, \code{\link[SpatialExtremes]{pgev}} and \code{\link[SpatialExtremes]{qgev}}.}} 

\value{
If not specified, "loc", "scale" and "shape" assume the default values of "0", "1", and "0", respectively.
}

\seealso{
\link[SpatialExtremes]{dgev}
}

\examples{
dgev.gevreg(seq(10,100,10))
dgev.gevreg(rep(seq(10,100,10),3),loc=1:3,scale=1:3,shape=rep(0,3))
}


