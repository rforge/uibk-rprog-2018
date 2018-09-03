\name{dgevr}
\alias{dgevr}
\alias{pgevr}
\alias{qgevr}


\title{GEV density, distribution and quantile function}

\description{Vectorized density, distribution and quantile function for the GEV distribution.}

\usage{
dgevr(x, loc = 0, scale = 1, shape = 0, log = FALSE)
pgevr(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
qgevr(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
}

\arguments{
 \item{x, q}{vector of quantiles for \code{dgevr} and \code{pgevr}.}
  \item{p}{vector of probabilitis for \code{qgevr}.}
  \item{loc, scale, shape}{vector of location, scale and shape parameters.}
  \item{log}{locigal; if \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are \eqn{Pr[X <= x]}, otherwise, \eqn{Pr[X > x]}.}
}

\details{
\code{dgevr,pgevr,qgevr}{ GEV density, distribution and quantile function which can take location, scale and shape parameters as vector. Simply modifed versions of \code{dgev}, \code{pgev} and \code{qgev} of package \code{\link{SpatialExtremes}}.}} 

\value{
If not specified, "loc", "scale" and "shape" assume the default values of "0", "1", and "0", respectively.
}

\seealso{
\link{SpatialExtremes}
}

\examples{
dgevr(seq(10,100,10))
dgevr(rep(seq(10,100,10),3),loc=1:3,scale=1:3,shape=rep(0,3))
}


