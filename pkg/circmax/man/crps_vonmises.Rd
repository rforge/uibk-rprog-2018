\name{crps_vonmises}
\alias{crps_vonmises}

\title{Circular CRPS}

\description{
Continuous Ranked Probability Score (CRPS) for a circular response
following the von Mises distribution.
}

\usage{
crps_vonmises(y, mu, kappa, sum = FALSE)
}

\arguments{
  \item{y}{numeric. Circular response.}
  \item{mu}{numeric. Location parameter of the von Mises distribution.}
  \item{kappa}{numeric. Concentration parameter of the von Mises distribution.}
  \item{sum}{logical. Should the sum of the CRPS-values over all response values be returned.}
}

\details{
  }

\value{
For \code{sum=TRUE} the sum of the CRPS-values of all response values is returned. Otherwise
a vector of the same length as \code{y} with the observation-wise CRPS-values is returned.
}

\examples{
require("CircStats")
y <- rvm(100, 0.5, 4)
y[y>pi] <- y[y>pi]-2*pi
testfit <- circfit(y)
testcrps <- crps_vonmises(y, mu = testfit$coefficients[1], kappa = testfit$coefficients[2])
}


\seealso{\code{\link{circtree}}}

\keyword{CPRS, circular}
