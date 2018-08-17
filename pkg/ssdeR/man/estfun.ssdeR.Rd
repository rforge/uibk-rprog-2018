\name{estfun.ssdeR}
\alias{estfun.ssdeR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{estfun.ssdeR}
\description{
\code{estfun} method for \code{ssdeR} objects}
\usage{
estfun.ssdeR(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a fitted model object.}
  \item{\dots}{
arguments passed to methods.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix containing the empirical estimating functions. Typically, this should be an \eqn{n x k} matrix corresponding to n observations and k parameters. The columns should be named as in coef or terms, respectively.

The estimating function (or score function) for a model is the derivative of the objective function with respect to the parameter vector. The empirical estimating functions is the evaluation of the estimating function at the observed data (n observations) and the estimated parameters (of dimension k).

}
\references{
Zeileis A (2006), Object-Oriented Computation of Sandwich Estimators. Journal of Statistical Software, 16(9), 1--16. URL \url{http://www.jstatsoft.org/v16/i09/}.
}

\author{
Michael Brottrager
}

\examples{
## data

## appl


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~estfun }% use one of  RShowDoc("KEYWORDS")
