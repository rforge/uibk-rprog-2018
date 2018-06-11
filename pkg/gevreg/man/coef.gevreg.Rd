\name{coef.gevreg}
\alias{coef.gevreg}
\alias{logLik.gevreg}
\alias{vcov.gevreg}
\alias{terms.gevreg}
\alias{model.frame.gevreg}
\alias{model.matrix.gevreg}
\alias{residuals.gevreg}


\title{Methods for gevreg objetcs}

\description{Various methods for extracting information from fitted \code{gevreg} objects.}

\usage{
\method{coef}{gevreg}(object, model = c("full", "location", "scale", "shape"), ...)
\method{vcov}{gevreg}(object, model = c("full", "location", "scale", "shape"), ...)
\method{terms}{gevreg}(x, model = c("location", "scale", "shape", "full"), ...)
\method{model.frame}{gevreg}(formula, ...)
\method{model.matrix}{gevreg}(object, model = c("location", "scale", "shape"), ...)
\method{residuals}{gevreg}(object, type = c("standardized", "pearson", "response"), ...)
}
\arguments{
\item{object, x}{an object of class \code{gevreg}.}
\item{model}{model for which coefficients shall be returned.}
\item{type}{type of residuals that shall be computed.}
\item{formula}{a model link{formula} formula or \link{terms} object or an R object.}
\item{...}{further arguments passed to or from other methods.}
}

\details{
In addition to the methods above, a set of standard extractor functions for "gevreg" objects is available, including methods to the generic functions \link{logLik} and \link{print}.} 


\seealso{
\link{gevreg} \link[gevreg]{print.gevreg}
}


