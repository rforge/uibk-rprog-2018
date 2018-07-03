\name{predict.gevreg}
\alias{predict.gevreg}
\title{Predicted values for fits obtained with gevreg}
\usage{
\method{predict}{gevreg}(object, newdata = NULL,
type = c("response", "location", "scale", "shape", 
"parameter", "probability", "quantile"),
na.action = na.pass, at = 0.5, ...)
}
\arguments{
\item{object}{object of class \emph{gevreg}.}
\item{newdata}{an optional data frame with variables with which to predict.}
\item{type}{type of prediction: \emph{location}, \emph{scale} and \emph{shape} return the location, scale and shape parameter of the predicted GEV distribution. \emph{parameter} returns a data frame with predicted location, scale and shape parameters. \emph{probability} and \emph{quantile} return a data frame of predicted probabilities and quantiles with probabilities and quantiles \emph{at}.}
\item{na.action}{a function which indicates what should happen when the data contains NAs. Default is na.pass}
\item{at}{a value or a vector of values to evaluate the predictive probability (type = "probability"), or quantile probabilities used for type = "quantile". Alternatively, with at = "function" a function is returned that takes at as an argument.}
\item{...}{further arguments passed to or from other methods.}
}

\description{
Provide various types of predictions for gevreg fits.
}

\value{
For types \emph{location}, \emph{scale} or \emph{shape} a vector with either the location or the scale or the shape of the predicted distribution.
For type \emph{parameter} a data frame with all GEV parameters.
For type \emph{probability} or \emph{quantile} a data frame of predicted probabilities or quantiles each column corresponding to an element of \emph{at}.
}

\seealso{\code{\link[gevreg]{gevreg}}}