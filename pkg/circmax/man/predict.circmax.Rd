\name{predict.circmax}
\alias{predict.circmax}

\title{Predicted/Fitted Values for CIRCMAX Fits}

\description{
Obtains various types of predictions for \code{circmax} models.
}

\usage{
\method{predict}{circmax}(object, newdata = NULL, type = c("location", "scale", 
  "parameter", "probability", "quantile"), 
  na.action = na.pass, at = 0.5, \dots)
}

\arguments{
  \item{object}{an object of class \code{"circmax"}.}
  \item{newdata}{an optional data frame in which to look for variables with
  which to predict.}
  \item{type}{type of prediction: \code{"location"} returns the location of 
  the predicted distribution.
  \code{"scale"} returns the scale of the predicted distribution.
  \code{"parameter"} returns a data frame with predicted location and scale 
  parameters.
  \code{"probability"} evaluates the predictive CDF at \code{at}.
  \code{"quantile"} returns a matrix of predicted quantiles with quantile
  probabilities \code{at}.}
  \item{na.action}{a function which indicates what should happen when the data
  contain \code{NA}s. Default is na.pass}
  \item{at}{a vector of values to evaluate the predictive probability (\code{type = "probability"}), 
  or a vector of quantile probabilities used for \code{type =  "quantile"}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
For type \code{"location"}, or \code{"scale"} a vector with
either the location or the scale of the predicted distribution.
}

\seealso{\code{\link{circmax}}}

\keyword{regression}
