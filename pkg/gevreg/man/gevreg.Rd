\name{gevreg}
\alias{gevreg}
\alias{gevreg_fit}
\alias{gevreg_control}

\title{Fit a GEV distribution}

\description{
  Estimate the parameters of the GEV distribution using maximum likelihood.
}
\usage{
gevreg(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE, z = FALSE, v = FALSE,
  control = gevreg_control(\dots), \dots)

gevreg_fit(x, y, z = NULL, v = NULL, control)

gevreg_control(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, ...)
}

\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x | z | v}  where
    \code{y} is the response and \code{x}, \code{z} and \code{v} are regressor variables
    for the location, scale and shape parameters of the Generalized Extreme value Distribution (GEV). 
	For details how to set up the formula see details and \code{\link{Formula}}}.
  \item{data}{a data frame containing the covariables used to fit the respons(es). Variables in the formula 
  must occur as column name in data.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y, z, v}{for \code{gevreg}: logical. If \code{TRUE} the model matrix and
    response matrix used for fitting are returned as part of the returned gevreg object.
    For \code{gevreg.fit}: \code{x, z, v} are design matrices with regressors for the
    location, scale and shape parameters of the GEV and \code{y} is a vector of observed parameters.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
  \item{control, maxit, start}{a list of control parameters passed to \code{\link{optim}} .}
  
  \item{grad}{should a gradient-based method used for optimization.}
  \item{hessian}{should a hessian be computed during the optimaztion process. Possible values are \code{"none"}, \code{"optim"} or \code{"numderic"}, which is the default.}
}

\details{
\code{gevreg}  Fit a GEV distribution simultaneously to a number of observation sites. 
The spatial component is achieved by maximisng the sum of the single log-likelihoods at the stations. 
The response y must be either a two or a three columned matrix. If it contains two columns, 
the first column is supposed to be the location and the second the scale parameter of the GEV. In that case, 
the  shape parameter is considered zero, in which case the Gumbel distribution as a special case of the GEV is fitted. 
If y contains three columns they are considered as location, scale and shape parameter, respectively. 
If 1 is given in either of the formulars right-hand-side, the respective parameter is taken as constant. 

\code{gevreg_fit} is the lower level function which provides the actual maximum likelihood fitting.
}

\value{
  An object of class \code{gevreg} which inherits from \code{optim} with components: 
  \item{coefficients}{Either all, or if specified e.g. as location, scale or shape the coefficients of the fitted model for the three GEV parameters are returned. }
  \item{loglik}{The maximised log-likelihood value.}
  
  
}

\references{
 	Blanchet J, Lehning M (2010). Mapping snow depth return levels: smooth spatial modeling versus station interpolation. 
	\emph{Hydrol. Earth Syst. Sci.}, \bold{14}, 2527--2544.
	\url{https://www.hydrol-earth-syst-sci.net/14/2527/2010/hess-14-2527-2010.pdf}.
}


\examples{
data("gevregdata")
d <- gevregdata

## Fit GEV to the snow depths of the first station in the dataset
d1 <- d[which(d$station == 1),]
gevreg(formula = hs ~ 1, data = d1)

## this is equivalent
gevreg(formula = hs ~ 1, data = d, subset = station == 1)

## print output in standard or extended format
m <- gevreg(formula = hs ~ 1, data = d, subset = station == 2)
print(m)
print(m, ext = TRUE)

## the hessian is computed by default 
m <- gevreg(hs ~ 1, d, station == 3)
m$hessian

## but you could also use the numderiv package after the optimzation
mnum <- gevreg(hs ~ 1, d, station == 3, hessian = "numderiv")
mnum$hessian

}

\keyword{gev, regression}
