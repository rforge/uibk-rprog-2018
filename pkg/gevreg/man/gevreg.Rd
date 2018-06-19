\name{gevreg}
\alias{gevreg}
\alias{gevreg_fit}
\alias{gevreg_control}


\title{Maximum Likelihood GEV Fitting}

\description{
  Fit a GEV distribution to observations using maximum likelihood.
}
\usage{
gevreg(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE, z = FALSE, v = FALSE,
  control = gevreg_control(\dots), \dots)

gevreg_fit(x, y, z = NULL, v = NULL, n.stats, control)

gevreg_control(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, ...)


}

\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x | z | v}  where
    \code{y} is the response and \code{x}, \code{z} and \code{v} are regressor variables
    for the location, scale and shape parameters of the Generalized Extreme value Distribution (GEV). 
	For details how to set up the formula see details and \code{\link{Formula}}}.
  \item{data}{a data frame containing the covariables used to fit the respons(es). Variables in the formula must occur as column name in data. A column named \code{station} is mandatory. }
  \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y, z, v}{for \code{gevreg}: logical. If \code{TRUE} the model matrix and
    response matrix used for fitting are returned as part of the returned gevreg object. For \code{gevreg.fit}: \code{x, z, v} are design matrices with regressors for the
    location, scale and shape parameters of the GEV and \code{y} is a vector of observed parameters. }
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
  \item{control, maxit, start}{a list of control parameters passed to \code{\link{optim}} .}
  
  \item{n.stats}{number of stations that have to be fitted simultaneously.}
  \item{grad}{should the gradient-based method "BFGS" used for optimization. If \code{FALSE} "Nelder-Mead" is used instead.}
  \item{hessian}{should a hessian be computed during the optimaztion process. Possible values are \code{"none"}, \code{"numderiv"} or \code{"optim"}, which is the default.}
  
 
}

\details{
\code{gevreg}{ Fit a GEV distribution simultaneously to a number of observation sites. The response y on the left-hand-side (LHS) of the formula must be available in the data data.frame. The right-hand-side (RHS) of the formual can consist of one, two or three parts, separated by a "|". If only one part is given it is supposed to be the regressor for the location parameter of the GEV distribution. The scale and shape parameter are extended to a regressor of 1. If the formula contains two parts, they are considered as the regressors for the location and scale parameter (from left to right). The formula is then extended to a shape regressor of 1.}

\code{gevreg_fit}{ is the lower level function which provides the actual maximum likelihood fitting.}

\code{gevreg_control}{ takes control parameters e.g. maxit which are transferred to \link{optim}. In addition \code{grad} specifies wether a gradient based optimization routine ("BFGS") or "Nelder-Mead" should be used for optimization. In the former case the gradient will be computed with the \code{\link[numDeriv]{grad}} function.}


}

\value{
  An object of class \code{gevreg} which inherits from \code{optim} with components: 
  \item{coefficients}{Either all, or if specified e.g. as location, scale or shape the coefficients of the fitted model for the three GEV parameters for al stations are returned.}
  \item{loglik}{The maximised log-likelihood value.}
  \item{vcov}{returns the variance-covariance matrix as the generic function \link{vcov} does.} 
  \item{x, y, z,v}{If \code{TRUE} in the call the model matrix and response matrix used for fitting are returned.}
  
  
  
}

\references{
 	Blanchet J, Lehning M (2010). Mapping snow depth return levels: smooth spatial modeling versus station interpolation. 
	\emph{Hydrol. Earth Syst. Sci.}, \bold{14}, 2527--2544.
	\url{https://www.hydrol-earth-syst-sci.net/14/2527/2010/hess-14-2527-2010.pdf}.
}


\examples{
data("gevregdata")

## Fit GEV to the snow depths of the first station in the dataset
gevreg(hs~station,gevregdata)

## Use Nelder-Mead as optmization method
gevreg(hs~station,gevregdata,grad=FALSE)

## Compute hessian afterwards with the grad function of the numDeriv package
m <- gevreg(hs~station,gevregdata,grad=FALSE,hessian="numderiv",maxit=6000)
print(m, ext = TRUE)

## check some S3-Methods
m <- gevreg(hs~station,gevregdata)
vcov(m,model="loc")
m$vcov
coef(m)
coef(m, model="shape")
logLik(m)


}

\keyword{gev, regression}
