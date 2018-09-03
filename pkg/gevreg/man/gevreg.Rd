\name{gevreg}
\alias{gevreg}
\alias{gevreg_fit}
\alias{gevreg_control}


\title{Smooth Spatial Maximum Likelihood GEV Fitting}

\description{
  Fit a GEV distribution spatially to observations using maximum likelihood.
}
\usage{
gevreg(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = TRUE, z = FALSE, v = FALSE, gev_params, 
  control = gevreg_control(\dots), \dots)

gevreg_fit(x, y, z = NULL, v = NULL, n.stats, n.years, gevp, control)

gevreg_control(maxit = 5000, start = NULL, grad = TRUE, ...)

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
    contain \code{NA}s. Currently na.action = NULL has to be specified.}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y, z, v}{for \code{gevreg}: logical. If \code{TRUE} the model matrix and
    response matrix used for fitting are returned as part of the returned gevreg object. For \code{gevreg.fit}: \code{x, z, v} are design matrices with regressors for the
    location, scale and shape parameters of the GEV and \code{y} is a vector of observed parameters. }
    \item{gev_params, gevp}{ GEV parameters at the stations. A data frame with three columns named \emph{loc}, \emph{scale} and \emph{shape}. Each row corresponds to one onservation site. If not provided, gev_params is created from the data with the function \code{\link{gevmle}} from package \code{\link{SpatialExtremes}}.}
    
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
  \item{control, maxit, start, grad}{a list of control parameters passed to \code{\link{optimx}}. If \code{grad=TRUE}, the gradient-based method "nlminb" is used for optimization, otherwise "Nelder-Mead" is used.}
  \item{n.stats,n.years}{for \code{gevreg_fi}: n.stats: number of stations that shall be fitted simultaneously. n.years: nuber of years per station. All stations must have the same number of years. }
  
 
}

\details{
\code{gevreg}{ Fit a GEV distribution spatially to a number of observation sites simultaneously. The spatial context is established by maximising the sum of the log-likelihoods of all stations. The response y on the left-hand-side (LHS) of the formula must be available in the data data.frame. The right-hand-side (RHS) of the formual can consist of one, two or three parts, separated by a "|". If only one part is given it is supposed to be the regressor for the location parameter of the GEV distribution. The scale and shape parameter are extended to a regressor of 1. If the formula contains two parts, they are considered as the regressors for the location and scale parameter (from left to right). The formula is then extended to a shape regressor of 1.}

\code{gevreg_fit}{ is the lower level function which provides the actual maximum likelihood fitting. For a more stable fitting procedure, a "squared log-likelihood" ist maximised if the standard likelihood fails.}

\code{gevreg_control}{ takes control parameters e.g. maxit which are transferred to \link{optimx}. In addition \code{grad} specifies wether a gradient based optimization routine ("nlminb") or "Nelder-Mead" should be used for optimization. }

}

\value{
  An object of class \code{gevreg} which inherits from \code{optimx} with components: 
  \item{coefficients}{ Either all, or if specified e.g. as location, scale or shape the coefficients of the fitted model for the three GEV parameters are returned.}
  \item{gev_params}{The GEV parameters used to fit the model.}
  \item{loglik}{ The maximised log-likelihood value.}
  \item{x, y, z,v}{ If \code{TRUE} in the call the model matrix and response matrix used for fitting are returned.}
  Refer to \code{\link{optimx}} for a list of return values inherited from there (\code{niter},\code{feval},\code{geval},\code{xtimes},\code{convcode},... ).
  
  
}

\references{
 	Blanchet J, Lehning M (2010). Mapping snow depth return levels: smooth spatial modeling versus station interpolation. 
	\emph{Hydrol. Earth Syst. Sci.}, \bold{14}, 2527--2544.
	\url{https://www.hydrol-earth-syst-sci.net/14/2527/2010/hess-14-2527-2010.pdf}.
}


\examples{

## Load example data
data("SnowAustria")
data("SnowAustriaGEV")

## explore dataset
with(SnowAustria, plot(lon,lat))
with(SnowAustria, hist(hs))
df <- cbind(SnowAustriaGEV,data.frame(lon = unique(SnowAustria$lon), alt = unique(SnowAustria$alt)))
plot(df)

## Fit GEV to the snow depths of all stations in the dataset
## by explaining the GEV parameters with smooth functions of spatial covariates 
## and the location parameter
SnowAustria$loc <- as.numeric(rep(SnowAustriaGEV$loc,
                                  each=nrow(SnowAustria)/length(levels(SnowAustria$station))))
m0 <- gevreg(formula = hs ~ lon + alt | lat + alt +loc | lon, 
             data = SnowAustria, gev_params = SnowAustriaGEV, na.action = NULL)

## Compare different models
m1 <- gevreg(hs ~ 1 | loc | alt, data = SnowAustria, gev_params = SnowAustriaGEV, na.action = NULL)
cbind(AIC(m0,m1),loglik=c(logLik(m0),logLik(m1)))

## compare observed and modeled return values for return period of 50 years
rl.modeled <- predict(m0, type = "quantile", at = 1-1/50)
rl.observed <- c()
if(require("extRemes")) {
  for (i in 1:length(levels(SnowAustria$station))){
  d <- as.numeric(na.omit(SnowAustria$hs[which(SnowAustria$station==i)]))
  fit <- fevd(d,type="GEV",method="MLE")
  rl <- return.level(fit,return.period = 50)
  rl.observed <- c(rl.observed,rl)
}
}
plot(rl.modeled,rl.observed)
abline(1,1)


## predict GEV parameters for new location
ungauged <- data.frame(lon=c(14.4,9.58,11.2),
                       lat=c(44.6,44.3,41.23),
                       alt=c(867,1600,2145),loc=c(50,60,80))
predict(m0, type = "parameter", newdata = ungauged)

}

\keyword{gev, regression, spatial}
