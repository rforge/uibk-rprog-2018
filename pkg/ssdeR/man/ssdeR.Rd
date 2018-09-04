\name{ssdeR}

\alias{ssdeR}
\alias{aic_ssdeR}
\alias{estfun_ssdeR}
\alias{fitted.ssdeR}
\alias{print.ssdeR}
\alias{summary.ssdeR}
\alias{print.summary.ssdeR}
\alias{model.matrix.ssdeR}
\alias{logLik.ssdeR}
\alias{model.frame.ssdeR}
\alias{nObs.ssdeR}
\alias{estfun.ssdeR}
\alias{predict.ssdeR}
\alias{residuals.ssdeR}
\alias{numDeriv.gmm}
\alias{ssdeR.cluster.fit}
\alias{ssdeR.control}
\alias{ssdeR.mfx}
\alias{ssdeR.vcv}
\alias{ssdeRgmm.fit}
\alias{terms.ssdeR}
\alias{vcov.ssdeR}
\alias{ssdeRbiprobit.fit}
\alias{nobs.ssdeR}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample selection models with a common dummy endogenous regressor in simultaneous equations: A simple two-step estimation}
\description{
Estimates sample selection models where a common dummy endogenous regressor appears both in the selection equation and in the censored equation. We interpret this model as an endogenous switching model and develop a simple two step estimation procedure. For model derivation and see Kim (2006) (\url{https://doi.org/10.1016/j.econlet.2005.12.003})
}
\usage{
ssdeR(formula, treatment, selection, data, subset,
      na.action = FALSE, weights, cluster = NULL,
      print.level = 0, control = ssdeR.control(...),
      model = TRUE, x = FALSE, y = FALSE, ...)

ssdeR.control(method = "BHHH", iterlim = NULL, start = NULL, robust = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  formula, Outcome equation (Continuous dependent variable).
}
  \item{treatment}{
  formula, Treatment Equation (Binary dependent variable).
}
  \item{selection}{
  formula Selection equation (Binary dependent variable).
}
  \item{data}{
  an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{ssdeR} is called.
}
  \item{subset}{
an optional index vector specifying a subset of observations to be used in the fitting process.
}
  \item{na.action}{
  Restricted to \code{na.pass}
}
  \item{weights}{
an optional vector of ‘prior weights’ to be used in the fitting process. Should be \code{NULL} or a numeric vector.
}
  \item{cluster}{
Character vector of up to 2 cluster variables.)
}
  \item{print.level}{
integer. Various debugging information, higher value gives more information. Not supported in this build.}
  \item{control}{
Further controls for \code{maxLik} maximization etc. (see \cite{maxLik})
}
  \item{model}{
keep \code{model.frame} if \code{model = T}.
}
  \item{x}{
keep independent second stage variables if \code{x = T}.
}
  \item{y}{
keep dependent outcome varibale if \code{y = T}.
}
\item{method}{
Maximisation method used in the bivariate probit model (First Stage). Default is "BHHH"
}
\item{iterlim}{
User specified maximal number of iterations. Default is 5000.
}
\item{start}{
User specified vector of starting values.
}
\item{robust}{
Robust standard errors.
}

  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
This package provides a estimation function for sample selection models where a common dummy endogenous regressor appears both in the selection equation and in the censored equation. This model is analyzed in the framework of an endogenous switching model. Following Kim (2006), a simple two-step estimator is used for this model, which is easy to implement and numerically robust compared to other methods.
}

\value{
\code{ssdeR()} returns an object of class \code{ssdeR}.

The first stage model (\code{firststage}) is estimated by Maximum Likelihood, which has all the components of a 'maxLik' object, and \code{vcov}, which contains the (cluster-) robust variance-covariance matrix of the first stage model.

Furthermore, the returned 'ssdeR' object returns the following list of components:

\item{coefficients}{estimated coefficients, the outcome model.
    coefficient for the auxiliary parameters \eqn{\mu_{ij}} are treated as a parameters.}

\item{residuals}{estimated residuals, the outcome model.}

\item{fitted.values}{fitted values, the outcome model.}

\item{loglik}{log likelihood of the outcome model.}

\item{df.residual}{degrees of freedom of the outcome model.}

\item{vcov}{variance covariance matrix of the estimated coefficients.}

\item{n}{total number of used observations 1st and 2nd stage.}

\item{controls}{\code{List} of controls applied to the \code{fit} functions.}

\item{weights}{Vector of weights. If no weightes were supported, \code{weights} is just a vector of ones with length \code{n =} censored observations.}

\item{param}{\code{List} object. \code{nParam} is the number of covariates used in the outcome model. \code{nObs} is the total number of used observations, \code{NT1} and \code{NT0} refer to the number of treated and untreated, respectively. \code{NS1, NS2} refer to the number of censored and uncensored observations, respectively.}

\item{firststage}{object of class 'maxLik' that contains the results of the
    1st step (bivariate probit estimation) and the (cluster-) robust variance-covariance matrix (if requested).}
}

\references{
  Cameron, A. C. and Trivedi, P. K. (2005)
  \emph{Microeconometrics: Methods and Applications},
  Cambridge University Press.

  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.

  Heckman, J. (1976)
  The common structure of statistical models of
  truncation, sample selection and limited dependent
  variables and a simple estimator for such models,
  \emph{Annals of Economic and Social Measurement},
  5(4), p. 475-492.

  Johnston, J. and J. DiNardo (1997)
  \emph{Econometric Methods, Fourth Edition}, McGraw-Hill.

  Lee, L., G. Maddala and R. Trost (1980)
  Asymetric covariance matrices of two-stage probit and two-stage tobit
  methods for simultaneous equations models with selectivity.
  \emph{Econometrica}, 48, p. 491-503.

  Mullahy, J. (2017) Marginal effects in multivariate probit models.  \emph{Empircal Economics}, 52: 447.

  il Kim, K. (2006). Sample selection models with a common dummy endogenous regressor in simultaneous equations: A simple two-step estimation. \emph{Economics Letters}, 91(2), 280-286.


  Petersen, S., G. Henningsen and A. Henningsen (2017)
  \emph{Which Households Invest in Energy-Saving Home Improvements?
  Evidence From a Danish Policy Intervention}.
  Unpublished Manuscript.
  Department of Management Engineering, Technical University of Denmark.

  Toomet, O. and A. Henningsen, (2008)
  Sample Selection Models in R: Package sampleSelection.
  \emph{Journal of Statistical Software} 27(7),
  \url{http://www.jstatsoft.org/v27/i07/}

  Wooldridge, J. M. (2003)
  \emph{Introductory Econometrics: A Modern Approach, 2e}, Thomson South-Western.}

\author{
Michael Brottrager
}


\seealso{
See \code{\link[maxLik]{maxLik}} or \code{\link[sampleSelection]{selection}} for further information. }

\examples{
# ----------------------------------------------------------------------------------- #
# 1. Climate, Conflict, Crossborder Migration Data
# ----------------------------------------------------------------------------------- #
# This is the data used Abel, Crespo-Cuaresma, Brottrager, Muttarak (2018).
# As the paper is currently under revision, readers are recommended to directly
# contact <michael.brottrager@jku.at> for a current version of the paper including
# the detailed data description.
#
# Note that running the model on the data provided below takes some
# considerable amount of time.
# ----------------------------------------------------------------------------------- #
 data(ConflictMigration, package="ssdeR")

 Results <- ssdeR(formula = asylum_seekers_ij ~ stock_ij + dist_ij + I(dist_ij^2) +
                            comlang_ij + colony_ij  +
                            polity_i + pop_i + polity_i + pop_i +
                            gdp_j,
                  treatment = conflict_i ~ battledeaths_i + spei_i +
                            polity_i + I(polity_i^2) +
                            diaspora_i + ethMRQ_i ,
                  selection = isflow_ij ~ dist_ij + I(dist_ij^2) +
                            outmigration_i + inmigration_j ,
                            cluster = c("iso_i","iso_j"),
                 data = ConflictMigration)

 summary(Results)
# ----------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------- #
# 2. Simulation Data
# ----------------------------------------------------------------------------------- #
# library(MASS)
# set.seed(12072018)
# n <- 5000
# gamma1 <- c(-0.2, 0.2, 0.1)
# gamma2 <- c(0.2, -0.1, 0.2)
# gamma3 <- c(2, 0.5, 0.3)
# beta2 <- 1
# beta3 <- -2
# rho120 <- -0.5
# rho121 <- 0.5
#
# muZ1 <- c(1,-1)
# muZ2 <- c(1,1)
# muZ3 <- c(2,0)
#
# s_Z1 <- matrix(c(3, 0.1, 0.1, 3), nrow = 2, ncol = 2)
# s_Z2 <- matrix(c(3, 0.1, 0.1, 3), nrow = 2, ncol = 2)
# s_Z3 <- matrix(c(4, 0.3, 0.3, 4), nrow = 2, ncol = 2)
#
# Z1 <- matrix(c(rep(1, n),
#                mvrnorm(n = n,
#                              mu = muZ1,
#                              Sigma = s_Z1)),
#              nrow = n,
#              ncol = 3)
# Z2 <- matrix(c(rep(1, n),
#                mvrnorm(n = n,
#                              mu = muZ2,
#                              Sigma = s_Z2)),
#              nrow = n,
#              ncol = 3)
# Z3 <- matrix(c(rep(1, n),
#                mvrnorm(n = n,
#                              mu = muZ3,
#                              Sigma = s_Z3)),
#              nrow = n,
#              ncol = 3)
#
# eps <- rnorm(n, 0, 1)
#
#
# y1star <- Z1 %*% gamma1 + eps
# y1 <- as.numeric(y1star>0)
#
#
# eps0 <- rho120*eps + sqrt(1-rho120^2)*rnorm(n, 0, 1)
# eps1 <- rho121*eps + sqrt(1-rho121^2)*rnorm(n, 0, 1)
#
# y2star <- Z2 %*% gamma2 + y1*beta2 + y1*eps1 + (1-y1)*eps0
# y2 <- as.numeric(y2star>0)
#
# rho230 <- 0.5
# rho231 <- 0.3
#
# y3star <- Z3 %*% gamma3 + beta3*y1 + rnorm(n, 0, 1.5) +
#   y1*eps1*rho231 + (1-y1)*eps0*rho230
# cens <- ifelse(y2 == 1, 1, NA)
# y3 <- y3star*cens
#
# df <- matrix(cbind(y1, y2, y3, Z1[,c(2:3)],Z2[,c(2:3)],Z3[,c(2:3)]), nrow = n,
#              dimnames = list(c(1:n), c("y1","y2", "y3",
#                                        "d1", "d2", "s1",
#                                        "s2", "x1", "x2")))
#
# df <- as.data.frame(do.call("rbind", replicate(1, df, simplify = FALSE)))
# rownames(df) <- 1:(n)
#
# m1 <- ssdeR(formula = y3 ~ x1 + x2,
#             treatment = y1 ~ d1 +d2,
#             selection = y2 ~ s1 +s2,
#             data = df)
#
# summary(m1)
# ----------------------------------------------------------------------------------- #
}

\keyword{ regression }
