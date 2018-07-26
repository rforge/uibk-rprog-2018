\name{logitr}
\alias{logitr}
\alias{logitr_fit}
\alias{logitr_control}

\alias{print.logitr}
\alias{summary.logitr}
\alias{print.summary.logitr}
\alias{update.logitr}
\alias{logLik.logitr}
\alias{model.frame.logitr}

\title{Logistic Regression}

\description{
  Fitting logit regression models.
}
\usage{
logitr(formula, data, subset, na.action,
  model = TRUE, x = FALSE, y = TRUE,
  control = logitr_control(\dots), \dots)

logitr_fit(x, y, control)

logitr_control(maxit = 5000, start = NULL, hessian = TRUE, \dots)
}
\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x} where
    \code{y} is the response and \code{x} are the regressor variables.}
  \item{data}{an optional data frame containing the variables occurring in the
    formulas.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y}{for \code{logitr}: logical. If \code{TRUE} the model matrix and
    response vector used for fitting are returned as components of the returned value.
    For \code{logitr_fit}: \code{x} is a design matrix with the regressors
    and \code{y} is a vector of observations.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
  \item{control, maxit, start}{a list of control parameters passed to \code{\link[stats]{optim}} .}
  \item{hessian}{logical or character. Should a numeric approximation of the
    (negative) Hessian matrix be computed? Either \code{FALSE} (or equivalently
    \code{"none"}) or \code{TRUE}. Alternatively, in the latter case,
    \code{hessian = "numDeriv"} could be specified to signal that the Hessian should
    be approximated by \code{\link[numDeriv]{hessian}}. Another option is
    \code{hessian = "optim"} so that \code{\link[stats]{optim}} is used
    for computing the Hessian.}
}

\details{
  \code{logitr} fits logit regression models using maximum likelihood
  estimation. The model assumes an underlying latent binomial variable.

  \code{logitr_fit} is the lower level function where the actual fitting takes place.

  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"logitr"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}},
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}},
  \code{\link[stats]{predict}}, \code{\link[stats]{terms}},
  \code{\link[stats]{model.frame}}, \code{\link[stats]{model.matrix}}, \code{\link[stats]{update}}.

  This is a simpler implementation of \code{\link[stats]{glm}} only modeling
  the \code{"family = binomial"} used for logistic regression.
  \code{logitr} was written for the sole purpose of learning how to write packages
  and is not recommended to ever be used instead of \code{glm}.
}

\value{
  \code{logitr} returns an object of class \code{"logitr"}, i.e., a list with components as follows.
  \code{logitr_fit} returns an unclassed list with components up to \code{df}.
  \item{coefficients}{a list containing the coefficients,}
  \item{counts}{count of function and gradient evaluations from \code{optim},}
  \item{convergence}{convergence code from \code{optim},}
  \item{message}{optional further information from \code{optim},}
  \item{vcov}{covariance matrix of all parameters in the model,}
  \item{residuals.pearson}{a vector containing the Pearson Residuals of the model,}
  \item{fitted.values}{a list containing the latent fitted values,}
  \item{method}{the method argument passed to the \code{optim} call,}
  \item{nobs}{number of observations,}
  \item{df}{number of estimated parameters,}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}
  \item{terms}{a list containing the terms objects for the model,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the numeric response vector (if \code{y = TRUE}),}
  \item{x}{a list containing the model matrices
    (if \code{x = TRUE}).}
}

\examples{
##  quick randomly generated dataset. Comparison package vs glm
set.seed(123)
x1 <- rnorm(30,3,2) + 0.1 * c(1:30)
x2 <- rbinom(30,1,0.3)
x3 <- rpois(n=30,lambda = 4)
x3[16:30] <- x3[16:30] - rpois(n=15, lambda = 2)
xdat <- cbind(x1,x2,x3)
ydat <- c(rbinom(5,1,0.1), rbinom(10,1,0.25), rbinom(10,1,0.75), rbinom(5,1,0.9))
(m0 <- logitr(ydat~xdat))
(m1 <- glm(ydat~xdat, family = "binomial"))
## comparing AIC and BIC
AIC(m0)
AIC(m1)
BIC(m0)
BIC(m1)
}

\keyword{regression}
