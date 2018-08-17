\name{ssdeRgmm.fit}
\alias{ssdeRgmm.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ssdeRgmm.fit}
\description{
Fits the main outcome equation using a method of moments estimator.
}
\usage{
ssdeRgmm.fit(theta, y1, y2, y3, x1, x2, x3, V1, weights = weights, control = control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
First stage parameter estimates.
}
  \item{y1}{
\code{treatment} response
}
  \item{y2}{
\code{selection} response
}
  \item{y3}{
\code{outcome}
}
  \item{x1}{
\code{treatment} covariates}
  \item{x2}{
\code{selection} covariates}
  \item{x3}{
\code{outcome} covariates}
  \item{V1}{
First stage variance-covariance matrix
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michael Brottrager}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theta, y1, y2, y3, x1, x2, x3, V1, weights = weights,
    control = control)
{
    ocontrol <- control
    robust <- control$robust
    delta <- control$reltol
    control$method <- control$start <- control$robust <- NULL
    n = NROW(x3)
    k = NCOL(x3)
    n_s = NROW(x2)
    k_s = NCOL(x2)
    n_d = NROW(x1)
    k_d = NCOL(x1)
    if (is.null(weights))
        weights <- rep.int(1, length(y3[!is.na(y3)]))
    else {
        weights <- weights
    }
    w <- weights[names(weights) \%in\% names(y3[!is.na(y3)])]
    nobs <- sum(weights > 0)
    NT1 <- sum(y1 == 1)
    NT0 <- nobs - NT1
    NS1 <- sum(y2 == 1)
    NS0 <- nobs - NS1
    if (NT0 == 0 | NT1 == 0) {
        stop("No variance in the treatment variable")
    }
    if (NS0 == 0 | NS1 == 0) {
        stop("No variance in the selection variable")
    }
    Y2Levels <- levels(as.factor(y2))
    if (length(Y2Levels) != 2) {
        stop("the left hand side of the 'formula' has to contain",
            " exactly two levels (e.g. FALSE and TRUE)")
    }
    y2 <- as.integer(y2 == Y2Levels[2])
    Y1Levels <- levels(as.factor(y1))
    if (length(Y1Levels) != 2) {
        stop("the left hand side of the 'formula' has to contain",
            " exactly two levels (e.g. FALSE and TRUE)")
    }
    y1 <- as.integer(y1 == Y1Levels[2])
    Gamma1 = theta[1:k_d]
    Gamma2 = theta[(k_d + 1):(k_d + k_s)]
    Beta2 = theta[k_d + k_s + 1]
    rho120 = theta[k_d + k_s + 2]
    rho121 = theta[k_d + k_s + 3]
    ZG1 = x1 \%*\% Gamma1
    ZG2 = x2 \%*\% Gamma2
    C11d = pbinorm(ZG1, ZG2 + Beta2, mean1 = 0, mean2 = 0, cov12 = rho121)
    C11 = dnorm(ZG1) * pnorm((ZG2 + Beta2 - rho121 * ZG1)/sqrt(1 -
        rho121^2))/C11d
    C12d = pbinorm(ZG1, ZG2 + Beta2, mean1 = 0, mean2 = 0, cov12 = rho121)
    C12 = dnorm(ZG2 + Beta2) * pnorm((ZG1 - rho121 * (ZG2 + Beta2))/sqrt(1 -
        rho121^2))/C12d
    Co1d = pbinorm(-ZG1, ZG2, mean1 = 0, mean2 = 0, cov12 = -rho120)
    Co1 = -dnorm(ZG1) * pnorm((ZG2 - rho120 * ZG1)/sqrt(1 - rho120^2))/Co1d
    Co2d = pbinorm(-ZG1, ZG2, mean1 = 0, mean2 = 0, cov12 = -rho120)
    Co2 = dnorm(ZG2) * pnorm((-ZG1 + rho120 * ZG2)/sqrt(1 - rho120^2))/Co2d
    y3 = y3[y2 != 0]
    y1 = y1[y2 != 0]
    C11 = C11[y2 != 0]
    C12 = C12[y2 != 0]
    Co1 = Co1[y2 != 0]
    Co2 = Co2[y2 != 0]
    Z1 = x1[y2 != 0, ]
    Z2 = x2[y2 != 0, ]
    Z3 = x3[y2 != 0, ]
    X = as.matrix(cbind(Z3, y1, y1 * C11, y1 * C12, (1 - y1) *
        Co1, (1 - y1) * Co2))
    colnames(X)[(k + 2):NCOL(X)] <- c("m_11", "m_12", "m_01",
        "m_02")
    rw <- sqrt(w)
    X_tilde <- rw * X
    y_tilde <- rw * y3
    theta_2 = solve(t(X_tilde) \%*\% (X_tilde)) \%*\% t(X_tilde) \%*\%
        (y_tilde)
    fitfun <- function(par) {
        beta <- par[seq.int(length.out = k + 5)]
        names(beta) <- rownames(par)
        mu <- drop(X \%*\% beta)
        list(beta = beta, mu = mu)
    }
    fit <- fitfun(theta_2)
    vcov <- ssdeR.vcv(X = X, theta = theta, theta_2 = theta_2,
        V1 = V1, k_s = k_s, k_d = k_d, Z1 = Z1, Z2 = Z2, Z3 = Z3,
        y1 = y1, y2 = y2, y3 = y3, delta = delta, robust = robust)
    ll <- 0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) +
        log(sum(w * (y3 - fit$mu)^2))))
    rval <- list(coefficients = fit$beta, residuals = y3 - fit$mu,
        fitted.values = fit$mu, loglik = ll, X = X, df.residual = n -
            k - 5, vcov = vcov, n = n, control = ocontrol, weights = w,
        param = list(nParam = k, nObs = n, NT1 = NT1, NT0 = NT0,
            NS1 = NS1, NS0 = NS0))
    class(rval) <- "ssdeR"
    return(rval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
