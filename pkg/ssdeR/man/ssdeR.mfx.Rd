\name{ssdeR.mfx}
\alias{ssdeR.mfx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ssdeR.mfx}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ssdeR.mfx(x, model = c("outcome", "treatment", "selection"), by = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{ssdeR} object
}
  \item{model}{
Either outcome, treatemtn or selection model.}
  \item{by}{
%%     ~~Describe \code{by} here~~
}
}
\details{
Computes the marginal effects of the bivatiate probit model. If \code{model == "outcome"}, \code{ssdeR.mfx} returns the parameter estimates of the outcome model.}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Mullahy, J. Empir Econ (2017) 52: 447. https://doi.org/10.1007/s00181-016-1090-8
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
function (x, model = c("outcome", "treatment", "selection"),
    by = NULL)
{
    model <- match.arg(model)
    checkBinaryTrait = function(v, naVal = "NA") {
        if (!is.numeric(v))
            stop("Only numeric vectors are accepted.")
        vSet = unique(v)
        if (!missing(naVal))
            vSet[vSet == naVal] = NA
        vSet = vSet[!is.na(vSet)]
        if (any(as.integer(vSet) != vSet))
            F
        else if (length(vSet) > 2)
            F
        else T
    }
    if (model == "outcome") {
        warning("no marginal effects for outcome equation - Direct/Indirect Effect estimates are provided via ssdeR.effects()")
    }
    else {
        if (model \%in\% c("selection")) {
            k <- length(coef(x$firststage))
            betas <- coef(x$firststage)[colnames(model.matrix(x,
                model))]
            rho121 <- tail(coef(x$firststage), 1)
            mm <- model.matrix(x, model)
            mmt <- model.matrix(x, "t")
            m <- c()
            for (i in 1:length(betas)) {
                b <- betas[i]
                if (checkBinaryTrait((mm[, names(b)]))) {
                  mm[, names(b)] <- 0
                  pdf0 <- mean((dnorm(predict(x, "selection",
                    pnorm = F, newdata = as.data.frame(mm))) *
                    dnorm(predict(x, "treatment", pnorm = F,
                      newdata = as.data.frame(mmt)) - rho121 *
                      predict(x, "selection", pnorm = F, newdata = as.data.frame(mm)))/sqrt(1 -
                    rho121^2)), na.rm = T) + mean((dnorm(predict(x,
                    "treatment", pnorm = F, newdata = as.data.frame(mmt))) *
                    dnorm(predict(x, "selection", pnorm = F,
                      newdata = as.data.frame(mm)) - rho121 *
                      predict(x, "treatment", pnorm = F, newdata = as.data.frame(mmt)))/sqrt(1 -
                    rho121^2)), na.rm = T)
                  mm[, names(b)] <- 1
                  pdf1 <- mean((dnorm(predict(x, "selection",
                    pnorm = F, newdata = as.data.frame(mm))) *
                    dnorm(predict(x, "treatment", pnorm = F,
                      newdata = as.data.frame(mmt)) - rho121 *
                      predict(x, "selection", pnorm = F, newdata = as.data.frame(mm)))/sqrt(1 -
                    rho121^2)), na.rm = T) + mean((dnorm(predict(x,
                    "treatment", pnorm = F, newdata = as.data.frame(mmt))) *
                    dnorm(predict(x, "selection", pnorm = F,
                      newdata = as.data.frame(mm)) - rho121 *
                      predict(x, "treatment", pnorm = F, newdata = as.data.frame(mmt)))/sqrt(1 -
                    rho121^2)), na.rm = T)
                  m[names(b)] <- pdf1 * b - pdf0 * b
                }
                else {
                  pdf <- mean((dnorm(predict(x, "selection",
                    pnorm = F)) * dnorm(predict(x, "treatment",
                    pnorm = F) - rho121 * predict(x, "selection",
                    pnorm = F))/sqrt(1 - rho121^2)), na.rm = T) +
                    mean((dnorm(predict(x, "treatment", pnorm = F)) *
                      dnorm(predict(x, "selection", pnorm = F) -
                        rho121 * predict(x, "treatment", pnorm = F))/sqrt(1 -
                      rho121^2)), na.rm = T)
                  m[names(b)] <- pdf * b
                }
            }
        }
        if (model == "treatment") {
            betas <- coef(x$firststage)[colnames(model.matrix(x,
                model))]
            mm <- model.matrix(x, model)
            m <- c()
            for (i in 1:length(betas)) {
                b <- betas[i]
                if (checkBinaryTrait((mm[, names(b)]))) {
                  mm[, names(b)] <- 0
                  pdf0 <- mean(dnorm(predict(x, "treatment",
                    pnorm = F, newdata = as.data.frame(mm))),
                    na.rm = T)
                  mm[, names(b)] <- 1
                  pdf1 <- mean(dnorm(predict(x, "treatment",
                    pnorm = F, newdata = as.data.frame(mm))),
                    na.rm = T)
                  pdf <- pdf1 - pdf0
                }
                else {
                  pdf <- mean(dnorm(predict(x, "treatment", pnorm = F)),
                    na.rm = T)
                }
                m[names(b)] <- pdf * b
            }
        }
        ifelse("(Intercept)" \%in\% names(m), return(m[which(names(m) !=
            "(Intercept)")]), return(m))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
