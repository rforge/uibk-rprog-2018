\name{summary.ssdeR}
\alias{summary.ssdeR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary.ssdeR}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.ssdeR(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
%%  ~~who you are~~
}
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
function (object, ...)
{
    k2 <- length(object$coefficients)
    k1 <- length(object$firststage$estimate)
    ko <- NCOL(model.matrix(object, "o"))
    kt <- NCOL(model.matrix(object, "t"))
    ks <- NCOL(model.matrix(object, "s"))
    cf2 <- object$coefficients
    cf1 <- object$firststage$estimate
    mus <- as.vector(do.call("c", as.list(cf2[(ko + 1):(k2)])))
    rhos <- as.vector(do.call("c", as.list(cf1[(kt + ks + 1):(k1)])))
    cfo <- as.vector(do.call("c", as.list(cf2[1:ko])))
    cft <- as.vector(do.call("c", as.list(cf1[1:kt])))
    cfs <- as.vector(do.call("c", as.list(cf1[(kt + 1):(kt +
        ks)])))
    seo <- sqrt(diag(object$vcov))
    seo_covar <- seo[1:(ko)]
    seo_mus <- seo[(ko + 1):(k2)]
    set <- sqrt(diag(object$firststage$vcov))
    set_t <- set[1:(kt)]
    set_s <- set[(kt + 1):(kt + ks)]
    set_rho <- set[(kt + ks + 1):(k1)]
    cf <- c(cft, cfs, cfo)
    se <- c(set_t, set_s, seo_covar)
    cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
    AP <- c(rhos, mus)
    se <- c(set_rho, seo_mus)
    AP <- cbind(AP, se, AP/se, 2 * pnorm(-abs(AP/se)))
    listDF <- list(Parameters = cf, Auxiliaries = AP)
    new_col_name <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    listDF <- lapply(listDF, function(i) {
        colnames(i) <- new_col_name
        i
    })
    cf <- listDF[[1]]
    AP <- listDF[[2]]
    if (length(object$coefficients)) {
        cf <- list(treatment = cf[seq.int(length.out = kt), ,
            drop = FALSE], selection = cf[seq.int(length.out = ks) +
            kt, , drop = FALSE], outcome = cf[seq.int(length.out = ko) +
            kt + ks, , drop = FALSE], Aux.Param = AP[, , drop = FALSE])
    }
    else {
        warning("\nNo Parameters found\n")
    }
    object$AIC <- AIC(object)
    object$BIC <- AIC(object, k = sqrt(nObs(object)))
    object$DF1 <- k2
    object$DF2 <- k1
    object$coefficients <- cf
    object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL
    class(object) <- "summary.ssdeR"
    object
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
