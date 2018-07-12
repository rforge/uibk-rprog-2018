\name{ssdeR.cluster.fit}
\alias{ssdeR.cluster.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ssdeR.cluster.fit}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ssdeR.cluster.fit(fm, dfcw, cluster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fm}{
%%     ~~Describe \code{fm} here~~
}
  \item{dfcw}{
%%     ~~Describe \code{dfcw} here~~
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
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
function (fm, dfcw, cluster)
{
    if (NCOL(cluster) == 1) {
        M <- length(unique(cluster))
        N <- length(cluster)
        dfc <- (M/(M - 1)) * ((N - 1)/(N - length(coef(fm))))
        u <- apply(estfun(fm), 2, function(x) tapply(x, cluster,
            sum))
        vcovCL <- dfc * sandwich(fm, meat = crossprod(u)/N) *
            dfcw
        return(vcovCL)
    }
    else {
        if (NCOL(cluster) > 2) {
            warning("max 2-way clustering")
        }
        cluster1 = cluster[rownames(cluster) \%in\% names(weights(fm)),
            1]
        cluster2 = cluster[rownames(cluster) \%in\% names(weights(fm)),
            2]
        cluster12 = paste(cluster1, cluster2, sep = "")
        M1 <- length(unique(cluster1))
        M2 <- length(unique(cluster2))
        M12 <- length(unique(cluster12))
        N <- length(cluster1)
        K <- length(coef(fm))
        dfc1 <- (M1/(M1 - 1)) * ((N - 1)/(N - K))
        dfc2 <- (M2/(M2 - 1)) * ((N - 1)/(N - K))
        dfc12 <- (M12/(M12 - 1)) * ((N - 1)/(N - K))
        u1 <- apply(estfun(fm), 2, function(x) tapply(x, cluster1,
            sum))
        u2 <- apply(estfun(fm), 2, function(x) tapply(x, cluster2,
            sum))
        u12 <- apply(estfun(fm), 2, function(x) tapply(x, cluster12,
            sum))
        vc1 <- dfc1 * sandwich(fm, meat = crossprod(u1)/N)
        vc2 <- dfc2 * sandwich(fm, meat = crossprod(u2)/N)
        vc12 <- dfc12 * sandwich(fm, meat = crossprod(u12)/N)
        vcovMCL <- (vc1 + vc2 - vc12) * dfcw
        return(vcovMCL)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
