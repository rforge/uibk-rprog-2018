\name{ssdeR.control}
\alias{ssdeR.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ssdeR.control}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ssdeR.control(method = "BHHH", iterlim = NULL, start = NULL, robust = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{iterlim}{
%%     ~~Describe \code{iterlim} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{robust}{
%%     ~~Describe \code{robust} here~~
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
function (method = "BHHH", iterlim = NULL, start = NULL, robust = F,
    ...)
{
    if (is.null(iterlim))
        iterlim <- 5000
    rval <- list(method = method, iterlim = iterlim, start = start,
        robust = F, biprob.fit = "ssdeRbiprob.fit", gmm.fit = "ssdeRgmm.fit")
    rval <- c(rval, list(...))
    if (is.null(rval$reltol))
        rval$reltol <- .Machine$double.eps^(1/1.2)
    rval
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
