\name{nObs.ssdeR}
\alias{nObs.ssdeR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
nObs.ssdeR}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nObs.ssdeR(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object for which the extraction of the number of observations is meaningful.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a scalar (number of observations)
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
function (object, ...)
{
    object$param$nObs
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
