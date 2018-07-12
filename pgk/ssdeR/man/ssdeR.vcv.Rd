\name{ssdeR.vcv}
\alias{ssdeR.vcv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ssdeR.vcv}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ssdeR.vcv(X, fm, theta, theta_2, V1, k_s, k_d, Z1, Z2, Z3, y1, y2, y3, delta, robust)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{fm}{
%%     ~~Describe \code{fm} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{theta_2}{
%%     ~~Describe \code{theta_2} here~~
}
  \item{V1}{
%%     ~~Describe \code{V1} here~~
}
  \item{k_s}{
%%     ~~Describe \code{k_s} here~~
}
  \item{k_d}{
%%     ~~Describe \code{k_d} here~~
}
  \item{Z1}{
%%     ~~Describe \code{Z1} here~~
}
  \item{Z2}{
%%     ~~Describe \code{Z2} here~~
}
  \item{Z3}{
%%     ~~Describe \code{Z3} here~~
}
  \item{y1}{
%%     ~~Describe \code{y1} here~~
}
  \item{y2}{
%%     ~~Describe \code{y2} here~~
}
  \item{y3}{
%%     ~~Describe \code{y3} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{robust}{
%%     ~~Describe \code{robust} here~~
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
function (X, fm, theta, theta_2, V1, k_s, k_d, Z1, Z2, Z3, y1,
    y2, y3, delta, robust)
{
    n = nrow(X)
    n1 = length(y1)
    k = ncol(X)
    H = 1/(n) * t(X) \%*\% X
    e3 = y3 - X \%*\% theta_2
    VCV = 1/(n) * as.numeric(t(e3) \%*\% e3) * solve(t(X) \%*\% X)
    vcov = numDeriv.gmm(theta = theta, theta_2 = theta_2, H = H,
        V1 = V1, vcv = VCV, X = X, Z1 = Z1, Z2 = Z2, Z3 = Z3,
        y1 = y1, delta = delta, n = n, n1 = n1, k = k, k_s = k_s,
        k_d = k_d)
    return(vcov)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
