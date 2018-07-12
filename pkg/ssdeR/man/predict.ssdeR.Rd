\name{predict.ssdeR}
\alias{predict.ssdeR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
predict.ssdeR}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.ssdeR(object, model = c("outcome", "treatment", "selection"), pnorm = F, newdata = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{pnorm}{
%%     ~~Describe \code{pnorm} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
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
function (object, model = c("outcome", "treatment", "selection"),
    pnorm = F, newdata = NULL, ...)
{
    model <- match.arg(model)
    if (is.null(newdata)) {
        mm <- model.matrix(object, model)
    }
    else {
        tt <- terms(object, model)
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, xlev = object$xlevels)
        mm <- model.matrix(Terms, m, model)
    }
    if (model == "outcome") {
        result <- drop(object$X \%*\% coef(object))
    }
    else {
        nm <- colnames(mm)
        result <- drop(mm \%*\% coef(object$firststage)[nm])
    }
    if (model \%in\% c("treatment", "selection", pnrom = T)) {
        result <- pnorm(result)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
