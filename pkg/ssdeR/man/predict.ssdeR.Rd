\name{predict.ssdeR}
\alias{predict.ssdeR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
predict.ssdeR}
\description{
Computes predicted values of \code{ssdeR} objects. If no data is supplied, \code{predict.ssder()} simply uses the model data.}
\usage{
predict.ssdeR(object, model = c("outcome", "treatment", "selection"), pnorm = F, newdata = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
\code{ssdeR} object}
  \item{model}{
Either outcome, treatment or selection model.}
  \item{pnorm}{
If \code{pnorm == T}, \code{predict.ssdeR} uses the \pnorm.
}
  \item{newdata}{
User supplied new data used for prediction
}
  \item{\dots}{
further arguments}
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
