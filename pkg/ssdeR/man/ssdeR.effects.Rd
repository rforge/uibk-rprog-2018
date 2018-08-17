\name{ssdeR.effects}
\alias{ssdeR.effects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ssdeR.effects}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ssdeR.effects(object, model = c("outcome", "treatment", "selection"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
\code{ssdeR} object
}
  \item{model}{
Either outcome, treatment or selection model.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the (marginal) effects of each model. If no \code{model} is chosen, \code{ssdeR.effects} just returns the parameter estimates of the outcome equation.
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
function (object, model = c("outcome", "treatment", "selection"))
{
    model <- match.arg(model)
    if (model \%in\% c("treatment", "selection")) {
        direct.effect <- ssdeR.mfx(object, model)
        std.err <- c()
        for (i in names(direct.effect)) {
            v1 <- object$firststage$vcov
            v1 <- v1[paste0(i), paste0(i)]
            v1 <- direct.effect[i] \%*\% v1 \%*\% direct.effect[i]
            std.err[paste0(i)] = sqrt(v1)
        }
    }
    else {
        direct.effect <- coef(object)[2:(length(coef(object)) -
            4)]
        std.err <- c()
        for (i in names(direct.effect)) {
            v1 <- vcov(object)
            v1 <- v1[paste0(i), paste0(i)]
            std.err[paste0(i)] = sqrt(v1)
        }
    }
    effects <- cbind(direct.effect, std.err)
    return(effects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
