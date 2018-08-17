\name{model.matrix.ssdeR}
\alias{model.matrix.ssdeR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model.matrix.ssdeR}
\description{
extracts the model matrix of \code{ssdeR} objects. If \code{model} is not specified, it returns the \code{model.matrix} of the outcome model. \code{model.matrix.ssdeR} returns the model matrices without the estimated auxiliary parameters. If you want to extract those as well, call \code{object$X}.
}
\usage{
model.matrix.ssdeR(object,
                  model = c("outcome",
                            "treatment",
                            "selection"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
\code{ssdeR} object
}
  \item{model}{
either \code{outcome, treamtment } or \code{selection}.
}
  \item{\dots}{
further controls
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
retuns \eqn{n x k} model matrices of the respective models.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michael Brottrager
}
\note{
note that model matrices contain no auxiliary parameters
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## mm
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
