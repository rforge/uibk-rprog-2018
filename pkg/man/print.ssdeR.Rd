\name{print.ssdeR}
\alias{print.ssdeR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
print.ssdeR}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.ssdeR(x, digits = max(3, getOption("digits") - 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{ssdeR} object
}
  \item{digits}{
Number of digits printed
}
  \item{\dots}{
Further arguments}
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
function (x, digits = max(3, getOption("digits") - 3), ...)
{
    k <- length(x$coefficients)
    kt <- length(x$firststage$estimate)
    ks <- length(x$firststage$estimate)
    mmo <- model.matrix(x, "o")
    mmt <- model.matrix(x, "t")
    mms <- model.matrix(x, "s")
    cfo <- x$coefficients
    cft <- x$firststage$estimate[1:NCOL(mmt)]
    cfs <- x$firststage$estimate[(NCOL(mmt) + 1):(NCOL(mmt) +
        (NCOL(mms) + 2))]
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") *
        0.85)), "", sep = "\n")
    if (k) {
        cat(paste("Outcome Coefficients:\n", sep = ""))
        print.default(format(cfo[], digits = digits), print.gap = 2,
            quote = FALSE)
        cat("\n")
        if (ks) {
            cat(paste("Selection Coefficients:\n", sep = ""))
            print.default(format(cfs[], digits = digits), print.gap = 2,
                quote = FALSE)
            cat("\n")
        }
        if (kt) {
            cat(paste("Treatment Coefficients:\n", sep = ""))
            print.default(format(cft[], digits = digits), print.gap = 2,
                quote = FALSE)
            cat("\n")
        }
    }
    cat("\n")
    if (is.numeric(x$df.residual)) {
        cat(paste("df.residual: ", format(x$df.residual, digits = digits),
            "\n", sep = ""))
    }
    cat("\n")
    invisible(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
