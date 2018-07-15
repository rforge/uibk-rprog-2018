\name{circfit_control}
\alias{circfit_control}

\title{Auxiliary Function for Controlling Circular Regression Tree Fitting,}

\description{
  Auxiliary function for \code{circfit} fitting. Specifies a list of values passed
  to \code{\link{optim}}.
}
\usage{
circfit_control(solve_kappa = solve_kappa_Newton_Fourier, ...)
}
\arguments{
  \item{solve_kappa}{Which kappa solver should be used for the starting values for kappa.
    By default a Newton Fourier is used (\code{solve_kappa_Newton_Fourier}). Alternatively,
    a uniroot provides a safe option (\code{solve_kappa_uniroot}) or code{solve_kappa_Banerjee_et_al_2005}
    provides a quick approximation.}
  \item{\dots}{additional parameters passed to \code{\link{optim}}. }
}
\value{
  A list with components named as the arguments.
}

\seealso{\code{\link{circfit}}, \code{\link{optim}}}

\keyword{distribution}
