## High-level convenience interface to mob() + circfit()
circtree <- function(formula, data, na.action, cor = FALSE, ...)
{
  ## Keep call
  cl <- match.call(expand.dots = TRUE)

  ## Use dots for setting up mob_control
  control <- partykit::mob_control(...)
  control$ytype <- "matrix"

  ## Control options for circfit
  circcontrol <- list(cor = cor)

  ## Formula
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula::Formula(formula(formula, rhs = 2L))
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }

  ## Call mob
  m <- match.call(expand.dots = FALSE)
  m$fit <- dist_family_fit
  #m$fit <- circfi
  m$formula <- formula
  m$control <- control
  for(n in names(circcontrol)) if(!is.null(circcontrol[[n]])) m[[n]] <- circcontrol[[n]]
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.call(quote(partykit::mob))
  rval <- eval(m, parent.frame())
  
  ## Extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("circtree", class(rval))
  return(rval)
}

## Methods
print.circtree <- function(x,
  title = "CIRCULAR tree", objfun = "negative log-likelihood", ...)
{
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}
