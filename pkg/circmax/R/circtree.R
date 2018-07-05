## High-level convenience interface to mob() + circfit()
circtree <- function(formula, data, na.action, control = partykit::mob_control(), ...){

  ## Keep call
  cl <- match.call(expand.dots = TRUE)

  ## Use dots for setting up circfit_control
  circcontrol <- circfit_control(...)

  ## Formula
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula::Formula(formula(formula, rhs = 2L))
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }

  ## Call mob
  m <- match.call(expand.dots = FALSE)
  #m$fit <- dist_family_fit
  m$fit <- circfit
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


## control function for circfit
circfit_control <- function(solve_kappa = solve_kappa_Newton_Fourier, ...) {
  ctrl <- c(
    list(solve_kappa = solve_kappa), 
    list(...)
  )
  ctrl
}

