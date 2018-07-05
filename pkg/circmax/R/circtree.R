## High-level convenience interface to mob() + circfit()
circtree <- function(formula, data, na.action, 
                     mob_control = partykit::mob_control(), control = circfit_control(...), ...){

  ## Keep call
  cl <- match.call(expand.dots = TRUE)

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
  m$control <- mob_control
  m$circfit_control <- control
  #for(n in names(control)) if(!is.null(control[[n]])) m[[n]] <- control[[n]]
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
  title = "Circular distribution tree, employing the Von Mises Distribution", 
  objfun = "negative log-likelihood", ...)
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

