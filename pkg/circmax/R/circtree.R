## High-level convenience interface to mob() + circfit()
circtree <- function(formula, data, start, subset, na.action, weights, offset,
                     control = partykit::mob_control(), fit_control = circfit_control(...), ...){

  ## Keep call
  cl <- match.call(expand.dots = TRUE)

  ## Data
  if(missing(data)) data <- environment(formula)

  ## Formula
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula::Formula(formula(formula, rhs = 2L))
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }

  ## Call mob (with reordered arguments)
  m <- match.call(expand.dots = FALSE)
  morder <- match(c("formula", "data", "start", "subset", "na.action", "weights", "offset"), names(m), 0L)
  m <- m[c(1L, morder)]
  m$fit <- circfit
  m$formula <- formula
  m$control <- control
  m$fit_control <- fit_control
  #for(n in names(control)) if(!is.null(control[[n]])) m[[n]] <- control[[n]]
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.call(quote(partykit::mob))
  rval <- eval(m, parent.frame())
  
  ## Extend class and keep original call
  rval$info$call <- cl
  rval$info$fit_control <- fit_control
  class(rval) <- c("circtree", class(rval))
  return(rval)
}


## Control function for circfit
circfit_control <- function(solve_kappa = solve_kappa_Newton_Fourier, ...) {
  ctrl <- c(
    list(solve_kappa = solve_kappa), 
    list(...)
  )
  ctrl
}


## Print method
print.circtree <- function(x, title = "Circular distribution tree, employing the Von Mises Distribution", 
                           objfun = "negative log-likelihood", ...) {
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}


## Predict method
predict.circtree <- function (object, newdata = NULL, type = c("parameter", "node"), OOB = FALSE, ...) {
  # FIXME: Implement prediciton for response, currently not working probably to not supported methods for circfit

  # set default to 'parameter'
  if(length(type)>1) type <- type[1]

  if(type == "node") {
      partykit::predict.modelparty(object = object, newdata = newdata, type = type, OOB = OOB, ...)
  } else if(type == "parameter") {
    pred.subgroup <- partykit::predict.modelparty(object, newdata =  newdata, type = "node")
    groupcoef <- coef(object)
    if(is.vector(groupcoef)) {
      groupcoef <- t(as.data.frame(groupcoef))
      rownames(groupcoef) <- 1
    }
    pred.par <- groupcoef[paste(pred.subgroup),]
    rownames(pred.par) <- c(1: (NROW(pred.par)))
    pred.par <- as.data.frame(pred.par)
    return(pred.par)
  }
}

## Coef method
coef.circtree <- function(object, ...){
  partykit:::coef.modelparty(object)
}


## Loglik method
logLik.circtree <- function(object) {
  tmp <- partykit::nodeapply(object, ids = partykit::nodeids(object, terminal = TRUE), 
    FUN = function(n) partykit::info_node(n)$objfun)
  nll <- sum(unlist(tmp))
  return(structure(-nll, df = ncol(coef(object)) * width(object) + width(object) - 1 , class = "logLik"))
}

