circmax <- function(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE,
  control = circmax_control(...), ...)
{
  ## Call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## Formula
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1L:2L))
      warning("formula must not have more than two RHS parts")
    }
  }
  mf$formula <- formula

  ## Evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## Extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)

  ## Convert response to values between 0 and 2pi
  Y <- Y %% (2 * pi)

  ## Sanity check
  if(length(Y) < 1) stop("empty model")
  n <- length(Y)

  ## Call the actual workhorse: circmax_fit()
  rval <- circmax_fit(X, Y, Z, control)

  ## Further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(location = mtX, scale = mtZ, full = mt)
  rval$levels <- list(location = .getXlevels(mtX, mf), scale = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(location = attr(X, "contrasts"), scale = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(location = X, scale = Z)
  class(rval) <- "circmax"
  return(rval)
}

circmax_control <- function(maxit = 5000, start = NULL, ...)
{
  ctrl <- c(
    list(maxit = maxit,
    start = start), list(...)
  )
  if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
  ctrl$fnscale <- 1
  if(is.null(ctrl$reltol)) ctrl$reltol <- .Machine$double.eps^(1/1.2)
  ctrl
}

circmax_fit <- function(x, y, z = NULL, control)
{
  ## Dimensions
  n <- length(y)
  if(is.null(z)) matrix(1, n, 1, dimnames = list(rownames(x), "(Intercept)"))
  m <- ncol(x)  
  p <- ncol(z)
  stopifnot(n == nrow(x), n == nrow(z))

  ## Negative log-likelihood    
  nll <- function(par) {
    beta <- par[1:m]
    gamma <- par[m + (1:p)]
    mu <- x[, 1, drop = FALSE] %*% beta[1] + 2 * atan(x[, -1, drop = FALSE] %*% beta[-1])
    kappa <- exp(z %*% gamma)
    ll <- dvonmises(y, mu = mu, kappa = kappa, log = TRUE)
    -sum(ll)
  }
  
  ## Starting values (by default zeros)
  if(is.null(control$start)) {
    start <- rep(0, m + p)
  } else {
    start <- control$start
    stopifnot(length(start) == m + p)
  }
  control$start <- NULL
  
  ## Optimization
  opt <- optim(par = start, fn = nll, control = control)

  ## Collect information
  names(opt)[1:2] <- c("coefficients", "loglik")
  opt$coefficients <- list(
    location = opt$coefficients[1:m],
    scale = opt$coefficients[m + 1:p]
  )
  names(opt$coefficients$location) <- colnames(x)
  names(opt$coefficients$scale) <- colnames(z)
  opt$loglik <- -opt$loglik
  opt$nobs <- n
  opt$df <- m + p
  
  return(opt)
}

logLik.circmax <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

coef.circmax <- function(object, model = c("full", "location", "scale"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
    "location" = {
      cf$location
    },
    "scale" = {
      cf$scale
    },
    "full" = {
      structure(c(cf$location, cf$scale),
        .Names = c(names(cf$location), paste("(scale)", names(cf$scale), sep = "_")))
    }
  )
}

print.circmax <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("Maximum likelihood estimation for the von Mises distribution\n\n")
  if(x$convergence > 0) {
    cat("Model did not converge\n")
  } else {
    if(length(x$coefficients$location)) {
      cat("Coefficients (location model with tanhalf link):\n")
      print.default(format(x$coefficients$location, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in location model)\n\n")
    }
    if(length(x$coefficients$scale)) {
      cat("Coefficients (scale model (density kappa) with log link):\n")
      print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in scale model)\n\n")
    }
    cat(paste("Log-likelihood: ", format(x$loglik, digits = digits), "\n", sep = ""))
    if(length(x$df)) {
      cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
    }
    cat("\n")
  }

  invisible(x)
}

terms.circmax <- function(x, model = c("location", "scale", "full"), ...) x$terms[[match.arg(model)]]

model.frame.circmax <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
} 

model.matrix.circmax <- function(object, model = c("location", "scale"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
    else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(rval)
}

predict.circmax <- function(object, newdata = NULL,
  type = c("location", "scale", "parameter", "probability", "quantile"),
  na.action = na.pass, at = 0.5, ...)
{
  ## Types of prediction
  type <- match.arg(type)

  ## Obtain model.frame/model.matrix
  tnam <- switch(type,
    "location" = "location",
    "scale" = "scale",
    "full")  
  if(is.null(newdata)) {
    X <- model.matrix(object, model = "location")
    Z <- model.matrix(object, model = "scale")
  } else {
    mf <- model.frame(delete.response(object$terms[[tnam]]), newdata, na.action = na.action, xlev = object$levels[[tnam]])
    if(type != "scale") X <- model.matrix(delete.response(object$terms$location), mf, contrasts = object$contrasts$location)
    if(type != "location") Z <- model.matrix(object$terms$scale, mf, contrasts = object$contrasts$scale)
  }

  ## Predicted parameters
  if(type != "scale") location <- drop(X[, 1, drop = FALSE] %*% object$coefficients$location[1] 
                                       + 2 * atan(X[, -1, drop = FALSE] %*% object$coefficients$location[-1]))
  if(type != "location") scale <- exp(drop(Z %*% object$coefficients$scale))

  ## Convert location to values between 0 and 2pi
  if (type != "scale") {
    location <- location %% (2 * pi)
  }

  ## Compute result
  rval <- switch(type,
    "location" = location,
    "scale" = scale,
    "parameter" = data.frame(location, scale),
    "probability" = circular::pvonmises(at, mean = location, sd = scale),
    "quantile" = circular::qvonmises(at, mean = location, sd = scale)
  )
  return(rval)
}
