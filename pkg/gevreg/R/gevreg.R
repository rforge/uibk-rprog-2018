gevreg <- function(formula, data, subset, na.action,
                model = TRUE, y = TRUE, x = FALSE, z = FALSE, v = FALSE,
                control = gevreg_control(...), ...)
{
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] < 3L) {
    formula <- as.Formula(formula(formula), ~ 1)
  } else {
    if(length(formula)[2L] > 3L) {
      formula <- Formula(formula(formula, rhs = 1L:3L))
      warning("formula must not have more than three RHS parts")
    }
  }
  mf$formula <- formula
  
  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  mtV <- delete.response(terms(formula, data = data, rhs = 3L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
  V <- model.matrix(mtV, mf)
  
  ## sanity check
  if(length(Y) < 1) stop("empty model")
  n <- length(Y)
  
  ## call the actual workhorse: htobit_fit()
  rval <- gevreg_fit(X, Y, Z, V, control)
  
  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(location = mtX, scale = mtZ, shape = mtV, full = mt)
  rval$levels <- list(location = .getXlevels(mtX, mf), scale = .getXlevels(mtZ, mf), shape = .getXlevels(mtV, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(location = attr(X, "contrasts"), scale = attr(Z, "contrasts"), shape = attr(V, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(location = X, scale = Z, shape = V)
  class(rval) <- "gevreg"
  return(rval)
}





gevreg_control <- function(maxit = 5000, start = NULL, ...)
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


gevreg_fit <- function(x, y, z = NULL, v = NULL, control){
  
  ## dimensions
  #n <- length(y)
  n <- nrow(y)
  if(is.null(z)) matrix(1, n, 1, dimnames = list(rownames(x), "(Intercept)"))
  if(is.null(v)) matrix(1, n, 1, dimnames = list(rownames(x), "(Intercept)"))
  m <- ncol(x)  
  p <- ncol(z)
  q <- ncol(v)
  stopifnot(n == nrow(x), n == nrow(z), n == nrow(v))
  
  nll = function(par) {
    
    # scale the coefficients back to original scaling
    par = par*parscale
    
    # define location, scale and shape coefficients
    loccoeff <- par[1:m]
    scalecoeff <- par[m + (1:p)]
    shapecoeff <- par[m + p + (1:q)]
    
    # compute the GEV parameters according to the given coefficients
    locs <- x %*% loccoeff
    scales <- z %*% scalecoeff
    shapes <- v %*% shapecoeff
    
    # compute log likelihood
    sum = 0
    for (k in 1:n) {
      s       = as.numeric(y[k,!is.na(y[k,])])
      summand = sum(dgev(s, loc = locs[k], scale = scales[k], shape = shapes[k], log = TRUE))
      sum     = sum + summand
    }
    -sum
    # ll <- sum(dgev(y, loc = locs, scale = scales, shape = shapes, log = TRUE))
    # -sum(ll)
  }
  
  # define the 'squared' log likelihood function
  nll2 = function(par) {
    
    # scale the coefficients back to original scaling
    par = par*parscale
    
    # define location, scale and shape coefficients
    loccoeff <- par[1:m]
    scalecoeff <- par[m + (1:p)]
    shapecoeff <- par[m + p + (1:q)]
    
    # compute the GEV parameters according to the given coefficients
    locs <- x %*% loccoeff
    scales <- z %*% scalecoeff
    shapes <- v %*% shapecoeff
    
    # compute the 'squared' log likelihood
    sum = 0
    for (k in 1:n) {
      s   = as.numeric(y[k,!is.na(y[k,])])
      val = 1 + shapes[k]*((s - locs[k])/scales[k])
      summand = -n*log(scales[k]) - (1 + 1/shapes[k])*sum(1/2*log(val^2)) - 
        sum(exp((-1/shapes[k])*1/2*log(val^2)))
      sum = sum + summand
    }
    -sum
  }
  
  ## starting values (by default via OLS)
  if(is.null(control$start)) {
    start.loc <- lm.fit(x, y[,1])
    start.scale <- lm.fit(z, y[,2])
    start.shape <- lm.fit(v, y[,3])
    start <- c(start.loc$coefficients,start.scale$coefficients,start.shape$coefficients)
  } else {
    start <- control$start
    stopifnot(length(start) == m + p + q)
  }
  control$start <- NULL
  
  ## rescale starting values for more stable optimization results
  parscale = 10^(floor(log10(abs(start))))
  start    = start/parscale
  
  ## if nll(start) is infinite, try to optimize the 'squared' log likelihood
  if (is.infinite(nll(start))) {
    cat("\nUsing squared log-likelihood for optimisation")
    nll.func = nll2
  } else {
    nll.func = nll
  }
  
  ## optimization
  opt <- optim(par = start, fn = nll.func, control = control, hessian = TRUE)
  #opt <- optimx(par = start, fn = nll.func, control = control, hessian = TRUE)
  
  ## enhance labels
  names(opt)[1:2] <- c("coefficients", "loglik")
  opt$coefficients <- opt$coefficients * parscale # scale back coefficients
  opt$coefficients <- list(
    location = opt$coefficients[1:m],
    scale = opt$coefficients[m + (1:p)],
    shape = opt$coefficients[m + p + (1:q)]
  )
  opt$loglik <- -opt$loglik
  opt$nobs <- n
  opt$df <- m + p + q
  
  l<- opt$value
  H<- -opt$hessian
  J<-jacobian(nll,x=unlist(opt$coefficients),method="Richardson")
  M<-(solve(H)*as.numeric(J%*%t(J)))
  opt$TIC<- -2*opt$loglik+2*sum(diag(M))
  
  return(opt)                                          
}



logLik.gevreg <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "gevreg")
}

coef.gevreg <- function(object, model = c("full", "location", "scale", "shape"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
         "location" = cf$location,
         "scale" = cf$scale,
         "shape" = cf$shape,
         "full" = c(cf$location, cf$scale, cf$shape)
  )
}


print.gevreg <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("Smooth spatial extreme value model\n\n")
  cat("Distribution: GEV\n")
  cat("Estimator: Maximum Likelihood\n")
  cat(sprintf("AIC: %f\n",AIC(x)))
  cat(sprintf("BIC: %f\n",BIC(x)))
  cat(sprintf("TIC: %f\n",x$TIC))
  
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  cat("\nCoefficients (location model):\n")
  print.default(format(x$coefficients$location, digits = digits), print.gap = 2, quote = FALSE)
  cat("\nCoefficients (scale model):\n")
  print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
  cat("\nCoefficients (shape model):\n")
  print.default(format(x$coefficients$shape, digits = digits), print.gap = 2, quote = FALSE)
  cat(sprintf("\nLog-likelihood: %s on %s Df\n", format(x$loglik, digits = digits), x$df))
  
  invisible(x)
}


terms.gevreg <- function(x, model = c("location", "scale", "shape", "full"), ...) x$terms[[match.arg(model)]]

model.frame.gevreg <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
} 

model.matrix.gevreg <- function(object, model = c("location", "scale", "shape"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
  else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(rval)
}

predict.gevreg <- function(object, newdata = NULL,
                           type = c("response", "location", "scale", "shape", "parameter", "probability", "quantile"),
                           na.action = na.pass, at = 0.5, ...)
{
  ## types of prediction
  ## response/location are synonymous
  type <- match.arg(type)
  if(type == "location") type <- "response"
  
  ## obtain model.frame/model.matrix
  tnam <- switch(type,
                 "response" = "location",
                 "scale" = "scale",
                 "shape" = "shape",
                 "full")  
  if(is.null(newdata)) {
    X <- model.matrix(object, model = "location")
    Z <- model.matrix(object, model = "scale")
    V <- model.matrix(object, model = "shape")
  } else {
    mf <- model.frame(delete.response(object$terms[[tnam]]), newdata, na.action = na.action, xlev = object$levels[[tnam]])
    if(type != "scale") X <- model.matrix(delete.response(object$terms$location), mf, contrasts = object$contrasts$location)
    if(type != "response") {
      Z <- model.matrix(object$terms$scale, mf, contrasts = object$contrasts$scale)
      V <- model.matrix(object$terms$shape, mf, contrasts = object$contrasts$shape)
    }
  }
  
  ## predicted parameters
  if(type != "scale") location <- drop(X %*% object$coefficients$location)
  if(type != "response") {
    scale <- drop(Z %*% object$coefficients$scale)
    shape <- drop(V %*% object$coefficients$shape)
  }
  
  ## compute result
  rval <- switch(type,
                 "response" = location,
                 "scale" = scale,
                 "shape" = shape,
                 "parameter" = data.frame(location, scale, shape),
                 "probability" = pgev(at, loc = location, scale = scale, shape = shape),
                 "quantile" = pmax(0, qgev(at, loc = location, scale = scale, shape = shape))
  )
  return(rval)
}
