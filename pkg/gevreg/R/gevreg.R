gevreg <- function (formula, data, subset, na.action, model = TRUE, y = TRUE, 
                    x = TRUE, z = FALSE, v = FALSE, gev_params, control = gevreg_control(...), 
                    ...) 
{
  cl <- match.call()
  if (missing(data)) 
    data <- environment(formula)
  if(!missing(gev_params)){
    if(!is.data.frame(gev_params)) stop("gev_params must be data frame")
    if (!all(names(gev_params) %in% c("loc","scale","shape"))) stop("columns in data farme gev_params must be named: loc, scale and shape")
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if (length(formula)[2L] == 1L) {
      formula <- as.Formula(formula(formula), ~1 | 1)
  } else if (length(formula)[2L] == 2L){
    formula <- as.Formula(formula(formula), ~1)
  } else {
    if (length(formula)[2L] > 3L) {
      formula <- Formula(formula(formula, rhs = 1L:3L))
      warning("formula must not have more than three RHS parts")
    }
  }
  
  n.stats <- length(levels(data$station))
  n.years <- nrow(data)/n.stats
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  mtV <- delete.response(terms(formula, data = data, rhs = 3L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
  V <- model.matrix(mtV, mf)
  if (length(Y) < 1) 
    stop("empty model")
  n <- length(Y)
  rval <- gevreg_fit(X, Y, Z, V, n.stats, n.years, gev_params, control)
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(location = mtX, scale = mtZ, shape = mtV, 
                     full = mt)
  rval$levels <- list(location = .getXlevels(mtX, mf), scale = .getXlevels(mtZ, 
                                                                           mf), shape = .getXlevels(mtV, mf), full = .getXlevels(mt, 
                                                                                                                                 mf))
  rval$contrasts <- list(location = attr(X, "contrasts"), scale = attr(Z, 
                                                                       "contrasts"), shape = attr(V, "contrasts"))
  if (model) 
    rval$model <- mf
  if (y) 
    rval$y <- Y
  if (x) 
    rval$x <- list(location = X, scale = Z, shape = V)
  class(rval) <- "gevreg"
  return(rval)
}



# gevreg_control <- function(maxit = 5000, start = NULL, grad = TRUE, hessian = FALSE, ...)
# {
#   if(is.logical(hessian)) hessian <- if(hessian) "optim" else "none"
#   if(is.character(hessian)) hessian <- match.arg(tolower(hessian), c("numderiv", "optim", "none"))
#   ctrl <- c(
#     list(maxit = maxit, start = start, grad = grad, hessian = hessian),
#     list(...)
#   )
#   if(is.null(ctrl$method)) {
#     ctrl$method <- if(grad) "nlminb" else "Nelder-Mead"
#   }
#   if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
#   ctrl$fnscale <- 1
#   if(is.null(ctrl$reltol)) ctrl$reltol <- .Machine$double.eps^(1/1.2)
#   ctrl
# }
gevreg_control <- function(maxit = 5000, start = NULL, grad = TRUE, ...)
{
  #if(is.logical(hessian)) hessian <- if(hessian) "numderiv" else "none"
  #if(is.character(hessian)) hessian <- match.arg(tolower(hessian), c("numderiv"))
  ctrl <- c(
    list(maxit = maxit, start = start, grad = grad),
    list(...)
  )
  if(is.null(ctrl$method)) {
    ctrl$method <- if(grad) "nlminb" else "Nelder-Mead"
  }
  #if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
  #ctrl$fnscale <- 1
  #if(is.null(ctrl$reltol)) ctrl$reltol <- .Machine$double.eps^(1/1.2)
  ctrl
}


gevreg_fit <- function (x, y, z = NULL, v = NULL, n.stats, n.years, gevp, control) 
{
  n <- length(y)
  if (is.null(z)) matrix(1, n, 1, dimnames = list(rownames(x), "(Intercept)"))
  if (is.null(v)) matrix(1, n, 1, dimnames = list(rownames(x), "(Intercept)"))
  m <- ncol(x)
  p <- ncol(z)
  q <- ncol(v)
  stopifnot(n == nrow(x), n == nrow(z), n == nrow(v))
  stopifnot(n.stats == nrow(gevp))
  
  nll <- function(par) {
    par <- par * parscale
    loccoeff <- par[1:m]
    scalecoeff <- par[m + (1:p)]
    shapecoeff <- par[m + p + (1:q)]
    locs <- x %*% loccoeff
    scales <- exp(z %*% scalecoeff)
    shapes <- v %*% shapecoeff
    i <- which(is.na(y))
    -sum(dgev.gevreg(na.omit(y), loc = locs[-i], scale = scales[-i], shape = shapes[-i], log = TRUE))
  }
  
  nll2 <- function(par) {
    par <- par * parscale
    loccoeff <- par[1:m]
    scalecoeff <- par[m + (1:p)]
    shapecoeff <- par[m + p + (1:q)]
    locs <- x %*% loccoeff
    scales <- exp(z %*% scalecoeff)
    shapes <- v %*% shapecoeff
    sum = 0
    for (k in 1:n.stats) {
      #z = as.numeric(max_data[k, !is.na(max_data[k, ])])
      idx <- (k-1)*n.years + 1:n.years
      z <- na.omit(y[idx])
      iloc <- locs[idx][1]
      iscale <- scales[idx][1]
      ishape <- shapes[idx][1]
      val = 1 + ishape * ((z - iloc)/iscale)
      n.obs = length(z)
      summand = -n.obs * log(iscale) - (1 + 1/iscale) * 
        sum(1/2 * log(val^2)) - sum(exp((-1/ishape) * 
                                          1/2 * log(val^2)))
      if(is.na(summand))  print(iscale)
      sum = sum + summand[1]
    }
    -sum
  }
  
  ngr <- function(par) {
    return(grad(nll, par))
  }
  
  grad <- control$grad
  hess <- control$hessian
  meth <- control$method
  control$grad <- control$hessian <- control$method <- NULL
  if (is.null(control$start)) {
    ploc <- pscale <- pshape <- c()
    for (k in 1:n.stats) {
      idx <- (k-1)*n.years + 1:n.years
      ploc[idx] <- gevp$loc[k]
      pscale[idx] <- gevp$scale[k]
      pshape[idx] <- gevp$shape[k]
    }
    start.loc <- glm.fit(x, ploc)
    start.scale <- glm.fit(z, pscale)
    start.shape <- glm.fit(v, pshape)
    start <- c(start.loc$coefficients,start.scale$coefficients,start.shape$coefficients)
  }
  else {
    start <- control$start
    stopifnot(length(start) == m + p + q)
  }
  control$start <- NULL
  parscale = 10^(floor(log10(abs(start))))
  start = start/parscale
  
  if (is.infinite(nll(start))) {
    optloglik = nll2
  } else {
    optloglik = nll
  }
  
  opt <- optimx(start, optloglik, method = meth, 
                itnmax = NULL, control = list(follow.on = FALSE, kkt = FALSE))
  
  # if (hess == "none") {
  #   opt <- c(opt, list(hessian = NULL))
  # } else if (hess == "numderiv") {
  #   opt$hessian <- numDeriv::hessian(nll, opt$par)
  # }
  if (!is.null(opt$hessian)) {
    rownames(opt$hessian) <- colnames(opt$hessian) <- c(colnames(x),
                                                        paste("(scale)", colnames(z), sep = "_"), paste("(shape)",
                                                                                                        colnames(v), sep = "_"))
    opt$vcov <- solve(opt$hessian)
  }
  
  res <- list()
  res$coefficients <- as.numeric(opt[1:(m + p + q)]) * parscale
  res$coefficients <- list(location = res$coefficients[1:m], 
                           scale = res$coefficients[m + (1:p)], 
                           shape = res$coefficients[m + p + (1:q)],
                           all = res$coefficients[1:(m + p + q)])
  mu <- drop(x %*% res$coefficients$location)
  sigma <- exp(drop(z %*% res$coefficients$scale))
  xi <- drop(v %*% res$coefficients$shape)
  res$residuals <- y - mu
  res$fitted.values <- list(location = mu, scale = sigma, shape = xi)
  res$method <- meth
  res$loglik <- -opt$value
  res$nobs <- n
  res$df <- m + p + q
  res$convergence <- opt$convcode
  res$niter <- opt$niter
  res$feval <- opt$fevals
  res$geval <- opt$gevals
  res$xtimes <- opt$xtimes
  res$hessian <- opt$hessian
  res$vcov <- opt$vcov
  return(res)
}

logLik.gevreg <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
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


## nams of covriats from gr$terms$location...
print.gevreg <- function(x, digits = max(3, getOption("digits") - 3), ...){
    cat("Smooth Spatial GEV Fitting\n\n")
    cat("Estimator: Maximum Likelihood\n")
    cat(sprintf("AIC: %f\n",AIC(x)))
    cat(sprintf("BIC: %f\n",BIC(x)))
     
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
    
    cat("\nConvergence:", ifelse(x$convergence == 0, paste("successful, method:",x$method,"\n"), 
                                 paste("not successful: code = ",x$convergence, "with",
                                       ifelse(is.na(x$niter),"",paste(x$niter, "iterations")),  
                                       ifelse(is.na(x$feval),"",paste(x$feval, "function evaluations")), 
                                       ifelse(is.na(x$geval),"",paste(x$geval, "gradient evaluations.\n")), "", sep = " ")))
    
    cat("\nCoefficients (location model):\n")
    locs <- x$coefficients$location
    print.default(format(locs, digits = digits), print.gap = 2, quote = FALSE)
    cat("\nCoefficients (scale model):\n")
    scales <- x$coefficients$scale
    print.default(format(scales, digits = digits), print.gap = 2, quote = FALSE)
    cat("\nCoefficients (shape model):\n")
    shapes <- x$coefficients$shape
    print.default(format(shapes, digits = digits), print.gap = 2, quote = FALSE)
    cat(sprintf("\nLog-likelihood: %s on %s Df\n", format(x$loglik, digits = digits), x$df))
  invisible(x)
}

# vcov.gevreg <- function(object, model = c("full", "location", "scale", "shape"), ...)
# {
#   vc <- object$vcov
#   k <- length(object$coefficients$location)
#   m <- length(object$coefficients$scale)
#   n <- length(object$coefficients$shape)
#   model <-  match.arg(model)
#   switch(model,
#          "location" = {
#            vc <- vc[seq.int(length.out = k), seq.int(length.out = k), drop = FALSE]
#            #colnames(vc) <- rownames(vc) <- names(object$coefficients$location)
#            colnames(vc) <- rownames(vc) <- paste0("(location)_station",1:3)
#            vc
#          },
#          "scale" = {
#            vc <- vc[seq.int(length.out = m) + k, seq.int(length.out = m) + k, drop = FALSE]
#            vc
#          },
#          "shape" = {
#            vc <- vc[seq.int(length.out = n) + k, seq.int(length.out = n) + k, drop = FALSE]
#            #colnames(vc) <- rownames(vc) <- names(object$coefficients$shape)
#            #colnames(vc) <- rownames(vc) <- paste0("(scale)_station",1:3)
#            vc
#          },
#          "full" = {
#            vc
#          }
#   )
# }

dgev.gevreg <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE)
{
  if (min(scale) <= 0)
    stop("invalid scale parameter")
  x <- (x - loc)/scale
  n <- length(x)
  loc <- rep_len(loc, n)
  shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)
  xs <- 1 + shape * x
  pos <- xs > 0 | is.na(xs)
  s0 <- abs(shape) < .Machine$double.eps^0.5
  dns <- rep_len(-Inf, n)
  dns[pos] <- -log(scale[pos]) - xs[pos]^(-1/shape[pos]) - (1/shape[pos] + 1) * log(xs[pos])
  dns[s0] <- -log(scale[s0]) - x[s0] - exp(-x[s0])
  if (!log) dns <- exp(dns)
  return(dns)
}

pgev.gevreg <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE) 
{
  if (min(scale) <= 0) 
    stop("invalid scale parameter")
  q <- (q - loc)/scale
  n <- length(loc)
  
  s0 <- abs(shape) < .Machine$double.eps^0.5
  pos <- abs(shape) > 0 
  p <- rep_len(NA, n)
  p[s0] <- p <- exp(-exp(-q))
  p[pos] <- exp(-pmax(1 + shape * q, 0)^(-1/shape))
  if (!lower.tail) 
    p <- 1 - p
  return(p)
}

qgev.gevreg <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE) 
{
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1) 
    stop("`p' must contain probabilities in (0,1)")
  if (min(scale) < 0) 
    stop("invalid scale parameter")
  n <- length(loc)
  if (!lower.tail) 
    p <- 1 - p
  s0 <- abs(shape) < .Machine$double.eps^0.5
  pos <- abs(shape) > 0
  q <- rep_len(NA, n)
  q[s0] <- loc - scale * log(-log(p))
  q[pos] <- loc + scale * ((-log(p))^(-shape) - 1)/shape
  #if (shape == 0) 
  #  return(loc - scale * log(-log(p)))
  #else return(loc + scale * ((-log(p))^(-shape) - 1)/shape)
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
    scale <- exp( drop(Z %*% object$coefficients$scale) )
    shape <- drop(V %*% object$coefficients$shape)
  }
  
  deconstruct <- function(x,at){
    if(length(at) == 1) return(x)
    n <- length(x)/length(at)
    res <- matrix(NA,n,length(at))
    for(i in 1:length(at)){
      idx <- (i-1)*n + 1:n
      res[,i] <- x[idx]
    }
    res <- data.frame(res)
    colnames(res) <- paste0(type,at)
    res
  }

  ## compute result
  rval <- switch(type,
                 "response" = unique(location),
                 "scale" = unique(scale),
                 "shape" = unique(shape),
                 "parameter" = unique(data.frame(location, scale, shape)),
                 "probability" = deconstruct(unique(pgev.gevreg(at, loc = location, scale = scale, shape = shape)),at),
                 "quantile" = deconstruct(unique(pmax(0, qgev.gevreg(at, loc = location, scale = scale, shape = shape))),at)
  )
  return(rval)
}


residuals.gevreg <- function(object, type = c("standardized", "pearson", "response"), ...) {
  if(match.arg(type) == "response") {
    object$residuals
  } else {
    object$residuals/object$fitted.values$scale
  }
}
#
#print.summary.gevreg

# summary.gevreg <- function (object, ...) 
# {
#   object$residuals <- object$residuals/object$fitted.values$scale
#   k <- length(object$coefficients$location)
#   m <- length(object$coefficients$scale)
#   n <- length(object$coefficients$shape)
#   cf <- as.vector(do.call("c", object$coefficients))
#   se <- sqrt(diag(object$vcov))
#   cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
#   colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
#   cf <- list(location = cf[seq.int(length.out = k), , drop = FALSE], 
#              scale = cf[seq.int(length.out = m) + k, , drop = FALSE],
#              shape = cf[seq.int(length.out = n) + m + k, , drop = FALSE])
#   rownames(cf$location) <- names(object$coefficients$location)
#   rownames(cf$scale) <- names(object$coefficients$scale)
#   rownames(cf$shape) <- names(object$coefficients$shape)
#   object$coefficients <- cf
#   object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL
#   class(object) <- "summary.gevreg"
#   object
# }
# 
# print.summary.gevreg <- function (x, digits = max(3, getOption("digits") - 3), ...) 
# {
#   cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 
#                                                         0.85)), "", sep = "\n")
#   if (x$convergence > 0L) {
#     cat("model did not converge\n")
#   }
#   else {
#     cat(paste("Standardized residuals:\n", sep = ""))
#     print(structure(round(as.vector(quantile(x$residuals)), 
#                           digits = digits), .Names = c("Min", "1Q", "Median", 
#                                                        "3Q", "Max")))
#     if (NROW(x$coefficients$location)) {
#       cat(paste("\nCoefficients (location model):\n", sep = ""))
#       printCoefmat(x$coefficients$location, digits = digits, 
#                    signif.legend = FALSE)
#     } else cat("\nNo coefficients (in location model)\n")
#     if (NROW(x$coefficients$scale)) {
#       cat(paste("\nCoefficients (scale model with log link):\n", 
#                 sep = ""))
#       printCoefmat(x$coefficients$scale, digits = digits, 
#                    signif.legend = FALSE)
#     } else cat("\nNo coefficients ( in scale model)\n")
#     if (NROW(x$coefficients$shape)) {
#       cat(paste("\nCoefficients (shape model):\n", 
#                 sep = ""))
#       printCoefmat(x$coefficients$shape, digits = digits, 
#                    signif.legend = FALSE)
#     } else cat("\nNo coefficients ( in shape model)\n")
#     if (getOption("show.signif.stars") & any(do.call("rbind", 
#                                                      x$coefficients)[, 4L] < 0.1, na.rm = TRUE)) 
#       cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", 
#           "\n")
#     cat("\nLog-likelihood:", formatC(x$loglik, digits = digits), 
#         "on", sum(sapply(x$coefficients, NROW)), "Df\n")
#     cat(paste("Number of iterations in", x$method, "optimization:", 
#               x$count[2L], "\n"))
#   }
#   invisible(x)
# }
