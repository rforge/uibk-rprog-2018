gevreg <- function (formula, data, subset, na.action, model = TRUE, y = TRUE, 
                    x = FALSE, z = FALSE, v = FALSE, control = gevreg_control(...), 
                    ...) 
{
  cl <- match.call()
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 
             0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  a <- attr(terms(formula), "term.labels")
  if (length(formula)[2L] == 1L) {
    if (length(a) == 0) {
      formula <- as.Formula(formula(formula), ~1 | 1)
    }
    else if (a[1] == "station") {
      formula <- update(formula(formula), . ~ . - 1)
      formula <- as.Formula(formula(formula), ~-1 + station | 
                              -1 + station)
    }
    else {
      formula <- update(formula(formula), . ~ . - 1)
      formula <- as.Formula(formula(formula), ~1 | 1)
    }
  }
  else if (length(formula)[2L] == 2L) {
    formula.l <- formula(formula, rhs = 0)
    formula <- eval(parse(text = paste0(formula.l, "-1 + ", 
                                        a)))
    formula <- as.Formula(formula(formula), ~1)
  }
  else {
    if (length(formula)[2L] > 3L) {
      formula <- Formula(formula(formula, rhs = 1L:3L))
      warning("formula must not have more than three RHS parts")
    }
  }
  mf$formula <- formula
  if (!is.factor(data$station)) 
    stop("station column must be a factor")
  n.stats <- length(levels(data$station))
  if (m[3] != 0) {
    data$station <- as.numeric(data$station)
    n.stats <- 1
  }
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
  rval <- gevreg_fit(X, Y, Z, V, n.stats, control)
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


gevreg_control <- function(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, ...)
{
  if(is.logical(hessian)) hessian <- if(hessian) "optim" else "none"
  if(is.character(hessian)) hessian <- match.arg(tolower(hessian), c("numderiv", "optim", "none"))
  ctrl <- c(
    list(maxit = maxit, start = start, grad = grad, hessian = hessian),
    list(...)
  )
  if(is.null(ctrl$method)) {
    ctrl$method <- if(grad) "BFGS" else "Nelder-Mead"
  }
  if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
  ctrl$fnscale <- 1
  if(is.null(ctrl$reltol)) ctrl$reltol <- .Machine$double.eps^(1/1.2)
  ctrl
}


gevreg_fit <- function (x, y, z = NULL, v = NULL, n.stats, control) 
{
  n <- length(y)
  if (is.null(z)) 
    matrix(1, n, 1, dimnames = list(rownames(x), "(Intercept)"))
  if (is.null(v)) 
    matrix(1, n, 1, dimnames = list(rownames(x), "(Intercept)"))
  m <- ncol(x)
  p <- ncol(z)
  q <- ncol(v)
  stopifnot(n == nrow(x), n == nrow(z), n == nrow(v))
  nll <- function(par, dat = NULL) {
    loccoeff <- par[1:m]
    scalecoeff <- par[m + (1:p)]
    shapecoeff <- par[m + p + (1:q)]
    locs <- x %*% loccoeff
    scales <- exp(z %*% scalecoeff)
    shapes <- v %*% shapecoeff
    -sum(dgev.gevreg(y, loc = locs, scale = scales, shape = shapes, 
                     log = TRUE))
  }
  ngr <- function(par, dat) {
    return(grad(nll, par))
    # loccoeff <- par[1:m]
    # scalecoeff <- par[m + (1:p)]
    # shapecoeff <- par[m + p + (1:q)]
    # mu <- x %*% loccoeff
    # sigma <- exp(z %*% scalecoeff)
    # xi <- v %*% shapecoeff
    # rval <- matrix(0, nrow = nrow(x), ncol = ncol(x) + ncol(z) + 
    #                  ncol(v))
    # if (!isTRUE(all.equal(xi, 0))) {
    #   rval <- c(sum(-(-(mu - dat) * xi/sigma + 1)^(-1/xi - 
    #                                                  1)/sigma - xi * (1/xi + 1)/(sigma * ((mu - dat) * 
    #                                                                                         xi/sigma - 1))), sum(-(dat - mu) * ((dat - mu) * 
    #                                                                                                                               xi/sigma + 1)^(-1/xi - 1)/sigma^2 + (dat - mu) * 
    #                                                                                                                (xi + 1)/(sigma^2 * ((dat - mu) * xi/sigma + 
    #                                                                                                                                       1)) - 1/sigma), sum(-(mu - dat) * (1/xi + 1)/(sigma * 
    #                                                                                                                                                                                       ((mu - dat) * xi/sigma - 1)) - (log(-(mu - dat) * 
    #                                                                                                                                                                                                                             xi/sigma + 1)/xi^2 - (mu - dat)/(sigma * ((mu - 
    #                                                                                                                                                                                                                                                                          dat) * xi/sigma - 1) * xi))/(-(mu - dat) * xi/sigma + 
    #                                                                                                                                                                                                                                                                                                         1)^(1/xi) + log(-(mu - dat) * xi/sigma + 1)/xi^2))
    # }
    # else {
    #   rval <- c(sum(-exp(mu/sigma - dat/sigma)/sigma + 
    #                   1/sigma), sum(mu * exp(mu/sigma - dat/sigma)/sigma^2 - 
    #                                   dat * exp(mu/sigma - dat/sigma)/sigma^2 - mu/sigma^2 - 
    #                                   1/sigma + dat/sigma^2), 0)
    # }
    # return(rval)
  }
  grad <- control$grad
  hess <- control$hessian
  meth <- control$method
  control$grad <- control$hessian <- control$method <- NULL
  if (is.null(control$start)) {
    start.scale <- sqrt(6 * var(y, na.rm = TRUE))/pi
    start.loc <- mean(y, na.rm = TRUE) - 0.58 * start.scale
    start.shape <- 0
    start <- c(rep(start.loc, n.stats), log(rep(start.scale, 
                                                n.stats)), rep(start.shape, n.stats))
  }
  else {
    start <- control$start
    stopifnot(length(start) == m + p + q)
  }
  control$start <- NULL
  opt <- if (grad) {
    optim(par = start, fn = nll, gr = ngr, dat = y, control = control, 
          method = meth, hessian = (hess == "optim"))
  }
  else {
    optim(par = start, fn = nll, control = control, method = meth, 
          hessian = (hess == "optim"))
  }
  if (hess == "none") {
    opt <- c(opt, list(hessian = NULL))
  }
  else if (hess == "numderiv") {
    opt$hessian <- numDeriv::hessian(nll, opt$par)
  }
  if (!is.null(opt$hessian)) {
    rownames(opt$hessian) <- colnames(opt$hessian) <- c(colnames(x), 
                                                        paste("(scale)", colnames(z), sep = "_"), paste("(shape)", 
                                                                                                        colnames(v), sep = "_"))
    opt$vcov <- solve(opt$hessian)
  }
  names(opt)[1:2] <- c("coefficients", "loglik")
  opt$coefficients <- opt$coefficients
  opt$coefficients <- list(location = opt$coefficients[1:m], 
                           scale = opt$coefficients[m + (1:p)], shape = opt$coefficients[m + 
                                                                                           p + (1:q)])
  mu <- drop(x %*% opt$coefficients$location)
  sigma <- exp(drop(z %*% opt$coefficients$scale))
  xi <- drop(v %*% opt$coefficients$shape)
  opt$residuals <- y - mu
  opt$fitted.values <- list(location = mu, scale = sigma, shape = xi)
  opt$method <- meth
  opt$loglik <- -opt$loglik
  opt$nobs <- n
  opt$df <- m + p + q
  return(opt)
}



logLik.gevreg <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

coef.gevreg <- function(object, model = c("full", "location", "scale", "shape"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
         "location" = cf$location,
         "scale" = exp(cf$scale),
         "shape" = cf$shape,
         "full" = c(cf$location, exp(cf$scale), cf$shape)
  )
}


print.gevreg <- function(x, digits = max(3, getOption("digits") - 3), ext = FALSE, ...){
  if (!ext){
    cf <- cbind(loc=x$coefficients$location,
            scale=exp(x$coefficients$scale),
            shape=x$coefficients$shape)
    rownames(cf) <- paste0("station",1:nrow(cf))
    print.default(format(cf, digits = digits, scientific = FALSE), print.gap = 2, quote = FALSE)
  } else {
    cat("Maximum Likelihood GEV Fitting\n\n")
    cat("Estimator: Maximum Likelihood\n")
    cat(sprintf("AIC: %f\n",AIC(x)))
    cat(sprintf("BIC: %f\n",BIC(x)))
    #cat(sprintf("TIC: %f\n",x$TIC))
    
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
    
    cat("\nConvergence:", ifelse(x$convergence == 0, "successful", paste("not successful: code = ",x$convergence, "with", x$counts, "function evaluations")), "", sep = "\n")
    
    cat("\nCoefficients (location model):\n")
    locs <- x$coefficients$location; names(locs) <- paste0("station_",1:length(locs))
    print.default(format(locs, digits = digits), print.gap = 2, quote = FALSE)
    cat("\nCoefficients (scale model with log link):\n")
    scales <- x$coefficients$scale; names(scales) <- paste0("station_",1:length(scales))
    print.default(format(exp(scales), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nCoefficients (shape model):\n")
    shapes <- x$coefficients$shape; names(shapes) <- paste0("station_",1:length(shapes))
    print.default(format(shapes, digits = digits), print.gap = 2, quote = FALSE)
    cat(sprintf("\nLog-likelihood: %s on %s Df\n", format(x$loglik, digits = digits), x$df))
  }
  invisible(x)
}

vcov.gevreg <- function(object, model = c("full", "location", "scale", "shape"), ...)
{
  vc <- object$vcov
  k <- length(object$coefficients$location)
  m <- length(object$coefficients$scale)
  n <- length(object$coefficients$shape)
  model <-  match.arg(model)
  switch(model,
         "location" = {
           vc <- vc[seq.int(length.out = k), seq.int(length.out = k), drop = FALSE]
           #colnames(vc) <- rownames(vc) <- names(object$coefficients$location)
           colnames(vc) <- rownames(vc) <- paste0("(location)_station",1:3)
           vc
         },
         "scale" = {
           vc <- vc[seq.int(length.out = m) + k, seq.int(length.out = m) + k, drop = FALSE]
           vc
         },
         "shape" = {
           vc <- vc[seq.int(length.out = n) + k, seq.int(length.out = n) + k, drop = FALSE]
           #colnames(vc) <- rownames(vc) <- names(object$coefficients$shape)
           #colnames(vc) <- rownames(vc) <- paste0("(scale)_station",1:3)
           vc
         },
         "full" = {
           vc
         }
  )
}

## need pgev and qgev as vector version as well
dgev.gevreg <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE) 
{
  #print(scale)
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

pgev.gevreg <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE){}
qgev.gevreg <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE){}
  


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

  ## compute result
  rval <- switch(type,
                 "response" = location,
                 "scale" = scale,
                 "shape" = shape,
                 "parameter" = data.frame(location, scale, shape),
                 "probability" = pgev.gevreg(at, loc = location, scale = scale, shape = shape),
                 "quantile" = pmax(0, qgev.gevreg(at, loc = location, scale = scale, shape = shape))
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

summary.gevreg <- function (object, ...) 
{
  object$residuals <- object$residuals/object$fitted.values$scale
  k <- length(object$coefficients$location)
  m <- length(object$coefficients$scale)
  n <- length(object$coefficients$shape)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(location = cf[seq.int(length.out = k), , drop = FALSE], 
             scale = cf[seq.int(length.out = m) + k, , drop = FALSE],
             shape = cf[seq.int(length.out = n) + m + k, , drop = FALSE])
  rownames(cf$location) <- names(object$coefficients$location)
  rownames(cf$scale) <- names(object$coefficients$scale)
  rownames(cf$shape) <- names(object$coefficients$shape)
  object$coefficients <- cf
  object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL
  class(object) <- "summary.gevreg"
  object
}

print.summary.gevreg <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 
                                                        0.85)), "", sep = "\n")
  if (x$convergence > 0L) {
    cat("model did not converge\n")
  }
  else {
    cat(paste("Standardized residuals:\n", sep = ""))
    print(structure(round(as.vector(quantile(x$residuals)), 
                          digits = digits), .Names = c("Min", "1Q", "Median", 
                                                       "3Q", "Max")))
    if (NROW(x$coefficients$location)) {
      cat(paste("\nCoefficients (location model):\n", sep = ""))
      printCoefmat(x$coefficients$location, digits = digits, 
                   signif.legend = FALSE)
    } else cat("\nNo coefficients (in location model)\n")
    if (NROW(x$coefficients$scale)) {
      cat(paste("\nCoefficients (scale model with log link):\n", 
                sep = ""))
      printCoefmat(x$coefficients$scale, digits = digits, 
                   signif.legend = FALSE)
    } else cat("\nNo coefficients ( in scale model)\n")
    if (NROW(x$coefficients$shape)) {
      cat(paste("\nCoefficients (shape model):\n", 
                sep = ""))
      printCoefmat(x$coefficients$shape, digits = digits, 
                   signif.legend = FALSE)
    } else cat("\nNo coefficients ( in shape model)\n")
    if (getOption("show.signif.stars") & any(do.call("rbind", 
                                                     x$coefficients)[, 4L] < 0.1, na.rm = TRUE)) 
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", 
          "\n")
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits), 
        "on", sum(sapply(x$coefficients, NROW)), "Df\n")
    cat(paste("Number of iterations in", x$method, "optimization:", 
              x$count[2L], "\n"))
  }
  invisible(x)
}
