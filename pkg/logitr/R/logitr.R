#### new logit try

logitr <- function(formula, data, subset, na.action,
                   model = TRUE, x = FALSE, y = TRUE,
                   control = logitr_control(...), ...)
{
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ##
  mf$formula <- as.formula(formula)

  ## model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## terms, model matrix, response
  mt <- terms(formula, data = data)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf) # maybe add contrasts later

  ##
  if(length(Y) < 1) stop("empty model")
  n <- length(Y)

  ## call workhorse function
  rval <- logitr_fit(X, Y, control)

  ## info
  rval$call <- cl
  rval$formula <- formula
  rval$terms <- mt
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- X
  class(rval) <- "logitr"
  return(rval)
}

logitr_control <- function(maxit = 5000, start = NULL, hessian = TRUE, ...)
{
  if(is.logical(hessian)) hessian <- if(hessian) "numderiv" else "none"
  if(is.character(hessian)) hessian <- match.arg(tolower(hessian), c("numderiv", "optim", "none"))
  ctrl <- c(list(maxit = maxit,start = start),
            hessian = hessian, list(...))
  if(is.null(ctrl$method)) ctrl$method <- "BFGS"
  ctrl
}

logitr_fit <- function(x, y, control)
{
  ## dim
  if(any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
  n <- length(y)
  m <- ncol(x)
  l <- nrow(x)
  if(is.null(l)) l <- length(x)

  ##negative loglik
  nll <- function(par){
    eta <- x%*%par
    pi <- exp(eta)*((1+exp(eta))^-1)
    #ll <- dbinom(y, size = 1, prob = pi, log = TRUE) ## y * log(pi) + (1 - y) * log(1 - pi)
    ll <- y * log(pi) + (1 - y) * log(1 - pi)
    -sum(ll)
  }

  ##
  hess <- control$hessian
  meth <- control$method
  control$hessian <- control$method <- NULL

  ## starting values (by default via OLS)
  if(is.null(control$start)){
    #start <- rep(0,m)
    start <- lm.fit(x, y) ## OLS
    start <- start$coefficients
  } else {
    stopifnot(length(start) == m)
  }
  control$start <- NULL

  ## optimization v1
  opt <- optim(par = start, fn = nll, control = control, method = meth, hessian = (hess == "optim"))

  ## hessian
  if(hess == "none") {
    opt <- c(opt, list(hessian = NULL))
  } else if(hess == "numderiv") {
    opt$hessian <- numDeriv::hessian(nll, opt$par)
  }
  if(!is.null(opt$hessian)) {
    rownames(opt$hessian) <- colnames(opt$hessian) <- colnames(x)
    opt$vcov <- solve(opt$hessian)
    opt$hessian <- NULL
  }


  ##
  names(opt)[1:2] <- c("coefficients", "loglik")
  names(opt$coefficients) <- colnames(x)

  ## residuals
  opt$fitted.values <- exp(x%*%opt$coefficients)/(1+exp(x%*%opt$coefficients))
  opt$residuals.pearson <- (y-opt$fitted.values)/(sqrt(opt$fitted.values*(1-opt$fitted.values)))
  opt$residuals.pearson <- opt$residuals.pearson[,1]

  ##
  opt$method <- meth
  opt$loglik <- -opt$loglik
  opt$nobs <- n
  opt$df <- m

  #class(opt) <- "logitr"

  return(opt)
}

coef.logitr <- function(object, ...) { ## add probit later
  cf <- object$coefficients
  cf
}

logLik.logitr <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

print.logitr <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("Probit model\n\n")
  if(x$convergence > 0) {
    cat("Model did not converge\n")
  } else {
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
    cat(sprintf("\nLog-likelihood: %s on %s Df\n", format(x$loglik, digits = digits), x$df))
  }
  cat("\n")

  invisible(x)
}

terms.logitr <- function(x, ...) x$terms

model.frame.logitr <- function(formula, ...)
{
  #if(!is.null(formula$model)) return(formula$model)
  formula$terms <- formula$terms
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
}

model.matrix.logitr <- function(object, ...) {
  rval <- model.matrix(object$terms, model.frame(object), contrasts = object$contrasts)
  return(rval)
}

fitted.logitr <- function(object, ...) object$fitted.values

predict.logitr <- function(object, newdata = NULL, na.action = na.pass, ...)
{
  ## obtain model.frame and model.matrix
  if(is.null(newdata)) {
    X <- model.matrix(object)
  } else {
    mf <- model.frame(delete.response(object), newdata, na.action = na.action, xlev = object$levels)
  }

  ## compute result
  rval <- drop(X %*% object$coefficients)
  return(rval)
}

vcov.logitr <- function(object, ...)
{
  vc <- object$vcov
  vc
}

summary.logitr <- function(object, ...)
{

  ## delete stuff
  object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL

  ## return
  class(object) <- "summary.logitr"
  object
}

print.summary.logitr <- function(x, digits = max(3, getOption("digits")-3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(x$convergence > 0L){
    cat("model did not converge\n")
  } else {
    cat(paste("Residuals:\n", sep = ""))
    print(structure(round(as.vector(quantile(x$residuals.pearson)), digits = digits),
                    .Names = c("Min", "1Q", "Median", "3Q", "Max")))

    cat(paste("\nCoefficients:\n", sep = ""))
    print(x$coefficients, digits = digits, signif.legend = FALSE)
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
        "on", sum(sapply(x$coefficients, NROW)), "Df\n")
    cat(paste("Number of iterations in", x$method, "optimization:", x$count[2L], "\n"))

  }
  invisible(x)
}

residuals.logitr <- function(object, ...) {
  cat(paste("Pearson Residuals:\n\n", sep = ""))
  print(object$residuals.pearson)
}

update.logitr <- function (object, formula., ..., evaluate = TRUE)
{
  call <- object$call
  if(is.null(call)) stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if(!missing(formula.)) call$formula <- formula.
  if(length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if(evaluate) eval(call, parent.frame())
  else call
}
