circmax <- function(y, x, z = x, start = NULL, ...)
{
   ## Dimensions
   n <- length(y)
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
   
   ## Starting values
   if(is.null(start)) {
      start <- rep(0, m + p)
   } else {
      stopifnot(length(start) == m + p)
   }
   
   ## Optimization
   opt <- stats::optim(par = start, fn = nll, ...)
   
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
   class(opt) <- "circmax"
   
   return(opt)
   }
   
   logLik.circmax <- function(object, ...) {
      structure(object$loglik, df = object$df, class = "logLik")
}

coef.circmax <- function(object, model = c("full", "location", "scale"), ...) {
   model <- match.arg(model)
   cf <- object$coefficients
   switch(model,
      "location" = cf$location,
      "scale" = cf$scale,
      "full" = c(cf$location, cf$scale),
   )
}

print.circmax <- function(x, digits = max(3, getOption("digits") - 3), ...) {
   cat("Maximum likelihood estimation of the von Mises distribution\n\n")
   cat("Coefficients (location model with 2-tan link):\n")
   print.default(format(x$coefficients$location, digits = digits), print.gap = 2, quote = FALSE)
   cat("\nCoefficients (scale model with log link):\n")
   print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
   cat(sprintf("\nLog-likelihood: %s on %s Df\n", format(x$loglik, digits = digits), x$df))
 
   invisible(x)
}
