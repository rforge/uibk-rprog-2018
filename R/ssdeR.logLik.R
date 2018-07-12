#######################
# ssdeR.logLik
########################
logLik.ssdeR <- function (object, ...)
{
  structure(object$loglik, df = sum(sapply(object$coefficients,
                                           length)), class = "logLik")
}
