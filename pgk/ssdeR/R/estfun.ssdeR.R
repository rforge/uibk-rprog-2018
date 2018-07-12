#############
# estfun.ssdeR
#############
estfun.ssdeR <- function(object, ...)
{
  xmat <- object$X
  xmat <- naresid(object$na.action, xmat)
  if(any(alias <- is.na(coef(object)))) xmat <- xmat[, !alias, drop = FALSE]
  wts <- weights(object)
  if(is.null(wts)) wts <- 1
  res <- residuals(object)
  rval <- as.vector(res) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  if(is.zoo(res)) rval <- zoo(rval, index(res), attr(res, "frequency"))
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  return(rval)
}
