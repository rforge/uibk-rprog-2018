ssdeR.control <- function (method = "BHHH", iterlim = NULL,
                           start = NULL, robust = FALSE, ...)
{
  if (is.null(iterlim))
    iterlim <- 5000

  rval <- list(method = method, iterlim = iterlim,
               start = start,
               robust = F,
               biprob.fit = "ssdeRbiprob.fit",
               gmm.fit = "ssdeRgmm.fit")
  rval <- c(rval, list(...))

  if (is.null(rval$reltol))
    rval$reltol <- .Machine$double.eps^(1/1.2)
  rval
}
