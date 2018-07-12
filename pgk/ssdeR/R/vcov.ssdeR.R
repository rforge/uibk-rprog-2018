
#################
# vcov
#################
vcov.ssdeR <- function (object, ...)
{
  vc <- object$vcov
  k <- length(object$coefficients)
  vc[seq.int(length.out = k), seq.int(length.out = k),
     drop = FALSE]
  colnames(vc) <- rownames(vc) <- names(object$coefficients)
  vc
}
