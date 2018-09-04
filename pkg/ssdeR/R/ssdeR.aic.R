#######################
# ssdeR.aic
########################
aic.ssdeR <- function (object, ..., k = 2)
{

  return(structure(- 2*logLik(object) + k*object$param$nParam, df = sum(sapply(object$coefficients,
                                        length)), class = "logLik"))
}


