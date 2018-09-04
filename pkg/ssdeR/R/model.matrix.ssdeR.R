model.matrix.ssdeR <- function (object, model = c("outcome", "treatment", "selection"), ...)
{
  model <- match.arg(model)
  if(model == "outcome"){
    rval <- if (!is.null(object$x[[model]]))
      object$x[[model]]
    else model.matrix(object$terms[[model]], model.frame(object),
                      contrasts = object$contrasts[[model]])
  }else{
    rval <- if (!is.null(object$firststage$x[[model]]))
      object$firststage$x[[model]]
    else model.matrix(object$firststage$terms[[model]], model.frame(object$firststage),
                      contrasts = object$firststage$contrasts[[model]])
  }

  return(rval)
}
