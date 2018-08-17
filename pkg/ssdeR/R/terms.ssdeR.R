terms.ssdeR <- function (x, model = c("outcome", "treatment", "selection"), ...)
{
  model <- match.arg(model)
  if(model == "outcome"){x$terms[[model]]}else{x$firststage$terms[[model]]}
}

