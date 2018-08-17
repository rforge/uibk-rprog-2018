predict.ssdeR <- function ( object, model = c("outcome", "treatment", "selection"), pnorm = F ,newdata = NULL, ... ) {
  model <- match.arg(model)

  if( is.null( newdata ) ) {
    mm <- model.matrix( object, model )
  } else {  # modified copy from predict.lm()
    tt <- terms( object , model)
    Terms <- delete.response( tt )
    m <- model.frame( Terms, newdata, xlev = object$xlevels )
    mm <- model.matrix( Terms, m , model)
  }

  if(model == "outcome"){
    result <- drop(object$X %*% coef( object ) )
  }else{
    nm <- colnames(mm)
    result <- drop( mm %*% coef( object$firststage )[nm] )
  }



  if( model  %in% c("treatment","selection", pnrom = T )) {
    result <- pnorm( result )
  }


}
