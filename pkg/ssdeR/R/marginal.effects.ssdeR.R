marginal.effects.ssdeR <- function(object, model = c("outcome", "treatment", "selection"))
{
  model <- match.arg(model)

  if(model %in% c("treatment", "selection") ){
    direct.effect <-  ssdeR.mfx(object, model)
    # delta method std.err
    std.err <- c()
    for (i in names(direct.effect)) {
      v1 <- object$firststage$vcov
      v1 <- v1[paste0(i), paste0(i)]
      v1 <- direct.effect[i] %*% v1 %*% direct.effect[i]
      std.err[paste0(i)] = sqrt(v1)
    }
  }else{
    direct.effect <- coef(object)[2:(length(coef(object))-4)]
    std.err <- c()
    for (i in names(direct.effect)) {
      v1 <- vcov(object)
      v1 <- v1[paste0(i), paste0(i)]
      std.err[paste0(i)] = sqrt(v1)
    }
  }

  effects <- cbind(direct.effect, std.err)
  return(effects)

}


