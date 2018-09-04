ssdeR.mfx <- function(x, model=c("outcome","treatment", "selection") , by = NULL)
{
  model <- match.arg(model)

  checkBinaryTrait = function(v, naVal="NA") {
    if (!is.numeric(v)) stop("Only numeric vectors are accepted.")
    vSet = unique(v)
    if (!missing(naVal)) vSet[vSet == naVal] = NA
    vSet = vSet[!is.na(vSet)]

    if (any(as.integer(vSet) != vSet)) FALSE
    else if (length(vSet) > 2) FALSE
    else TRUE
  }

  if(model == "outcome"){
    warning("no marginal effects for outcome equation - Direct/Indirect Effect estimates are provided via ssdeR.effects()")
  }else{
    if(model %in% c("selection")){

      k <- length( coef(x$firststage))
      betas <- coef(x$firststage)[colnames(model.matrix(x, model))]
      rho121 <- tail(coef(x$firststage), 1)
      mm <- model.matrix(x, model)
      mmt <- model.matrix(x, "t")

      m <- c()
      for (i in 1:length(betas)) {
        b <- betas[i]

        if(checkBinaryTrait((mm[, names(b)]))){
          mm[, names(b)] <- 0

          pdf0 <-  mean((dnorm(predict(x, "selection", pnorm = FALSE, newdata = as.data.frame(mm))) *
                           dnorm(predict(x, "treatment", pnorm = FALSE, newdata = as.data.frame(mmt)) - rho121*predict(x, "selection", pnorm = FALSE, newdata = as.data.frame(mm))) /sqrt(1-rho121^2)), na.rm = TRUE) +
            mean((dnorm(predict(x, "treatment", pnorm = FALSE, newdata = as.data.frame(mmt))) *
                    dnorm(predict(x, "selection", pnorm = FALSE, newdata = as.data.frame(mm)) - rho121*predict(x, "treatment", pnorm = FALSE, newdata = as.data.frame(mmt))) /sqrt(1-rho121^2)), na.rm = TRUE)

          mm[, names(b)] <- 1
          pdf1 <-  mean((dnorm(predict(x, "selection", pnorm = FALSE, newdata = as.data.frame(mm))) *
                           dnorm(predict(x, "treatment", pnorm = FALSE, newdata = as.data.frame(mmt)) - rho121*predict(x, "selection", pnorm = FALSE, newdata = as.data.frame(mm))) /sqrt(1-rho121^2)), na.rm = TRUE) +
            mean((dnorm(predict(x, "treatment", pnorm = FALSE, newdata = as.data.frame(mmt))) *
                    dnorm(predict(x, "selection", pnorm = FALSE, newdata = as.data.frame(mm)) - rho121*predict(x, "treatment", pnorm = FALSE, newdata = as.data.frame(mmt))) /sqrt(1-rho121^2)), na.rm = TRUE)

          m[names(b)] <- pdf0 - pdf1
        }
        else{
          pdf <-  mean((dnorm(predict(x, "selection", pnorm = FALSE)) *
                          dnorm(predict(x, "treatment", pnorm = FALSE) - rho121*predict(x, "selection", pnorm = FALSE)) /sqrt(1-rho121^2)), na.rm = TRUE) +
            mean((dnorm(predict(x, "treatment", pnorm = FALSE)) *
                    dnorm(predict(x, "selection", pnorm = FALSE) - rho121*predict(x, "treatment", pnorm = FALSE)) /sqrt(1-rho121^2)), na.rm = TRUE)
          m[names(b)] <- pdf*b

        }
      }
    }

    if (model == "treatment") {
      betas <- coef(x$firststage)[colnames(model.matrix(x, model))]
      mm <- model.matrix(x, model)

      m <- c()
      for (i in 1:length(betas)) {
        b <- betas[i]
        if(checkBinaryTrait((mm[, names(b)]))){

          mm[, names(b)] <- 0
          pdf0 <- mean(dnorm(predict(x, "treatment", pnorm = FALSE, newdata=as.data.frame(mm))), na.rm = TRUE)
          mm[, names(b)] <- 1
          pdf1 <- mean(dnorm(predict(x, "treatment", pnorm = FALSE, newdata=as.data.frame(mm))), na.rm = TRUE)

          pdf <- pdf0 - pdf1
        }else{
          pdf <- mean(dnorm(predict(x, "treatment", pnorm = FALSE)), na.rm = TRUE)
        }
        m[names(b)]  <- pdf*b
      }
    }


    ifelse("(Intercept)" %in% names(m),
           return(m[which(names(m)!= "(Intercept)" )]),return(m))
  }
}
