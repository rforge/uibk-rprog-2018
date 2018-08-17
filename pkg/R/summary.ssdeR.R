##################
# summary.ssdeR
##################
summary.ssdeR <- function (object, ...)
{

  k2 <- length(object$coefficients)
  k1 <- length(object$firststage$estimate)

  ko <- NCOL(model.matrix(object, "o"))
  kt <- NCOL(model.matrix(object, "t"))
  ks <- NCOL(model.matrix(object, "s"))

  cf2 <- object$coefficients
  cf1 <- object$firststage$estimate

  mus <- as.vector(do.call("c", as.list(cf2[(ko+1):(k2)])))
  rhos <- as.vector(do.call("c", as.list(cf1[(kt+ks+1):(k1)])))

  cfo <- as.vector(do.call("c", as.list(cf2[1:ko])))            #; names(cfo) <- paste(names(cfo), "o", sep = "_")
  cft <- as.vector(do.call("c", as.list(cf1[1:kt])))            #; names(cfs) <- paste(names(cfs), "t", sep = "_")
  cfs <- as.vector(do.call("c", as.list(cf1[(kt+1):(kt+ks)])))  #; names(cfs) <- paste(names(cfs), "s", sep = "_")


  seo <- sqrt(diag(object$vcov))
  seo_covar <- seo[1:(ko)]                                      #; names(seo_covar) <- paste(names(seo_covar), "o", sep = "_")
  seo_mus <- seo[(ko+1):(k2)]

  set <- sqrt(diag(object$firststage$vcov))
  set_t <- set[1:(kt)]                                          #; names(set_t) <- paste(names(set_t), "o", sep = "_")
  set_s <- set[(kt+1):(kt + ks)]                                #; names(set_t) <- paste(names(set_t), "o", sep = "_")

  set_rho <- set[(kt+ks+1):(k1)]

  cf <- c(cft, cfs, cfo)
  se <- c(set_t, set_s, seo_covar)
  cf <- cbind(cf, se, cf/se,  2 * pnorm(-abs(cf/se)))

  AP <- c(rhos, mus)
  se <- c(set_rho, seo_mus)
  AP <- cbind(AP, se, AP/se,  2 * pnorm(-abs(AP/se)))

  listDF <- list(Parameters=cf, Auxiliaries=AP)
  new_col_name <-  c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  listDF <- lapply(listDF, function(i){
    colnames(i) <- new_col_name
    i
  })

  cf <- listDF[[1]]
  AP <- listDF[[2]]

  if (length(object$coefficients)) {
    cf <- list(treatment = cf[seq.int(length.out = kt), , drop = FALSE],
               selection = cf[seq.int(length.out = ks) + kt, , drop = FALSE],
               outcome   = cf[seq.int(length.out = ko) + kt + ks, , drop = FALSE],
               Aux.Param = AP[, , drop = FALSE])
  }
  else {
   warning("\nNo Parameters found\n")
  }


  object$AIC <- AIC(object)
  object$BIC <- AIC(object, k=sqrt(nObs(object)))

  object$DF1 <- k2
  object$DF2 <- k1

  object$coefficients <- cf
  object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL
  class(object) <- "summary.ssdeR"
  object
}
