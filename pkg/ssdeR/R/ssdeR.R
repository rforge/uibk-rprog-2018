ssdeR <-function(formula,
                 treatment,
                 selection,
                 data, subset, na.action = FALSE, weights,
                 cluster = NULL,
                 print.level = 0,
                 control = ssdeR.control(...),
                 model = TRUE, x = FALSE, y = FALSE, ...)
{
  # require(Formula)
  # require(maxLik)
  #
  # require(VGAM)    # for pbinorm()
  # require(mvtnorm) # for rmvnorm()
  # require(fBasics) # for makePositiveDefinite()


  cl <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## What is the role of na.action here?  We cannot use na.omit -- we must not omit the observation
  # where outcome is not observed.  na-s cannot be passed either.
  # However, we can (and should?) omit the na-s in explanatory and probit outcomes.  This needs
  # a bit of refinement.
  mf$na.action <- na.pass


  oformula <- as.formula(formula)
  tmt <- print(treatment[[2]])



  formula <- as.formula(paste(paste(formula[2]), paste(formula[3], tmt, sep = " + "), sep=" ~ "))
  selection <- as.formula(paste(paste(selection[2]), paste(selection[3], tmt, sep = " + "), sep=" ~ "))
  formula <- as.Formula(formula, selection, treatment)

  if (length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~1+1)
    simple_formula <- TRUE
  }
  if (length(formula)[2L] < 3L) {
    formula <- as.Formula(formula(formula), ~1)
    simple_formula <- TRUE
  }
  else {
    if (length(formula)[2L] > 3L) {
      formula <- Formula(formula(formula, rhs = 1:3))
      warning("formula must not have more than 3 RHS parts")
    }
    simple_formula <- FALSE
  }

  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  mt  <- terms(formula, data = data, dot = control$dot)
  mtX3 <- terms(formula, data = data, lhs = 1L, rhs = 1L, dot = control$dot)
  mtX2 <- terms(formula, data = data, lhs = 2L, rhs = 2L, dot = control$dot)
  mtX1 <- terms(formula, data = data, lhs = 3L, rhs = 3L, dot = control$dot)

  Y3 <- model.part(formula, data = mf, lhs = 1, drop = TRUE)
  Y2 <- model.part(formula, data = mf, lhs = 2, drop = TRUE)
  Y1 <- model.part(formula, data = mf, lhs = 3, drop = TRUE)


  X3 <- model.matrix(mtX3, mf)
  X2 <- model.matrix(mtX2, mf)
  X1 <- model.matrix(mtX1, mf)

  # get rid of the endogenous vars
  X3 <- X3[, 1:(NCOL(X3)-1)]
  X2 <- X2[, 1:(NCOL(X2)-1)]

  badRow <- is.na(Y1) | is.na(Y2)
  badRow <- badRow | apply(X1, 1, function(v) any(is.na(v)))
  badRow <- badRow | apply(X2, 1, function(v) any(is.na(v)))
  # check for NA-s.  Because we have to find NA-s in several
  # frames, we cannot use the standard na.* functions here.
  # Find bad rows and remove them later.
  ## Remove NA observations
  badRow <- badRow | (is.na(Y3) & (!is.na(Y2) & Y2 == 1))
  badRow <- badRow | (apply(X3, 1, function(v) any(is.na(v))) & (!is.na(Y2) & Y2 == 1))
  # rows in outcome, which contain NA and are observable -> bad too


  if(print.level > 0) {
    cat(sum(badRow), "invalid observations\n")
  }

  X1 <- X1[!badRow,,drop=FALSE]
  Y1 <- Y1[!badRow]
  X2 <- X2[!badRow,,drop=FALSE]
  Y2 <- Y2[!badRow]
  X3 <- X3[!badRow,,drop=FALSE]
  Y3 <- Y3[!badRow]



  .add_predvars_and_dataClasses <- function(terms, model.frame) {
    rval <- terms
    nval <- if (inherits(model.frame, "terms"))
      model.frame
    else terms(model.frame, dot = control$dot)
    ovar <- sapply(as.list(attr(rval, "variables")), deparse)[-1]
    nvar <- sapply(as.list(attr(nval, "variables")), deparse)[-1]
    if (!all(ovar %in% nvar))
      stop(paste("The following terms variables are not part of the model.frame:",
                 paste(ovar[!(ovar %in% nvar)], collapse = ", ")))
    ix <- match(ovar, nvar)
    if (!is.null(attr(rval, "predvars")))
      warning("terms already had 'predvars' attribute, now replaced")
    attr(rval, "predvars") <- attr(nval, "predvars")[1L +
                                                       c(0L, ix)]
    if (!is.null(attr(rval, "dataClasses")))
      warning("terms already had 'dataClasses' attribute, now replaced")
    attr(rval, "dataClasses") <- attr(nval, "dataClasses")[ix]
    return(rval)
  }

  .add_depvar_names <- function(depvar,formula)
  {
    getResponse <- function(formula) {
      tt <- terms(formula)
      vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
      response <- attr(tt, "response") # index of response var
      vars[response]
    }

    rval <- getResponse(formula)
    depvar <- structure(depvar, varname = rval)
  }

  Y1 <- .add_depvar_names(Y1, treatment)
  Y2 <- .add_depvar_names(Y2, selection)
  Y3 <- .add_depvar_names(Y3, oformula)

  mt   <- .add_predvars_and_dataClasses(mt, mf)
  mtX1 <- .add_predvars_and_dataClasses(mtX1, mf)
  mtX2 <- .add_predvars_and_dataClasses(mtX2, mf)
  mtX3 <- .add_predvars_and_dataClasses(mtX3, mf)


  if (length(Y3) < 1)
    stop("empty model")

  n <- length(Y3)

  weights <- model.weights(mf)

  if( !is.null( weights ) ) {
    if( length( weights ) != length( badRow ) ) {
      stop( "number of weights (", length( weights ), ") is not equal",
            " to the number of observations (", length( badRow ), ")" )
    }
    badRow <- badRow | is.na( weights )
  }else{weights <- 1}
  if (length(weights) == 1)
    weights <- rep.int(weights, length( badRow ))

  bifit <- control$biprob.fit
  gmmfit <- control$gmm.fit
  control$biprob.fit <- control$gmm.fit  <- NULL



  biprob <- do.call(bifit, list(x1 = X1, x2 = X2, x3 = X3,
                                y1 = Y1, y2 = Y2, y3 = Y3,
                                weights = weights,
                                control = control))

  names(biprob$estimate) <- names(biprob$gradient) <- names(biprob$hessian) <- c(colnames(X1), colnames(X2), attr(Y1, "varname"), "rho120", "rho121")


  if(!is.null(cluster)){
    #############################
    nms <- names(Y1)
    V1 <- ssdeR.cluster.fit(fm = biprob, dfcw = 1, cluster = cluster, badRow = badRow, data = data, nms=nms)
  }else{
    V1 <- vcov(biprob)
  }

  # V1 <- vcov(biprob)
  theta_1 = coef(biprob)

  gmm <- do.call(gmmfit, list(x1 = X1, x2 = X2, x3 = X3,
                              y1 = Y1, y2 = Y2, y3 = Y3,
                              V1 = V1,
                              theta = theta_1,
                              weights = weights,
                              control = control))
  #### retrieve the results
  rval <- gmm

  # rval$AIC <- nrow(Y3)*(log(2*pi)+1+log((sum(gmm$residuals^2)/nrow(Y3))))+
  #   ((length(lm_mtcars$coefficients)+1)*2)
  # rval$R2 <- rSquared(Y3, sqrt(diag(gmm$vcov)))
  # rval$adj.R2 <- 1 - (1 - R2) * ((gmm$param$NS1 - 1)/(gmm$param$NS1-length(coef(gmm))))
  #

  # has to be done yet...

  if(!is.null(cluster)){
    nms <- names(Y3[Y2 != 0])
    V2 <- ssdeR.cluster.fit(fm = gmm, dfcw = 1, cluster = cluster, badRow = badRow, data = data, nms=nms)
  }else{
    V2 <- vcov(gmm)
  }

  rval$gmm$vcov <- V2



  rval$firststage <- biprob
  rval$firststage$vcov <- V1

  rval$call <- if (length(control$call))
    control$call
  else cl
  rval$formula <- oformula
  rval$terms <-  list(outcome = mtX3)
  rval$levels <- list(outcome = .getXlevels(mtX3, mf))
  rval$contrasts <- list(outcome = attr(X3, "contrasts"))

  rval$firststage$formula <- formula
  rval$firststage$terms <- list(treatment = mtX1,
                                selection = mtX2)
  rval$firststage$levels <-  list(treatment = .getXlevels(mtX1, mf),
                                  selection = .getXlevels(mtX2, mf))
  rval$firststage$contrasts <- list( treatment =  attr(X1, "contrasts"),
                                     selection =  attr(X2, "contrasts"))

  if (model)
    rval$model <- mf
  rval$firststage$model <- mf
  if (y)
    rval$y <- Y3
  if (x)
    rval$x <- X3
  return(rval)


}
