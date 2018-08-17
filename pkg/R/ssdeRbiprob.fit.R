#########################
# 3. ssdeR.fit
########################
ssdeRbiprob.fit <- function(y1, y2, y3,
                            x1, x2, x3,
                            weights = weights,
                            control = control)
{
  ocontrol <- control
  n = NROW(x3)
  k = NCOL(x3)
  n_s = NROW(x2)
  k_s = NCOL(x2)
  n_d = NROW(x1)
  k_d = NCOL(x1)

  if (is.null(weights))
    weights <- rep.int(1, n)
  nobs <- sum(weights > 0)


  NT1 <- sum(y1 == 1)  # number of treated
  NT0 <- nobs - NT1    # number of untreated
  NS1 <- sum(y2 == 1)  # number of selected
  NS0 <- nobs - NS1    # number of untreated

  if(NT0 == 0 | NT1 == 0 ) {
    stop("No variance in the treatment variable")
  }
  if(NS0 == 0 | NS1 == 0 ) {
    stop("No variance in the selection variable")
  }

  bilvls <- function(x){
    ylvl <- levels(as.factor(x))
    if( length( ylvl ) != 2 ) {      # restrict to binary case!
      stop( "the left hand side of the 'formula' has to contain",
            " exactly two levels (e.g. FALSE and TRUE)" )
    }
    x <- as.integer(x == ylvl[ 2 ])    # ... convert to integers F==0, T==1
  }

  y1 <- bilvls(y1)
  y2 <- bilvls(y2)

  ll <- function(pars){

    Gamma1 = pars[1:k_d]
    Gamma2 = pars[(k_d+1) : (k_d+k_s)]
    beta2 = pars[k_d+k_s +1]
    rho120 = (exp(2*pars[k_d+k_s +2]) -1)/(exp(2*pars[k_d+k_s +2]) +1)
    rho121 = (exp(2*pars[k_d+k_s +3]) -1)/(exp(2*pars[k_d+k_s +3]) +1)

    ZG1 = x1 %*% Gamma1
    ZG2 = cbind(x2, y1) %*% c(Gamma2, beta2)

    #  pbinorm gives the cumulative distribution function
    Loo = pmax(pbinorm(-ZG1,   -ZG2      , mean1 = 0, mean2 = 0 ,  cov12 =  rho120), delta)
    Lol = pmax(pbinorm(-ZG1,    ZG2      , mean1 = 0, mean2 = 0 ,  cov12 = -rho120), delta)
    Llo = pmax(pbinorm( ZG1,   -ZG2, mean1 = 0, mean2 = 0 ,  cov12 = -rho121), delta)
    Lll = pmax(pbinorm( ZG1,    ZG2, mean1 = 0, mean2 = 0 ,  cov12 =  rho121), delta)

    LL = (1-y1)*(1-y2)*log(Loo) +
      (1-y1)*(  y2)*log(Lol) +
      (  y1)*(1-y2)*log(Llo) +
      (  y1)*(  y2)*log(Lll)

    return(LL)
  }

  method <- control$method
  start <- control$start
  robust <- control$robust
  control$method <- control$start <- control$robust <-  NULL
  delta <- control$reltol

  if(is.null(start)){
    start <- c(coef(glm.fit(x1,y1, family = binomial(link = "probit"))),
               coef(glm.fit(x2,y2),family = binomial(link = "probit")),
               rep(0,3))
  } else if(length(start)<(k_s + k_d + 3)){warning(paste("too few entries in start! Please provide",
                                                         (k_s + k_d + 3), "starting values"))}
  if (is.list(start))
    start <- do.call("c", start)
  if (length(start) > (k_s + k_d + 3)) {
    warning(paste("too many entries in start! only first",
                  (k_s + k_d + 3), "entries are considered"))
    start <- start[1:(k_s + k_d + 3)]
  }


  opt <- suppressWarnings(maxLik(logLik=ll, method = method, start=start, control=control))


  return(opt)
}
