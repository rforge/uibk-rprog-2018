circfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...,
                    estfun = TRUE, object = FALSE, circfit_control = circfit_control()) {

    # TODO: Why is circfit quite often called with object TRUE?!
    # if(object) cat("now object = TRUE\n")

    if(!(is.null(x) || NCOL(x) == 0L)) warning("no regression coefficients 
      are currently taken into account..")
    if(!is.null(offset)) warning("offset not used")

    ny <- NROW(y)
    allequy <- (length(unique(y)) == 1)
    if(is.null(weights) || (length(weights)==0L)) weights <- as.vector(rep.int(1, ny))

    ## control parameters
    control <- circfit_control
    ocontrol <-    control
    solve_kappa <- control$solve_kappa
    control$solve_kappa <- NULL

    # FIXME: change eta to par in all functions and use dvonmises not ddist 

    ## MLE according to Bettina Gruen
    eta <- circmax:::startfun(y, weights = weights, solve_kappa = solve_kappa)
    par <- circmax:::linkinv(eta)

    ## Compute negative loglik
    nll <- -circmax:::ddist(y, eta, log = TRUE, weights = weights, sum = TRUE)

    if(estfun) {
      if(allequy) {
        ef <- matrix(0, ncol = length(eta), nrow = ny)
      } else {
        ef <- as.matrix(weights * circmax:::sdist(y, eta, sum = FALSE))
        n <- NROW(ef)
        ef <- ef/sqrt(n) # TODO: Is this really necessary ?!
      }
    } else {
      ef <- NULL
    }

    # FIXME: Calculate vc or vcov ?!

    list(coefficients = par,
         objfun = nll,
         estfun = ef,
         object = if(object) eta else NULL)
}

ddist <- function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {

  par <- circmax:::linkinv(eta)

  if (any(par[2] < 0)) {
    stop("kappa must be non-negative")
  }
  be <- besselI(par[2], nu = 0, expon.scaled = TRUE)
  val <- - log(2 * pi * be) + par[2] * (cos(y - par[1]) - 1)
  if (!log) {
    val <- exp(val)
  }

  if(sum) {
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
    val <- sum(weights * val, na.rm = TRUE)
  }
  
  return(val)
}

sdist <- function(y, eta, weights = NULL, sum = FALSE) {
  #par <- linkinv(eta)
  par <- c(2 * atan(eta[1]), exp(eta[2]))
  score <- cbind(drop(2 * par[2] * sin(y - par[1]) / ((tan(par[1]/2))^2 + 1) ),
                 drop(par[2] * (cos(y - par[1])
                  - besselI(par[2], nu = 1, expon.scaled = TRUE) / besselI(par[2], nu = 0, expon.scaled = TRUE))))

  score <- as.matrix(score)
  colnames(score) <- c("tan(mu/2)", "log(kappa)")
  if(sum) {
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
    # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN -> gradient is NaN
    score[score==Inf] = 1.7e308
    score <- colSums(weights * score, na.rm = TRUE)
  }
  return(score)
}

startfun <- function(y, weights = NULL, solve_kappa = solve_kappa_Newton_Fourier) {
  x <- cbind(cos(y), sin(y))
  if (is.null(weights) || (length(weights)==0L)) {
      xbar <- colMeans(x)
  } else {
      xbar <- colMeans(weights * x) / mean(weights)
  }
  mu <- atan(xbar[2] /xbar[1]) + (xbar[1] < 0) * sign(xbar[2]) * pi
  rbar <- sqrt(sum(xbar^2))
  kappa <- solve_kappa(rbar)

  starteta <- c(tan(mu / 2), log(kappa))
  names(starteta) <- c("tan(mu/2)", "log(kappa)")
  return(starteta)
}

linkfun <- function(par) {
  eta <- c(tan(par[1] / 2), log(par[2]))
  names(eta) <- c("tan(mu/2)", "log(kappa)")
  return(eta)
}

linkinv <- function(eta) {
  par <- c(2 * atan(eta[1]), exp(eta[2]))
  names(par) <- c("mu", "kappa")
  return(par)
}

## Package movMF implements different kappa solvers:
## o Newton Fourier is used by default.
## o Uniroot seems to provide a safe option.
## o Banerjee_et_al_2005 provides a quick approximation.
solve_kappa_Newton_Fourier <- function (r, tol = 1e-06, maxiter = 100L) {
    lower <- movMF:::Rinv_lower_Amos_bound(r, 0)
    upper <- movMF:::Rinv_upper_Amos_bound(r, 0)
    iter <- 1L
    while (iter <= maxiter) {
        A <- movMF:::A(lower, 2)
        Aprime <- movMF:::Aprime(lower, 2, A = A)
        lower <- lower - (A - r)/Aprime
        A <- movMF:::A(upper, 2)
        upper <- upper - (A - r)/Aprime
        if ((upper - lower) < tol * (lower + upper)) {
            if ((upper - lower) < -tol * (lower + upper)) 
                stop("no convergence")
            break
        }
        iter <- iter + 1L
    }
    (lower + upper)/2
}

solve_kappa_uniroot <- function(r, tol = 1e-06) {
    interval <- c(movMF:::Rinv_lower_Amos_bound(r, 0),
                  movMF:::Rinv_upper_Amos_bound(r, 0))
    if (abs(diff(interval)) < tol) 
        mean(interval)
    else uniroot(function(kappa) movMF:::A(kappa, 2) - r,
                 interval = interval, tol = tol)$root
}

solve_kappa_Banerjee_et_al_2005 <- function(r) {
   r * (2 - r^2)/(1 - r^2)
}


