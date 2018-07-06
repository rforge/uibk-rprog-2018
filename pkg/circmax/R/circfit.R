circfit <- function(y, x = NULL, start = NULL, subset = NULL, na.action = NULL, 
                    weights = NULL, offset = NULL, ...,
                    estfun = TRUE, object = FALSE, fit_control = circfit_control()) {

    # TODO: Why is circfit quite often called with object TRUE and what are the start values used for?!
    # if(object) cat("now object = TRUE\n")

    # Check unsupported arguments
    if(!(is.null(x) || NCOL(x) == 0L)) warning("no regression coefficients 
      are currently taken into account..")
    if(!is.null(subset)) warning("'subset' currently not supported and therefore not used")
    if(!is.null(na.action)) warning("'na.action' currently not supported and therefore not used")
    if(!is.null(offset)) warning("'offset' currently not supported and therefore not used")

    ## Convenience variables
    n <- NROW(y)
    allequy <- (length(unique(y)) == 1)
    
    ## Weights
    if(is.null(weights) || (length(weights)==0L)) weights <- as.vector(rep.int(1, n))
    if(unique(weights) == 0L) stop("weights are not allowed to be all zero")
    if(length(weights) != n) stop("number of observations and length of weights are not equal")

    ## Control parameters
    solve_kappa <- fit_control$solve_kappa
    fit_control$solve_kappa <- NULL

    ## Get Von Mises family functions (for interchangeability the same distlist is used)
    circfam <- dist_vonmises() 

    ## MLE according to Bettina Gruen
    eta <- circfam$startfun(y, weights = weights, solve_kappa = solve_kappa)
    par <- circfam$linkinv(eta)

    ## Compute negative loglik
    nll <- -circfam$ddist(y, eta, log = TRUE, weights = weights, sum = TRUE)

    if(estfun) {
      if(allequy) {
        ef <- matrix(0, ncol = length(eta), nrow = n)
      } else {
        ef <- as.matrix(weights * circfam$sdist(y, eta, sum = FALSE))
        n <- NROW(ef)
        ef <- ef/sqrt(n) # TODO: Is this really necessary ?!
      }
    } else {
      ef <- NULL
    }

    # FIXME: Calculate vc or vcov ?!
    model <- list()
    
    list(coefficients = par,
         objfun = nll,
         estfun = ef,
         object = if(object) model else NULL)
}
