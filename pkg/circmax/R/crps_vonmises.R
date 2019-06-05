# -------------------------------------------------------------------
# - NAME:   crps_vonmises.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-05-21
# -------------------------------------------------------------------
# - PURPOSE: Circular CRPS (von Mises) based on numeric integration using
#            the charististic equation
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-06-05 on thinkmoritz
# -------------------------------------------------------------------

### Function
#crps_vonmises <- function(y, mu, kappa) {
#
#  require(CharFun)
#  
#  int_fun <- function(x) { 
#    1 / pi * abs(cfC_vonMises(x, mu, kappa) - exp((0+1i) * x * y))^2 / x^2
#  } 
#    
#  rval <- integrate(int_fun, 0, Inf)$value
#  return(rval)
#}

## Function
crps_vonmises <- function(y, mu, kappa, sum = FALSE) {

  if(any(y < -pi) || any(y > pi) || any(mu < -pi) || any(mu > pi) || any(kappa < 0 ))
    stop("y and mu must be in the interval of [-pi, pi], and kappa must be non negative!") 

  require(CharFun)

  if(!inherits(y, c("numeric", "integer")) || !inherits(mu, c("numeric", "integer")) ||
    !inherits(kappa, c("numeric", "integer"))) {
    stop("Input 'y', 'mu', and 'kappa' must be numeric vectors...")
  }

  dat <- data.frame("y" = y, "mu" = mu, "kappa" = kappa)

  idx <- which.min(abs(dat$mu + c(-2 * pi, 0, 2 * pi) - dat$y))
  dat$mu <- dat$mu + c(-2 * pi, 0, 2 * pi)[idx]
  
  rval <- sapply(1:nrow(dat), function(i){

    int_fun <- function(x) {
      1 / pi * abs(CharFun::cfC_vonMises(x, dat[i, "mu"], dat[i, "kappa"]) - 
        exp((0+1i) * x * dat[i, "y"]))^2 / x^2
    }

    crps <- integrate(int_fun, 0, Inf)$value
    return(crps)
  })

  if(sum) rval <- sum(rval)
  return(rval)
}

### Testing
#grid <- expand.grid(mu = seq(0, 2 * pi, by = pi / 16), kappa = seq(0.1,10, by = 0.2))
#
#for(i in 1:nrow(grid)){
#  grid$crps[i] <- crps_vonmises(pi/2, mu = grid[i, "mu"], kappa = grid[i, "kappa"])
#}
##bamlss::plot3d(grid, image = TRUE)
#mgrid <- matrix(grid$crps, nrow = 33)
#image(unique(grid$mu), unique(grid$kapp), mgrid, col = colorspace::sequential_hcl(31, "Oslo"),
#  xlab = "mu [rad]", ylab = "kappa")
#
#test <- crps_vonmises(rep(pi/2, nrow(grid)), mu = grid[, "mu"], kappa = grid[,"kappa"])


