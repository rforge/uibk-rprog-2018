# -------------------------------------------------------------------
# - NAME:   crps_vonmises.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-05-21
# -------------------------------------------------------------------
# - PURPOSE: Circular CRPS (von Mises) based on numeric integration using
#            the charististic equation
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-05-22 on thinkmoritz
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

  require(CharFun)

  if(!is.atomic(c(y, mu, kappa)) && is.numeric(c(y, mu, kappa))) 
    stop("Input 'y', 'mu', and 'kappa' must be numeric vectors...")
  if(!(length(y) == length(mu) && length(mu) == length(kappa))) 
    stop("Lengths of 'y', 'mu', and 'kappa' do not match...")
  
  rval <- sapply(seq_along(y), function(i){

    int_fun <- function(x) {
      1 / pi * abs(cfC_vonMises(x, mu[i], kappa[i]) - exp((0+1i) * x * y[i]))^2 / x^2
    }

    crps <- integrate(int_fun, 0, Inf)$value
    return(crps)
  })

  if(sum) rval <- sum(rval)
  return(rval)
}

## Testing
grid <- expand.grid(mu = seq(0, 2 * pi, by = pi / 16), kappa = seq(0.1,10, by = 0.2))

for(i in 1:nrow(grid)){
  grid$crps[i] <- crps_vonmises(pi/2, mu = grid[i, "mu"], kappa = grid[i, "kappa"])
}
#bamlss::plot3d(grid, image = TRUE)
mgrid <- matrix(grid$crps, nrow = 33)
image(unique(grid$mu), unique(grid$kapp), mgrid, col = colorspace::sequential_hcl(31, "Oslo"),
  xlab = "mu [rad]", ylab = "kappa")

#test <- crps_vonmises(rep(pi/2, nrow(grid)), mu = grid[, "mu"], kappa = grid[,"kappa"])


