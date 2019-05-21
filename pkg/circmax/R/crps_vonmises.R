# -------------------------------------------------------------------
# - NAME:   crps_vonmises.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-05-21
# -------------------------------------------------------------------
# - PURPOSE: Circular CRPS (von Mises) based on numeric integration using
#            the charististic equation
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-05-21 on thinkmoritz
# -------------------------------------------------------------------

## Function
crps_vonmises <- function(y, mu, kappa) {
  
  require(CharFun)

  int_fun <- function(x) {
    1 / pi * abs(cfC_vonMises(x, mu, kappa) - exp((0+1i) * x * y))^2 / x^2
  }

  rval <- integrate(int_fun, 0, Inf)$value
  return(rval)
}



## Testing
grid <- expand.grid(mu = seq(0, 2 * pi, by = pi / 4), kappa = seq(0.1,10, by = 0.5))

for(i in 1:nrow(grid)){
  grid$crps[i] <- crps_vonmises(pi/2, mu = grid[i, "mu"], kappa = grid[i, "kappa"])
}
bamlss::plot3d(grid, image = TRUE)
