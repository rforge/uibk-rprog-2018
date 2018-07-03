## Helper functions
plot_circ <- function(x, ...){
  require("circular", quietly = TRUE)
  suppressWarnings(x <- as.circular(x))
  plot(x, cex = 1., bin = 720, stack = TRUE, sep = 0.035, shrink = 1.3, ...)
  rose.diag(x, bins = 16, col = "darkgrey", cex = 1., prop = 1.3, add = TRUE)
  lines(density.circular(x, bw = 1), lty = 2)
}

## Create observations
n <- 100
alpha0 <- 0
beta0 <- 0
beta1 <- 1
#gamma0 <- 0

#x1 <- runif(n, -1, 1)
x1 <- seq(-1, 1, length.out = n)

mu <- 2 * atan(alpha0) + 2 * atan(beta0 + beta1 %*% x1)
#kappa <- rep(exp(gamma), n)

#obs <- NULL
#for( i in 1:n){
#  obs <- c(obs, suppressWarnings(circular::rvonmises(1, mu = mu[i], kappa = kappa[i])))
#}

## Create Pseudo Forecasts
test_alpha0 <- tan(0.75 *pi)

mu.fcst <- 2 * atan(test_alpha0) + 2 * atan(beta0 + beta1 %*% x1)
#kappa.fcst <- rep(exp(gamma), n)

#fcst <- NULL
#for( i in 1:n){
#  fcst <- c(fcst, suppressWarnings(circular::rvonmises(1, mu = mu[i], kappa = kappa[i])))
#}

par(mfrow = c(1,2))
# PLot observations
plot_circ(as.numeric(mu), main = sprintf("VM(y | alpha0 = %d, beta0 = %d, beta1 = %d)", alpha0, beta0, beta1))

# PLot observations
plot_circ(as.numeric(mu.fcst), main = sprintf("VM(y | alpha0 = %d, beta0 = %d, beta1 = %d)", alpha0, beta0, beta1))

col_hcl <- c(gray(0, alpha = 0.2), colorspace::rainbow_hcl(2, alpha =0.5))
plot(seq(-50, 50, 0.1),  as.numeric(2 * atan(alpha0) + 2 * atan(beta0 + beta1 %*% seq(-50, 50, 0.1))), 
  type = "line", xlab = "x1", ylab = "mu", ylim = c(-2*pi, 2*pi))
lines(seq(-50, 50, 0.1),  as.numeric(2 * atan(test_alpha0) + 2 * atan(beta0 + beta1 %*% seq(-50, 50, 0.1)))
  , lty = 2, type = "line", xlab = "x1", ylab = "mu", ylim = c(-2*pi, 2*pi))
points(x1, as.numeric(mu), pch = 20, col = col_hcl[1], cex = 0.5)
points(x1, as.numeric(mu.fcst), pch = 20, col = col_hcl[3], cex = 0.5)
rug(as.numeric(mu), side = 2, col = col_hcl[1])
rug(as.numeric(mu.fcst), side = 2, col = col_hcl[3])
legend("topright", c("observation", "misspecified fcst"), pch = 20, col = col_hcl[c(1,3)])


ll <- sapply(test_eta.mu, function(x) sum(circmax::dvonmises(y, 2 * atan(eta.mu0) + 2 * atan(x %*% x1), kappa, log = TRUE)))
plot(ll ~ test_eta.mu, type = "l")
abline(v = eta.mu, lty = 2)
abline(v = tan((2*atan(eta.mu) - pi)/2), lty = 2)

plot(2*atan(test_eta.mu), ll, type = "l")

# PLot observations
plot_circ(y, main = sprintf("VM(y | beta0 = %d, beta1 = %d, beta2 = %d, gamma = %s)", beta0, beta1, beta2, gamma))
###
eta.mu0 <- 3
eta.mu <- 3
eta.kappa <- 1

mu <- 2 * atan(eta.mu0) + 2 * atan(eta.mu)
kappa <- exp(eta.kappa)

test_eta.mu <- seq(-20, 20, 0.1)
ll <- sapply(test_eta.mu, function(x) sum(circmax::dvonmises(0, 2 * atan(eta.mu0) + 2 * atan(x), kappa, log = TRUE)))
plot(ll ~ test_eta.mu, type = "l")

abline(v = eta.mu, lty = 2)
abline(v = tan((2*atan(eta.mu) - pi)/2), lty = 2)

plot(2*atan(test_eta.mu), ll, type = "l")
