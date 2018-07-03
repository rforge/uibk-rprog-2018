dgp <- function(n = 1000, beta = c(3, 5, 2), gamma = c(1), control = "version2") {
  m <- length(beta) - 1  # here: number of betas minus intercept
  p <- length(gamma) -1  # here: number of gammas minus intercept
  if(control == "version1") {
    d <- sapply(1:(m + p), function(x) {
      mu <- runif(1, 0, 2 * pi)
      rval <- rnorm(n, mu, 0.2)
      rval}
    )
  } else if(control == "version2") {
    d <- sapply(1:(m + p), function(x) runif(n, -1, 1))
  } else if(control == "version3") {
    d <- sapply(1:(m + p), function(x) rnorm(n, 0, 0.2))
  }
  colnames(d) <- paste0("x", 1:(m + p))

  mu <- 2 * atan(beta[1]) + 2 * atan(crossprod(t(d[, 1:m, drop = FALSE]), beta[-1]))
  if (p > 1) {
    kappa <- exp(gamma[1] + crossprod(t(d[, m + 1:p, drop = FALSE]), gamma[-1]))
  } else {
    kappa <- rep(exp(gamma[1]), n)
  }
  d <- data.frame(d)
  d$y <- NULL
  for(i in 1:n) {
    d[i, "y"] <- circular::rvonmises(1, mu = circular::circular(mu[i]), kappa = kappa[i])
  }
  return(list(data = d, mu = mu, kappa = kappa))
}

plot_circ <- function(x){
  plot(circular::as.circular(x), stack = TRUE, col = colorspace::rainbow_hcl(length(x)))
}

eta.mu0 <- 0
eta.mu <- 2
eta.kappa <- 2

sim_dat <- dgp(n = 1, beta = c(eta.mu0, eta.mu), gamma = eta.kappa)
d <- sim_dat[["data"]]
mu <- sim_dat[["mu"]]
kappa <- sim_dat[["kappa"]]

#par(mfrow = c(3, 1))
#plot(d$y ~ d$x1)
#plot_circ(d$y)

ll_fun <- function(eta.mu0, eta.mu, eta.kappa, y, log = FALSE) {
  mu0 <- 2 * atan(eta.mu0)
  mu <- 2 * atan(eta.mu)
  kappa <- exp(eta.kappa)
  ll <- circmax::dvonmises(y, mu = (mu0 + mu), kappa = kappa, log = TRUE)
  #ll <- NULL
  #for(i in 1:length(y)) {
  #  ll <- c(ll,circular::dvonmises(y[i], as.circular(mu), kappa, log = TRUE))
  #}
  if(log == TRUE) {
    return(sum(ll))
  } else {
    return(exp(sum(ll)))
  }
}


test_eta.mu <- seq(-20, 20, 0.1)
##ll <- sapply(plt_beta, function(x) sum(circmax::dvonmises(d$y, 2 * atan(test_mu) + 2 * atan(crossprod(t(d$x1), x)), kappa, log = TRUE)))
ll <- sapply(test_eta.mu, function(x) ll_fun(eta.mu0, x, eta.kappa, d$y, log = TRUE))
par(mfrow = c(2,1))
plot(ll ~ test_eta.mu, type = "l")
abline(v = 2 * atan(eta.mu0) + 2 * atan(eta.mu), lty = 2)
abline(v = 2 * atan(eta.mu0) + 2 * atan(eta.mu))

abline(v = test_beta - pi, lty = 2)
plot(2*atan(plt_beta/2), ll, type = "l")


ll_fun <- function(eta.mu, eta.kappa, y, log = FALSE) {
  mu <- 2 * atan(eta.mu)
  kappa <- exp(eta.kappa)
  ll <- dvonmises(y, mu = mu, kappa = kappa, log = TRUE)
  #ll <- NULL
  #for(i in 1:length(y)) {
  #  ll <- c(ll,circular::dvonmises(y[i], as.circular(mu), kappa, log = TRUE))
  #}
  if(log == TRUE) {
    return(sum(ll))
  } else {
    return(exp(sum(ll)))
  }
}

eta.mu <- seq(-100, 100, 0.1)
plot(eta.mu, sapply(eta.mu, ll_fun, y = d$y, eta.kappa = 2, log = TRUE), xlab = expression(mu), ylab = "logLik", type = "l")




## Example 2: Periwinkle Dataset of Fisher and Lee, 1992:
distance <- c(107, 46, 33, 67, 122, 69, 43, 30, 12, 25, 37, 69, 5, 83,
  68, 38, 21, 1, 71, 60, 71, 71, 57, 53, 38, 70, 7, 48, 7, 21, 27)
directdeg <- c(67, 66, 74, 61, 58, 60, 100, 89, 171, 166, 98, 60, 197,
  98, 86, 123, 165, 133, 101, 105, 71, 84, 75, 98, 83, 71, 74, 91, 38, 200, 56)
cdirect <- circular(directdeg * 2 * pi/360)
plot(as.numeric(cdirect) ~ distance, ylim = c(0, 4*pi), pch = 20)
points(as.numeric(cdirect) + 2*pi ~ distance, pch = 20)

pm1 <- lm.circular(type = "c-l", y = cdirect, x = distance, init = 0.0)
pm2 <- circmax(cdirect ~ distance, data = data.frame(cbind(distance, cdirect)),
  which.method = "version1", control = circmax_control(trace = 3))
pm3 <- circmax(cdirect ~ distance, data = data.frame(cbind(distance, cdirect)),
  which.method = "default", control = circmax_control(trace = 3))
#pm4 <- circmax(cdirect ~ distance, data = data.frame(cbind(distance, cdirect)), 
#  which.method = "L-BFGS-B", control = circmax_control(trace = 3))

ll_fun <- function(eta.mu, eta.kappa, y, log = FALSE) {
  mu <- 2 * atan(eta.mu)
  kappa <- exp(eta.kappa)
  ll <- dvonmises(y, mu = mu, kappa = kappa, log = TRUE)
  #ll <- NULL
  #for(i in 1:length(y)) {
  #  ll <- c(ll,circular::dvonmises(y[i], as.circular(mu), kappa, log = TRUE))
  #}
  if(log == TRUE) {
    return(sum(ll))
  } else {
    return(exp(sum(ll)))
  }
}

par(mfrow = c(1, 2))
eta.mu <- seq(-100, 100, 0.1)
plot(eta.mu, sapply(eta.mu, ll_fun, y = cdirect, eta.kappa = coef(pm2, model = "scale")[1],
  log = TRUE), xlab = expression(mu), ylab = "logLik", type = "l")

i.max <- which.max(sapply(eta.mu, ll_fun, y = cdirect, eta.kappa = coef(pm2, model = "scale")[1],
  log = TRUE))
eta.kappa <- seq(0, 200, 0.1)
plot(eta.kappa, sapply(eta.kappa, ll_fun, y = cdirect, eta.mu = eta.mu[i.max],
  log = TRUE), xlab = expression(kappa), ylab = "logLik", type = "l")


