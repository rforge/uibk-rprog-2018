dvonmises <- function(y, mu, kappa, log = FALSE) {
   if (any(kappa < 0)) {
      stop("kappa must be non-negative")
   }
   be <- besselI(kappa, nu = 0, expon.scaled = TRUE)
   out <- - log(2 * pi * be) + kappa * (cos(y - mu) - 1)
   if (!log) {
      out <- exp(out)
   }
   out
}

ll <- function(par) {
    beta <- par[1:m]
    gamma <- par[m + (1:p)]
    mu <- x[, 1, drop = FALSE] %*% beta[1] + 2 * atan(x[, 
        -1, drop = FALSE] %*% beta[-1])
    kappa <- exp(z %*% gamma)
    ll <- dvonmises(y, mu = mu, kappa = kappa, log = TRUE)
    sum(ll)
}

dvonmises.false <- function(y, mu, kappa, log = FALSE) {
   if (any(kappa < 0)) {
      stop("kappa must be non-negative")
   }
   be <- besselI(kappa, nu = 0, expon.scaled = FALSE)
   out <- - log(be) + kappa * (cos(y - mu))
   if (!log) {
      out <- exp(out)
   }
   out
}

ll.false <- function(par) {
    beta <- par[1:m]
    gamma <- par[m + (1:p)]
    mu <- x[, 1, drop = FALSE] %*% beta[1] + 2 * atan(x[, 
        -1, drop = FALSE] %*% beta[-1])
    kappa <- exp(z %*% gamma)
    ll <- dvonmises.false(y, mu = mu, kappa = kappa, log = TRUE)
    sum(ll)
}

ll.true <- function(par) {
    beta <- par[1:m]
    gamma <- par[m + (1:p)]
    mu <- x[, 1, drop = FALSE] %*% beta[1] + 2 * atan(x[, 
        -1, drop = FALSE] %*% beta[-1])
    kappa <- exp(z %*% gamma)
    ll <- NULL
    for(i in 1:length(y)){
      ll <- c(ll, circular:::DvonmisesRad(y[i], mu = mu[i], kappa = kappa[i], log = TRUE))
    }
    sum(ll)
}

distance <- c(107, 46, 33, 67, 122, 69, 43, 30, 12, 25, 37, 69, 5, 83,
  68, 38, 21, 1, 71, 60, 71, 71, 57, 53, 38, 70, 7, 48, 7, 21, 27)
directdeg <- c(67, 66, 74, 61, 58, 60, 100, 89, 171, 166, 98, 60, 197,
  98, 86, 123, 165, 133, 101, 105, 71, 84, 75, 98, 83, 71, 74, 91, 38, 200, 56)
cdirect <- circular::circular(directdeg * 2 * pi/360)

pm1 <- circular::lm.circular(type = "c-l", y = cdirect, x = distance, init = 0.0)
par1 <- c(as.numeric(pm1$mu), pm1$coefficients, log(pm1$kappa))
debugonce(circmax_fit)
pm3 <- circmax(cdirect ~ distance, data = data.frame(cbind(distance, cdirect)),
 which.method = "default", control = circmax_control(trace = 3))

ll(par1)
ll.true(par1)
ll.false(par1)
pm1$log.lik
