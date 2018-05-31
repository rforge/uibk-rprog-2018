#vonmises_bamlss <- function(...) {
#   f <- list(
#         "family" = "vonmises",
#         "names" = c("mu", "kappa"),
#         "links" = c("tanhalf", "log"),
#         "d" = function(y, par, log = FALSE) {
#            if (any(par$kappa < 0)) {
#               stop("kappa must be non-negative")
#            }
#            be <- besselI(par$kappa, nu = 0, expon.scaled = TRUE)
#            out <- - log(2 * pi * be) + par$kappa * (cos(y - par$mu) - 1)
#            if (!log) {
#               out <- exp(out)
#            }
#            out
#         },
#         "score" = list(
#            "mu" = function(y, par, ...) {
#               drop(2 * par$kappa * sin(y - par$mu) / ((tan(par$mu/2))^2 + 1) )
#            },
#            "kappa" = function(y, par, ...) {
#               drop(par$kappa * (cos(y - par$mu)
#                    - besselI(par$kappa, nu = 1, expon.scaled = TRUE) / besselI(par$kappa, nu = 0, expon.scaled = TRUE)))
#            }
#         ),
#         "hess" = list(
#            "mu" = function(y, par, ...) {
#               ta <- tan(par$mu/2)
#               drop(4 * par$kappa / (ta^2 + 1)^2 * (sin(y - par$mu) * ta + cos(y - par$mu)))
#            },
#            "sigma" = function(y, par, ...) {
#               be0 <- besselI(par$kappa, nu = 0, expon.scaled = TRUE)
#               be1 <- besselI(par$kappa, nu = 1, expon.scaled = TRUE)
#               be2 <- besselI(par$kappa, nu = 2, expon.scaled = TRUE)
#               drop(-par$kappa * (cos(y - par$mu) + (-2 * be0 * be1  + be0^2 + be2 * be0 - 2 * be1^2)/(2 * be0^2)))
#            }
#         )
#   )
#   class(f) <- "family.bamlss"
#   return(f)
#}

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
