#
#
# ###
# # simulate data
# ###

# require(MASS)
# set.seed(12072018)
# n <- 5000
# gamma1 <- c(-0.2, 0.2, 0.1)
# gamma2 <- c(0.2, -0.1, 0.2)
# gamma3 <- c(2, 0.5, 0.3)
# beta2 <- 1
# beta3 <- -2
# rho120 <- -0.5
# rho121 <- 0.5
#
# muZ1 <- c(1,-1)
# muZ2 <- c(1,1)
# muZ3 <- c(2,0)
#
# s_Z1 <- matrix(c(3, 0.1, 0.1, 3), nrow = 2, ncol = 2)
# s_Z2 <- matrix(c(3, 0.1, 0.1, 3), nrow = 2, ncol = 2)
# s_Z3 <- matrix(c(4, 0.3, 0.3, 4), nrow = 2, ncol = 2)
#
# Z1 <- matrix(c(rep(1, n),
#                mvrnorm(n = n,
#                        mu = muZ1,
#                        Sigma = s_Z1)),
#              nrow = n,
#              ncol = 3)
# Z2 <- matrix(c(rep(1, n),
#                mvrnorm(n = n,
#                        mu = muZ2,
#                        Sigma = s_Z2)),
#              nrow = n,
#              ncol = 3)
# Z3 <- matrix(c(rep(1, n),
#              mvrnorm(n = n,
#                      mu = muZ3,
#                      Sigma = s_Z3)),
#              nrow = n,
#              ncol = 3)
#
# eps <- rnorm(n, 0, 1)
#
#
# y1star <- Z1 %*% gamma1 + eps
# y1 <- as.numeric(y1star>0)
#
#
# eps0 <- rho120*eps + sqrt(1-rho120^2)*rnorm(n, 0, 1)
# eps1 <- rho121*eps + sqrt(1-rho121^2)*rnorm(n, 0, 1)
#
# y2star <- Z2 %*% gamma2 + y1*beta2 + y1*eps1 + (1-y1)*eps0
# y2 <- as.numeric(y2star>0)
#
# rho230 <- 0.5
# rho231 <- 0.3
#
# y3star <- Z3 %*% gamma3 + beta3*y1 + rnorm(n, 0, 1.5) +
#   y1*eps1*rho231 + (1-y1)*eps0*rho230
# cens <- ifelse(y2 == 1, 1, NA)
# y3 <- y3star*cens
#
# df <- matrix(cbind(y1, y2, y3, Z1[,c(2:3)],Z2[,c(2:3)],Z3[,c(2:3)]), nrow = n,
#              dimnames = list(c(1:n), c("y1","y2", "y3","d1", "d2", "s1", "s2", "x1", "x2")))
#
# df <- as.data.frame(do.call("rbind", replicate(1, df, simplify = FALSE)))
# rownames(df) <- 1:(n)
#
# #undebug(ssdeR)
# require(ssdeR)
#
# #undebug(ssdeR:::ssdeR.cluster.fit)
# m1 <- ssdeR(formula = y3 ~ x1 + x2,
#             treatment = y1 ~ d1 +d2,
#             selection = y2 ~ s1 +s2,
#             data = df)
#
# summary(m1)
#
#
# save(df, file="/Users/michaelbrottrager/Filr/Meine Dateien/PhD_UIBK/Programming in R/ssdeR/data/ssdeRdata.RData")
