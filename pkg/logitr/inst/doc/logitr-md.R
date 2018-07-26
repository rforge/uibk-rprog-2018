## ----preliminaries, echo = FALSE, message = FALSE, results = "hide"------
library("logitr")
library("numDeriv")

## ----logitr vs glm-------------------------------------------------------
## generating data
set.seed(123)
x1 <- rnorm(30,3,2) + 0.1 * c(1:30)
x2 <- rbinom(30,1,0.3)
x3 <- rpois(n=30,lambda = 4)
x3[16:30] <- x3[16:30] - rpois(n=15, lambda = 2)
xdat <- cbind(x1,x2,x3)
ydat <- c(rbinom(5,1,0.1), rbinom(10,1,0.25), rbinom(10,1,0.75), rbinom(5,1,0.9))
## comparison of model outputs
(m0 <- logitr(ydat~xdat))
(m1 <- glm(ydat~xdat, family = "binomial"))
## comparing AIC and BIC
AIC(m0)
AIC(m1)
BIC(m0)
BIC(m1)
  

