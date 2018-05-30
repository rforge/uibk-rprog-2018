ic_fit <- function(x, y, control = ic_control(), ...) {

    n <- length(y)

    ### control
    lambda <- control$lambda
    method <- control$method
    diagConst <- control$diagConst

    ### Fit margins
    dens <- density(y)      ## Wann stats::density ???
    F    <- ecdf(y)
    z    <- qnorm(F(y) * n / (n+1))

    MarginConst <- sum(log(dens$y)) - sum(pnorm(z, log = TRUE))

    ### Design matrix
    sc <- mgcv::smooth.construct(mgcv::s(times, bs = "ps"), mcycle, knots = NULL)
    B  <- sc$X
    P  <- sc$S[[1]] + diag(diagConst, ncol(B))

    ### Objective function
    MakeNllh <- function(z, B, P) {
        S    <- diag(1/sqrt(1 + diag(B %*% P %*% t(B))))
        SB   <- S %*% B
        rval <- function(par, lambda) {
            -sum(dnorm(S %*% z, SB %*% par, sqrt(diag(S)), log = TRUE)) +
                 lambda * (t(par) %*% P %*% par)
        }
        rval
    }
    
    ### Optimization
    init <- lm.fit(x = B, y = z)$coef
    nllh <- MakeNllh(z, B, P)
    rval <- optim(init, nllh, method = method, lambda = lambda)

    rval
}

ic_control <- function(lambda = 0, method = "BFGS", diagConst = 0.0001, ...) {

    rval <- list()
    rval$lambda <- lambda
    rval$method <- method
    rval$diagConst <- diagConst

    rval
}

