############################
# Robust VCV
###########################
ssdeR.vcv <- function(X,fm ,theta, theta_2, V1,k_s,k_d,
                      Z1,Z2,Z3,y1,y2,y3, delta, robust)
{
  ## Define n and k parameters
  n = nrow(X)
  n1 = length(y1)
  k = ncol(X)
  H = 1/(n)*t(X)%*%X

  #Calculate vector of residuals
  e3 = y3 - X%*%theta_2
  ## Calculate Variance-Covariance Matrix
  VCV = 1/(n) * as.numeric(t(e3)%*%e3) * solve(t(X)%*%X)

  vcov = numDeriv.gmm(theta = theta,
                      theta_2 = theta_2,
                      H = H,
                      V1 = V1,
                      vcv = VCV,
                      X = X,
                      Z1 = Z1,
                      Z2 = Z2,
                      Z3 = Z3,
                      y1 = y1,
                      delta = delta,
                      n = n,
                      n1 = n1,
                      k = k,
                      k_s = k_s,
                      k_d = k_d)
  return(vcov)
}

