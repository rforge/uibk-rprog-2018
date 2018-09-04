numDeriv.gmm <- function(theta,theta_2,H,
                         V1,vcv,X,
                         Z1,Z2,Z3,
                         y1, delta, n,n1, k, k_s, k_d)
{

  dvC11_1 = list()
  dvC12_1 = list()
  dvCo1_1 = list()
  dvCo2_1 = list()
  dvC11 = c()
  dvC12 = c()
  dvCo1 = c()
  dvCo2 = c()


  for(i in 1:(length(theta))) {

    Gamma1 = theta[1:k_d]
    Gamma2 = theta[(k_d+1) : (k_d+k_s)]
    Beta2  = theta[k_d+k_s +1]
    rho120 = theta[k_d+k_s +2]
    rho121 = theta[k_d+k_s +3]

    ele = rep(0, (length(theta)))
    ele[i] = 1

    dtheta = theta+theta[i]*ele*delta
    dgamma1 = dtheta[1:length(Gamma1)]
    dgamma2 = dtheta[(length(Gamma1)+1):(length(Gamma1)+length(Gamma2))]
    dbeta2 = dtheta[(length(Gamma1)+length(Gamma2)) + 1]
    drho120 = dtheta[(length(Gamma1)+length(Gamma2)) + 2]
    drho121 = dtheta[(length(Gamma1)+length(Gamma2)) + 3]

    ZG1 = Z1 %*% dgamma1
    ZG2 = Z2 %*% dgamma2

    C11d = pbinorm(ZG1,ZG2+dbeta2,mean1 = 0,mean2 = 0,cov12 = drho121)
    C11 =  dnorm(ZG1)*pnorm((ZG2+dbeta2-drho121*ZG1)/sqrt(1-drho121^2))/C11d

    C12d = pbinorm(ZG1,ZG2+dbeta2,mean1 = 0,mean2 = 0,cov12 = drho121)
    C12 =  dnorm(ZG2+dbeta2)*pnorm(ZG1-drho121*(ZG2+dbeta2)/sqrt(1-drho121^2))/C12d

    Co1d = pbinorm(-ZG1,ZG2,mean1 = 0,mean2 = 0,cov12 = -drho120)
    Co1 =  -dnorm(ZG1)*pnorm((ZG2-drho120*ZG1)/sqrt(1-drho120^2))/Co1d

    Co2d = pbinorm(-ZG1,ZG2,mean1 = 0,mean2 = 0,cov12 = -drho120)
    Co2 =  dnorm(ZG2)*pnorm((-ZG1+drho120*ZG2)/sqrt(1-drho120^2))/Co2d

    C11p=C11;
    C12p=C12;
    Co1p=Co1;
    Co2p=Co2;

    dtheta = theta-theta[i]*ele*delta
    dgamma1 = dtheta[1:length(Gamma1)]
    dgamma2 = dtheta[(length(Gamma1)+1):(length(Gamma1)+length(Gamma2))]
    dbeta2 = dtheta[(length(Gamma1)+length(Gamma2)) + 1]
    drho120 = dtheta[(length(Gamma1)+length(Gamma2)) + 2]
    drho121 = dtheta[(length(Gamma1)+length(Gamma2)) + 3]

    ZG1 = Z1 %*% dgamma1
    ZG2 = Z2 %*% dgamma2

    C11d = pbinorm(ZG1,ZG2+dbeta2,mean1 = 0,mean2 = 0,cov12 = drho121)
    C11 =  dnorm(ZG1)*pnorm((ZG2+dbeta2-drho121*ZG1)/sqrt(1-drho121^2))/C11d

    C12d = pbinorm(ZG1,ZG2+dbeta2,mean1 = 0,mean2 = 0,cov12 = drho121)
    C12 =  dnorm(ZG2+dbeta2)*pnorm(ZG1-drho121*(ZG2+dbeta2)/sqrt(1-drho121^2))/C12d

    Co1d = pbinorm(-ZG1,ZG2,mean1 = 0,mean2 = 0,cov12 = -drho120)
    Co1 =  -dnorm(ZG1)*pnorm((ZG2-drho120*ZG1)/sqrt(1-drho120^2))/Co1d

    Co2d = pbinorm(-ZG1,ZG2,mean1 = 0,mean2 = 0,cov12 = -drho120)
    Co2 =  dnorm(ZG2)*pnorm((-ZG1+drho120*ZG2)/sqrt(1-drho120^2))/Co2d

    C11m=C11;
    C12m=C12;
    Co1m=Co1;
    Co2m=Co2;

    dvC11_1[[i]] = (C11p-C11m)/(2*delta*theta[i])
    dvC12_1[[i]] = (C12p-C12m)/(2*delta*theta[i])
    dvCo1_1[[i]] = (Co1p-Co1m)/(2*delta*theta[i])
    dvCo2_1[[i]] = (Co2p-Co2m)/(2*delta*theta[i])

  }

  dvC11 <- matrix(unlist(dvC11_1), nrow=nrow(dvC11_1[[1]]), byrow=F)
  dvC12 <- matrix(unlist(dvC12_1), nrow=nrow(dvC12_1[[1]]), byrow=F)
  dvCo1 <- matrix(unlist(dvCo1_1), nrow=nrow(dvCo1_1[[1]]), byrow=F)
  dvCo2 <- matrix(unlist(dvCo2_1), nrow=nrow(dvCo2_1[[1]]), byrow=F)

  # derivatives of Z3'*gamma3+y1*beta3+mu11*y1*C11+mu12*y1*C12+mu01*(1-y1)*C01+mu02*(1-y1)*C02 w.r.t. theta_1

  mu11   = theta_2[ncol(Z3)+2,]
  mu12   = theta_2[ncol(Z3)+3,]
  muo1   = theta_2[ncol(Z3)+4,]
  muo2   = theta_2[ncol(Z3)+5,]

  D_1 =
    ((y1*mu11)%*%t(c(rep(1,ncol(dvC11))))) * dvC11 +
    ((y1*mu12)%*%t(c(rep(1,ncol(dvC12))))) * dvC12 +
    (((1-y1)*muo1)%*%t(c(rep(1,ncol(dvCo1))))) * dvCo1 +
    (((1-y1)*muo2)%*%t(c(rep(1,ncol(dvCo2))))) * dvCo2


  # get the variance matrix/robust
  V_12 = (t(X)%*%D_1/n)*sqrt(n/n1) # the last term accounts for the difference in FS Outcome n
  V1 =  matrix(unlist(V1), ncol = length(theta), byrow = TRUE)
  # RobustVar_2 = vcv + solve(H)%*%V_12%*%V1%*%t(V_12)*t(solve(H))
  RobustVar_2 = vcv + V_12%*%V1%*%t(V_12)

  return(RobustVar_2)
}
