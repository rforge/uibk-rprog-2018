\name{marginal.effects.ssdeR}
\alias{marginal.effects.ssdeR}

\title{
marginal.effects.ssdeR}
\description{
\code{marginal.effects.ssdeR} calculates and returns the marginal effects associated with either outcome, treatment or selection models of \code{ssdeR} objects.
}
\usage{
marginal.effects.ssdeR(object, model = c("outcome", "treatment", "selection"))
}
\arguments{
  \item{object}{
\code{ssdeR} object
}
  \item{model}{
Either outcome, treatment or selection model.}
}
\details{
\code{marginal.effects.ssdeR} computes the marginal effects based on Mullahy, J. (2017).

If \code{model = "selection"} is chosen, \code{marginal.effects.ssdeR} returns the marginal effects in the bivariate probit model.
In case, \code{model = "treatment"}, the marginal effects computation reduces to simple probit marginal effects and in case \code{model = "outcome"}, simple 3rd-step parameter estimates are returned.

As \code{ssdeR} estimates a bivariate probit model with structural shift, \code{selection} model indirect effects are just the \code{treatment} model's direct effects.

Standard errors are computed using the delta method.
}
\value{
Returns a n x 2 matrix containing the n (marginal) direct effects of each model. If no \code{model} is chosen, \code{marginal.effects.ssdeR} just returns the parameter estimates of the outcome equation.

}

\references{
  Cameron, A. C. and Trivedi, P. K. (2005)
  \emph{Microeconometrics: Methods and Applications},
  Cambridge University Press.

  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.

  Heckman, J. (1976)
  The common structure of statistical models of
  truncation, sample selection and limited dependent
  variables and a simple estimator for such models,
  \emph{Annals of Economic and Social Measurement},
  5(4), p. 475-492.

  Johnston, J. and J. DiNardo (1997)
  \emph{Econometric Methods, Fourth Edition}, McGraw-Hill.

  Lee, L., G. Maddala and R. Trost (1980)
  Asymetric covariance matrices of two-stage probit and two-stage tobit
  methods for simultaneous equations models with selectivity.
  \emph{Econometrica}, 48, p. 491-503.

  Mullahy, J. (2017) Marginal effects in multivariate probit models.  \emph{Empircal Economics}, 52: 447.

  il Kim, K. (2006). Sample selection models with a common dummy endogenous regressor in simultaneous equations: A simple two-step estimation. \emph{Economics Letters}, 91(2), 280-286.


  Petersen, S., G. Henningsen and A. Henningsen (2017)
  \emph{Which Households Invest in Energy-Saving Home Improvements?
  Evidence From a Danish Policy Intervention}.
  Unpublished Manuscript.
  Department of Management Engineering, Technical University of Denmark.

  Toomet, O. and A. Henningsen, (2008)
  Sample Selection Models in R: Package sampleSelection.
  \emph{Journal of Statistical Software} 27(7),
  \url{http://www.jstatsoft.org/v27/i07/}

  Wooldridge, J. M. (2003)
  \emph{Introductory Econometrics: A Modern Approach, 2e}, Thomson South-Western.}
\author{
Michael Brottrager}

\examples{
# ----------------------------------------------------------------------------------- #
# 1. Climate, Conflict, Crossborder Migration Data
# ----------------------------------------------------------------------------------- #
# This is the data used Abel, Crespo-Cuaresma, Brottrager, Muttarak (2018).
# As the paper is currently under revision, readers are recommended to directly
# contact <michael.brottrager@jku.at> for a current version of the paper including
# the detailed data description.
#
# Note that running the model on the data provided below takes some
# considerable amount of time.
#
# Please consider looking at the accompaning vignette.
#
# ----------------------------------------------------------------------------------- #
# data(ConflictMigration, package="ssdeR")
#
# Results <- ssdeR(formula = asylum_seekers_ij ~ stock_ij + dist_ij + I(dist_ij^2) +
#                    comlang_ij + colony_ij  +
#                    polity_i + pop_i + polity_i + pop_i +
#                    gdp_j,
#                  treatment = conflict_i ~ battledeaths_i + spei_i +
#                    polity_i + I(polity_i^2) +
#                    diaspora_i + ethMRQ_i ,
#                  selection = isflow_ij ~ dist_ij + I(dist_ij^2) +
#                    outmigration_i + inmigration_j ,
#                  cluster = c("iso_i","iso_j"),
#                  data = ConflictMigration)
#
# summary(Results)
# marginal.effects.ssdeR(Results, "treatment")
# marginal.effects.ssdeR(Results, "selection")
# marginal.effects.ssdeR(Results, "outcome")
# ----------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------- #
# 2.  Simulation Data
# ----------------------------------------------------------------------------------- #
library(MASS)
set.seed(12072018)
n <- 5000
gamma1 <- c(-0.2, 0.2, 0.1)
gamma2 <- c(0.2, -0.1, 0.2)
gamma3 <- c(2, 0.5, 0.3)
beta2 <- 1
beta3 <- -2
rho120 <- -0.5
rho121 <- 0.5

muZ1 <- c(1,-1)
muZ2 <- c(1,1)
muZ3 <- c(2,0)

s_Z1 <- matrix(c(3, 0.1, 0.1, 3), nrow = 2, ncol = 2)
s_Z2 <- matrix(c(3, 0.1, 0.1, 3), nrow = 2, ncol = 2)
s_Z3 <- matrix(c(4, 0.3, 0.3, 4), nrow = 2, ncol = 2)

Z1 <- matrix(c(rep(1, n), mvrnorm(n = n, mu = muZ1, Sigma = s_Z1)), nrow = n, ncol = 3)
Z2 <- matrix(c(rep(1, n), mvrnorm(n = n, mu = muZ2, Sigma = s_Z2)), nrow = n, ncol = 3)
Z3 <- matrix(c(rep(1, n), mvrnorm(n = n, mu = muZ3, Sigma = s_Z3)), nrow = n, ncol = 3)

eps <- rnorm(n, 0, 1)


y1star <- Z1 \%*\% gamma1 + eps
y1 <- as.numeric(y1star>0)


eps0 <- rho120*eps + sqrt(1-rho120^2)*rnorm(n, 0, 1)
eps1 <- rho121*eps + sqrt(1-rho121^2)*rnorm(n, 0, 1)

y2star <- Z2 \%*\% gamma2 + y1*beta2 + y1*eps1 + (1-y1)*eps0
y2 <- as.numeric(y2star>0)

rho230 <- 0.5
rho231 <- 0.3

y3star <- Z3 \%*\% gamma3 + beta3*y1 + rnorm(n, 0, 1.5) +
  y1*eps1*rho231 + (1-y1)*eps0*rho230
cens <- ifelse(y2 == 1, 1, NA)
y3 <- y3star*cens

df <- matrix(cbind(y1, y2, y3, Z1[,c(2:3)],Z2[,c(2:3)],Z3[,c(2:3)]), nrow = n,
             dimnames = list(c(1:n), c("y1","y2", "y3",
                                       "d1", "d2", "s1",
                                       "s2", "x1", "x2")))

df <- as.data.frame(do.call("rbind", replicate(1, df, simplify = FALSE)))
rownames(df) <- 1:(n)

m1 <- ssdeR(formula = y3 ~ x1 + x2,
            treatment = y1 ~ d1 +d2,
            selection = y2 ~ s1 +s2,
            data = df)

summary(m1)
marginal.effects.ssdeR(m1, "treatment")
marginal.effects.ssdeR(m1, "selection")
marginal.effects.ssdeR(m1, "outcome")
# ----------------------------------------------------------------------------------- #

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
