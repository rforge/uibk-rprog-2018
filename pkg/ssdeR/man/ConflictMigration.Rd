\name{ConflictMigration}
\alias{ConflictMigration}

\title{Excerpt of the Data used in Abel, Crespo-Cuaresma, Brottrager and Muttarak (2018)}

\description{
This is a data frame used by Abel, Crespo-Cuaresma, Brottrager, Muttarak (2018). As the paper is currently under revision, readers are recommended to directly contact <michael.brottrager@jku.at> for a current version of the paper including the detailed data description.

Note that running the model on the data provided below takes some considerable amount of time.
}

\format{
A data frame containing 24336 observations on 21 variables.
\describe{
  \item{iso_i}{ISO code of origin.}
  \item{iso_j}{ISO code of destination.}
  \item{asylum_seekers_ij}{log transformed number of asylum seekers from origin i in destination j.}
  \item{conflict_i}{Conflict in origin i indicated by any reported battle related deaths in that country.}
  \item{isflow_ij}{Non-zero flows between origin i and destination j.}
  \item{stock_ij}{log transformed stock of origin natives in destination j before observational period. (t-1)}
  \item{dist_ij}{Metric distance. (t-1)}
  \item{comlang_ij}{Common Language in both origin and destination (Indicator). (t-1)}
  \item{colony_ij}{Colonial relationship (Indicator). (t-1)}
  \item{polity_i}{normalized (0-1) PolityIV score. (t-1)}
  \item{polity_j}{normalized (0-1) PolityIV score. (t-1)}
  \item{pop_i}{log transformed origin population. (t-1)}
  \item{pop_j}{log transformed destination population. (t-1)}
  \item{gdp_j}{log transformed GDP in destination. (t-1)}
  \item{diaspora_i}{Origin diaspora outside. (t-1)}
  \item{ethMRQ_i}{Ethnic Fractionalization measurement. (t-1)}
  \item{outmigration_i}{log transformed total outmigration of of origin i. (t-1)}
  \item{inmigration_j}{log transformed total inmigration in to destination j. (t-1)}
  \item{spei_i}{12 month average SPEI index. (t-1)}
  \item{battledeaths_i}{log transformed battledeaths in i. (t-1)}
}
}


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
# ----------------------------------------------------------------------------------- #

 data(ConflictMigration, package="ssdeR")

 Results <- ssdeR(formula   =  asylum_seekers_ij ~ stock_ij + dist_ij + I(dist_ij^2) +
                               comlang_ij + colony_ij  +
                               polity_i + pop_i + polity_i + pop_i +
                               gdp_j,
                  treatment =  conflict_i ~ battledeaths_i + spei_i +
                               polity_i + I(polity_i^2) +
                               diaspora_i + ethMRQ_i ,
                  selection =  isflow_ij ~ dist_ij + I(dist_ij^2) +
                               outmigration_i + inmigration_j ,
                  cluster   =  c("iso_i","iso_j"),
                  data      =  ConflictMigration)


 summary(Results)
 marginal.effects.ssdeR(Results, "treatment")
 marginal.effects.ssdeR(Results, "selection")
# ----------------------------------------------------------------------------------- #
}

\keyword{datasets}



