\name{circmax}
\alias{circmax}
\title{Estimate the parameters of the von Mises distribution by maximum likelihood estimation.}
\usage{
circmax(y, x, z = x, start = NULL, ...)
}
\arguments{
\item{y}{Response vector.}

\item{x}{Design matrix for location.}

\item{z}{Design matrix for scale (density parameter kappa).}

\item{start}{Starting values.}

\item{...}{Additional parameters for optim.}
}
\description{
Estimate the parameters of the von Mises distribution by maximum likelihood estimation.
}
