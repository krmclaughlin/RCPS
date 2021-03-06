% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postplot.R
\name{postplot}
\alias{postplot}
\title{Methods to produce plots and diagnostics for posterior distribution fits.}
\usage{
postplot(
  out,
  rdsout,
  burn,
  int,
  covnames,
  N,
  type = c("dens"),
  parallel = FALSE,
  nobself = FALSE
)
}
\arguments{
\item{out}{The outcome of a \code{urInfR} function call.}

\item{rdsout}{The outcome of a \code{simrec} function call, including the RDS data frame, original alpha and beta, and sample size(s).}

\item{burn}{Desired burn-in period.}

\item{int}{Desired thinning interval.}

\item{covnames}{Vector of covariate names (not including self).}

\item{N}{Population size.}

\item{type}{The desired type of outcome plot. Defaults to density. Valid types are \code{dens}, \code{mcmc}, and \code{cor}.}

\item{parallel}{Logical: were multiple (parallel) simulations and inferences performed?}

\item{nobself}{Logical: do not allow peers to reject matches. Defaults to false (peers can reject matches). Not currently used.}
}
\value{
The plot specified by \code{type}.
}
\description{
Allows for the examination of MCMC chains, posterior predictive distributions, and correlation plots between alpha and beta.
}
