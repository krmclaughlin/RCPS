% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sumfuncs.R
\name{infEval}
\alias{infEval}
\title{Returns summary plots to assess inference}
\usage{
infEval(out, talpha, tbeta, parallel = FALSE, burn, interval, gtype = "dens")
}
\arguments{
\item{out}{The outcome of a \code{urInfR} call for inference.}

\item{talpha}{The true value of alpha to compare to.}

\item{tbeta}{The true value of beta to compare to.}

\item{parallel}{Logical: were multiple Markov chains used.}

\item{burn}{Desired burn-in (number of iterations at beginning to disregard).}

\item{interval}{Desired interval (use every xth iteration).}

\item{gtype}{Type of graph to produce; defaults to density curves for alpha and beta.  Types are \code{dens}, \code{seqval}, and \code{corplot}.}
}
\value{
The specified type of graph.
}
\description{
Variety of plots intended to assess performance of inference.  Can combine when rbmpi used for parallelization. Also deals with extra 0s added if Gibbs sampler aborted early.
}
