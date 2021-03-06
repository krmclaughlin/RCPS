% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makepopulation.R
\name{makePopulation}
\alias{makePopulation}
\title{Generate network with full population information}
\usage{
makePopulation(N, density)
}
\arguments{
\item{N}{Population size.}

\item{density}{Target density in the population network.}
}
\value{
A list containing:
 \item{adj}{Population adjacency matrix}
 \item{char}{Covariate values}
}
\description{
Generate a network with full population information.  This is then used to make computation faster for preferential recruitment simulations.  People can only recruit those they are tied to in the underlying network, and thus utilities only need to be updated for these people.
}
\examples{
makePopulation(1000, 0.2)
}
