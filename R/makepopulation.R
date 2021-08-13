#' Generate network with full population information
#'
#' Generate a network with full population information.  This is then used to make computation faster for preferential recruitment simulations.  People can only recruit those they are tied to in the underlying network, and thus utilities only need to be updated for these people.
#'
#' @param N Population size.
#' @param density Target density in the population network.
#'
#' @return A list containing:
#'  \item{adj}{Population adjacency matrix}
#'  \item{char}{Covariate values}
#'
#' @examples
#' makePopulation(1000, 0.2)

makePopulation <- function(N, density) {

  # Generate population covariate
  char <- sample(0:1, N, replace=TRUE)

  # Generate network, making sure there are no isolates
  drawnet <- 1
  while (drawnet == 1) {
    g.net <- network(N, density=density, directed=FALSE)

    if (any(degree(g.net) == 0) == FALSE) {
      drawnet <- 0
    }

  }
  g.net %v% "char" <- char

  adj <- as.matrix.network(g.net)

  return(list(adj=adj, char=char))

}
