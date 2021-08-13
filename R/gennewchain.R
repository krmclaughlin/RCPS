#' Generate a new recruitment chain from posterior medians.
#'
#' Use particular values of \code{alpha}, \code{beta}, \code{xi}, and \code{zeta} (e.g. resulting from fitting an RCPS model) to generate a new recruitment chain with the same seeds as a given RDS sample.
#'
#' @param alphaout Vector of \code{alpha} values.
#' @param betaout Vector of \code{beta} values.
#' @param usenet Logical: use the network size variable? Defaults to false.
#' @param xiout Vector of \code{xi} values.
#' @param zetaout Numerical value of \code{zeta}.
#' @param rdsdf An rds.data.frame object containing the full recruitment chain.
#' @param covname The name of the covariate to be used in the model.
#'
#' @return The outcome of a new \code{simrec} call, which is a list consisting of
#'  \item{net}{network object for the full recruitment chain}
#'  \item{rds}{rds.data.frame object for the full recruitment chain}
#'  \item{rdsfull}{rds.data.frame object including population members not in the sample. For simulation studies. These members of pop not in the sample have recruiter value NA, wave value 1000 (larger than would be encountered).}
#'  \item{seedsrec}{indicator of whether or not the seeds recruited anyone}
#'  \item{alphaorig}{original value of alpha used for simulation}
#'  \item{betaorig}{original value of beta used for simulation}
gennewchain <- function(alphaout,
                         betaout,
                         usenet = FALSE,
                         xiout,
                         zetaout,
                         rdsdf,
                         covname) {
  n <- dim(rdsdf)[1]
  nwav <- max(get.wave(rdsdf))
  stopifnot(class(rdsdf[, covname]) == "factor")
  isseed <- rep(0, n)
  for (i in 1:n) {
    if (get.rid(rdsdf)[i] == get.seed.rid(rdsdf)) {
      isseed[i] <- 1
    }
  }
  nseeds <- length(which(isseed == 1))
  ch.r <- as.matrix(data.frame(char=as.numeric(rdsdf[, covname])))
  degree <- get.net.size(rdsdf)

  if (usenet == TRUE) {
    alpha <- alphaout[1:(length(alphaout) - 1)]
    netr <- alphaout[length(alphaout)]
    beta <- betaout[1:(length(betaout) - 1)]
    netp <- betaout[length(betaout)]
    usenetsize <- TRUE
  } else {
    alpha <- alphaout
    netr <- 0
    beta <- betaout
    netp <- 0
    usenetsize <- FALSE
  }
  ksi <- xiout
  zeta <- zetaout

  tmp=simrec(nseeds, n, ch.r, alpha, beta, ksi, zeta, 10, specifyseeds=TRUE, seedID=which(isseed==1), usenetsize=TRUE, netsiz=degree, netr=netr, netp=netp)

  return(tmp)

}
