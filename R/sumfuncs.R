## Functions that perform various summary measures

#' Returns a table showing how many recruiters enrolled 0-3 people
#'
#' @param netstage An amended adjacency matrix for a particular wave.  The final row and final column must represent self-matches.
#'
#' @return A table showing how many recruiters enrolled 0, 1, 2, or 3 people into the study.
numRecCount <- function(netstage) {

    nrec  <- nrow(netstage) - 1
    npeer <- ncol(netstage) - 1

    table(rowSums(netstage[1:nrec, 1:npeer], na.rm=TRUE))

}

#' Returns summary plots to assess inference
#'
#' Variety of plots intended to assess performance of inference.  Can combine when rbmpi used for parallelization. Also deals with extra 0s added if Gibbs sampler aborted early.
#'
#' @param out The outcome of a \code{urInfR} call for inference.
#' @param talpha The true value of alpha to compare to.
#' @param tbeta The true value of beta to compare to.
#' @param parallel Logical: were multiple Markov chains used.
#' @param burn Desired burn-in (number of iterations at beginning to disregard).
#' @param interval Desired interval (use every xth iteration).
#' @param gtype Type of graph to produce; defaults to density curves for alpha and beta.  Types are \code{dens}, \code{seqval}, and \code{corplot}.
#'
#' @return The specified type of graph.
infEval <- function(out,
                    talpha,
                    tbeta,
                    parallel=FALSE,
                    burn,
                    interval,
                    gtype="dens") {

  if (parallel == FALSE) {
    nchain <- 1
    nout <- ifelse(is.na(out$earlystop) == FALSE, out$earlystop, length(out$alpha))
  } else {
    nchain <- length(out)
    nout <- ifelse(is.na(out[[1]]$earlystop) == FALSE, out[[1]]$earlystop, length(out[[1]]$alpha))
  }

  sq <- floor((nout - burn) / interval)

  alpha <- c()
  beta <- c()

  if (parallel == FALSE) {

    # Deal with tailing 0s because of early stopping
    if (is.na(out$earlystop) == FALSE) {
      alphastop <- out$alpha[1:out$earlystop, ]
      betastop <- out$beta[1:out$earlystop, ]
    } else {
      alphastop <- out$alpha
      betastop <- out$beta
    }

    alpha.unl <- c(t(alphastop))
    alpha <- alpha.unl[(0:(sq-1))*interval + burn + 1]

    beta.unl <- c(t(betastop))
    beta <- beta.unl[(0:(sq-1))*interval + burn + 1]

  } else {

    for (i in 1:nchain) {

      # Deal with tailing 0s because of early stopping
      if (is.na(out[[i]]$earlystop) == FALSE) {
        alphastop <- out[[i]]$alpha[1:out[[i]]$earlystop, ]
        betastop <- out[[i]]$beta[1:out[[i]]$earlystop, ]
      } else {
        alphastop <- out[[i]]$alpha
        betastop <- out[[i]]$beta
      }

      alpha.unl <- c(t(alphastop))
      tempa <- alpha.unl[(0:(sq-1))*interval + burn + 1]

      beta.unl <- c(t(betastop))
      tempb <- beta.unl[(0:(sq-1))*interval + burn + 1]

      alpha <- c(alpha, tempa)
      beta <- c(beta, tempb)

    }

  }

  # Density plot of alpha and beta
  if (gtype == "dens") {

    par(mfrow = c(1, 2))

    plot(density(alpha), main="Posterior Density of Alpha")
    abline(v=talpha, col="red")

    plot(density(beta), main="Posterior Density of Beta")
    abline(v=tbeta, col="red")

  }

  # Plot of alpha and beta over time
  if (gtype == "seqval") {

    par(mfrow = c(2, 1))

    plot(alpha, main="Posterior Draws of Alpha", type="l")
    abline(h=talpha, col="red")
    if (parallel == TRUE) {
      nlseg <- (0:(nchain))*sq
      abline(v=nlseg, col="gray80")
    }

    plot(beta, main="Posterior Draws of Beta", type="l")
    abline(h=tbeta, col="red")
    if (parallel == TRUE) {
      nlseg <- (0:(nchain))*sq
      abline(v=nlseg, col="gray80")
    }

  }

  # Looking at correlation
  if (gtype == "corplot") {

    par(mfrow = c(1, 1))

    plot(alpha, beta, main="Alpha vs. Beta Values")
    abline(v=talpha, col="red")
    abline(h=tbeta, col="red")

  }

}

#' Tabulate dyad matches.
#'
#' Tabulate if both members of the dyad have the same characteristic value.
#'
#' @param rdf An object of class \code{rds.data.frame}, with one column for the characteristic of interest, named char.
#'
#' @return Table with dyad match possibilities.
dyadmatch <- function(rdf) {

  rec.rdf <- rdf[rdf$wave > 0,]
  num.dyads <- dim(rec.rdf)[1]

  dyad.match <- c()
  for (i in 1:num.dyads) {

    peer.char <- rec.rdf$char[i]
    rec.id <- rec.rdf$recruiter[i]
    rec.char <- rdf$char[which(rdf$peer == rec.id)]

    dyad.match[i] <- peer.char == rec.char

  }

  table(dyad.match)

}

#' Return table of all dyad types
#'
#' @param rdf An object of class \code{rds.data.frame}, with one column for the characteristic of interest, named char.
#'
#' @return A table with all dyad match types
dyadtype <- function(rdf) {

  rec.rdf <- rdf[rdf$wave > 0,]
  num.dyads <- dim(rec.rdf)[1]

  peer.char <- c()
  rec.char <- c()
  for (i in 1:num.dyads) {

    peer.char[i] <- rec.rdf$char[i]
    rec.id <- rec.rdf$recruiter[i]
    rec.char[i] <- rdf$char[which(rdf$peer == rec.id)]

  }

  table(rec.char, peer.char)

}

#' Return information about recruitment homophily, with parallel option.
#'
#' For multiple chains, use \code{sapply(listname, rec.homoph)}.
#'
#' @param simout The output of a simrec function call in \code{prefrecruit}. Must have an object of class \code{rds.data.frame} called \code{rds}, with one column for the characteristic of interest, named char.
#' @param vname Name of outcome variable.
#'
#' @return The recruitment homophily and p-value for a chi-square test of independence.
rec.homoph <- function(simout, vname = "char") {

  rdf <- simout$rds

  stopifnot("rds.data.frame" %in% class(rdf))

  tab <- homophily.estimates(rdf, vname, recruitment = TRUE)@estimate
  pval <- suppressWarnings(chisq.test(tab, correct=FALSE))$p.value

  dg <- diag(tab)
  rs <- rowSums(tab)
  cs <- colSums(tab)
  all <- sum(tab)

  den <- 0
  for (i in 1:nrow(tab)) {
    den = den + rs[i]*cs[i]/all
  }

  homoph <- sum(dg) / den

  ret <- data.frame(homoph, pval)
  names(ret) <- c("homophily", "p-value")

  return(ret)

}
