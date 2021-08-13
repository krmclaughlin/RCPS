## Generate documentation using roxygen

#' Generate recruitment edgelist for a specific wave.
#'
#' Simulation mechanism: generate recruitment edgelist for a specific wave.  Internal function called by sim_mecC.R
#'
#' @param w Wave number.
#' @param N Population size.
#' @param alpha Preference coefficient for recruiters to peers.
#' @param beta Preference coefficient for peers to recruiters.
#' @param ksi Preference coefficient for abstention.
#' @param zeta Preference coefficient for declension.
#' @param elbuild Existing edgelist (for waves 0 to \code{w}-1).
#' @param pindex Index of peers for wave \code{w}.
#' @param ch Covariate(s).  Must have length \code{N}.
#' @param usenetsize Logical: use network size?
#' @param netsiz Vector of network sizes. Must be provided if \code{usenetsize} is true.
#' @param netr Must be provided if \code{usenetsize} is true.
#' @param netp Must be provided if \code{usenetsize} is true.
#'
#' @return Edgelist for wave \code{w}.

RaddEL <- function(w,
                   N,
                   alpha,
                   beta,
                   ksi,
                   zeta,
                   elbuild,
                   pindex,
                   ch,
                   usenetsize,
                   netsiz,
                   netr,
                   netp) {

	nrindex <- elbuild$peer[elbuild$wave == w-1]
	orindexns <- elbuild$peer[elbuild$wave < w]
	pindex1 <- pindex[-which(pindex %in% orindexns)]

	nrec  <- length(nrindex)
	npeer <- length(pindex1)

	crec  <- ch[nrindex, , drop=FALSE]
	cpeer <- ch[pindex1, , drop=FALSE]

	if (usenetsize == TRUE) {
	  recnet <- netsiz[nrindex]
	  recnet <- recnet / max(recnet)
	  peernet <- netsiz[pindex1]
	  peernet <- peernet / max(peernet)
	}

	## Generate utilities
	epsilon <- matrix(rnorm(nrec * (npeer+3), mean=0, sd=1),
					nrow=nrec,
					ncol=npeer+3)
	gamma   <- matrix(rnorm((nrec+1) * npeer, mean=0, sd=1),
					nrow=npeer,
					ncol=nrec+1)

	U <- makeU(crec, cpeer, alpha, ksi, epsilon)		#cpp function
	rownames(U) <- nrindex
	colnames(U) <- c(pindex1, "s1", "s2", "s3")

	V <- makeV(crec, cpeer, beta, zeta, gamma)		#cpp function
	rownames(V) <- pindex1
	colnames(V) <- c(nrindex, "0")

	if (usenetsize == TRUE) {
	  for (i in 1:nrec) {
	    U[i,1:npeer] <- U[i,1:npeer] + netr*peernet
	  }
	  for (j in 1:npeer) {
	    V[j,1:nrec] <- V[j,1:nrec] + netp*recnet
	  }
	}


	## Recruitment

	utiladjmat <- recfromUV3urR(U, V, siml=TRUE)			#r function that calls cpp function
	rownames(utiladjmat) <- nrindex
	colnames(utiladjmat) <- pindex1

	eltemp <- makeEL(utiladjmat, nrindex, pindex1, w)

	return(eltemp)

}
