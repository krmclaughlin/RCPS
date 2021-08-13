#' Make annotated adjacency matrix for a given wave.
#'
#' Create a "netstage" (amended adjacency matrix) matrix for a given wave of recruitment.  For \code{n_r} recruiters and \code{n_p} peers in wave \code{w}, this is a \eqn{(n_r+1) \times (n_p+3)} matrix, where the final row and final three columns contain self-match information.  Should be called in a urInfR run.
#'
#' @param rdsdf An rds.data.frame object containing the full recruitment chain.  Must also contain information about covariate(s) for each person.
#' @param covnames Vector of names of covariates in \code{rdf} to be used in the model.
#' @param wave The wave of the recruiter, so choose 0 for wave where seeds recruit wave 1.
#' @param restrict.local Logical: use additional information about the underlying network? Makes computation faster.  Must also specify \code{popadjmat}.
#' @param popadjmat Population adjacency matrix.  Must be supplied if \code{restrict.local=TRUE}.
#' @param usenetsize Logical: use network size? Defaults to false.
#'
#' @return A list containing
#'  \item{netstage}{the annotated adjacency matrix, containing observed recruitment information and self-match information for the given wave.  Assumed to take two values (0-1) if \code{restrict.local=FALSE}, or three values (0-1-NA) if \code{restrict.local=TRUE}.}
#'  \item{chrec}{the covariate values for this wave's recruiters}
#'  \item{chpeer}{the covariate values for this wave's peers}

makeNetStageR <- function(rdsdf,
                          covnames,
                          wave,
                          restrict.local=FALSE,
                          popadjmat=NULL,
                          usenetsize=FALSE) {

	dfpart <- rdsdf[rdsdf$wave >= wave, ]

	rec.id  <- dfpart$peer[dfpart$wave == wave]
	peer.id <- dfpart$peer[dfpart$wave > wave]

	nrec  <- length(rec.id)
	npeer <- length(peer.id)

	netstage <- matrix(0, nrow=nrec+1, ncol=npeer+3)
	rownames(netstage) <- c(rec.id, "0")
	colnames(netstage) <- c(peer.id, "s1", "s2", "s3")

	recID <- dfpart$recruiter[dfpart$wave > wave]
	recIDct <- rep(0, npeer)

	for (i in 1:npeer) {
	  idtemp <- which(rec.id == recID[i])
		recIDct[i] <- ifelse(length(idtemp) == 0, NA, idtemp)
	}

	peerIDct <- 1:npeer

	rpdf <- na.omit(data.frame(recIDct, peerIDct))
	recIDf  <- rpdf$recIDct
	peerIDf <- rpdf$peerIDct

	netstageF <- makeNetstageC(netstage, recIDf, peerIDf)		#cpp function

	# Return covariate information
	nchar <- length(covnames)
	chrec <- matrix(nrow=nrec, ncol=nchar)
	chpeer <- matrix(nrow=npeer, ncol=nchar)
	for (i in 1:nchar) {
	  colid <- which(names(rdsdf) == covnames[i])
	  chrec[, i] <- dfpart[dfpart$wave == wave, colid]
	  chpeer[, i] <- dfpart[dfpart$wave > wave, colid]
	}
	neid <- which(names(dfpart) == attr(dfpart, "network.size.variable"))
	netrec <- dfpart[dfpart$wave == wave, neid]
	netpeer <- dfpart[dfpart$wave > wave, neid]

	## If using additional population information
	if (restrict.local == TRUE) {

	  # Make sure have required information
	  stopifnot(is.null(popadjmat) == FALSE)

	  for (i in 1:nrec) {
	    for (j in 1:npeer) {
	      if (popadjmat[rec.id[i], peer.id[j]] == 0) {
	        netstageF[i, j] <- NA
	      }
	    }
	  }

	}

	## Set self-self match to NA
	netstageF[nrec+1, (npeer+1):(npeer+3)] <- NA

	## Return list with adjacency matrix, characteristics
	return(list(netstage=netstageF,
		chrec=chrec,
		chpeer=chpeer,
		netrec=netrec,
		netpeer=netpeer))

}
