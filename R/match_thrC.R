## Generate documentation using roxygen

#' Generate recruitment chain based on utilities.
#'
#' Generates a recruitment chain based on the given utilities \code{U} and \code{V}.  Uses the multiple matching process I developed.  This is the matching algorithm for the UNRANKED process (i.e., if a person recruits three peers, we do not make assumptions about those peers' relative utilities.)
#'
#' @param U Matrix of utilities from recruiters to peers.  It has dimension \eqn{n_r \times (n_p + n_c)}. Currently implemented for \eqn{n_c=3}.
#' @param V Matrix of utilities from peers to recruiters.  It has dimension \eqn{n_p \times (n_r + 1)}.
#' @param reslocal Logical: use additional information about the underlying network? Makes
#'   computation faster. If TRUE, assumes that \code{netstage} takes three values: 1 indicates a
#'   recruitment tie (and thus a tie in the underlying network); 0 indicates no recruitment tie,
#'   but a tie does exist in the underlying network; NA indicates no tie in the underlying network,
#'   and thus no possibility of recruitment.
#' @param siml Logical: is this a simulation (gives warning when < 3 peers).
#'
#' @return An annotated adjacency matrix with the generated recruitment chain.  It has dimension \eqn{(n_r+n_c) \times (n_p+1)} where \eqn{n_r} is the number of recruiters, \eqn{n_p} is the number of peers, and \eqn{n_c} is the number of coupons.
#'
#' @examples
#' U <- matrix(rnorm(50), nrow=5)
#' V <- matrix(rnorm(42), nrow=7)
#' recfromUV3urR(U, V)
recfromUV3urR <- function(U, V, reslocal=FALSE, siml=FALSE) {

	nrec  <- nrow(U)
	npeer <- nrow(V)

	# add in stop/error when have too few peers
	if(siml == TRUE & npeer < 3) {
	  stop("Not enough peers")
	}

	stopifnot(nrec == ncol(V)-1)
	stopifnot(npeer == ncol(U)-3)

	Urank <- matrix(nrow=nrec, ncol=npeer+3)
	Vrank <- matrix(nrow=npeer, ncol=nrec+1)

	for (i in 1:nrec) {
		Urank[i, ] <- rank(U[i, ], na.last="keep")		#1 is low
	}
	for (j in 1:npeer) {
		Vrank[j, ] <- rank(V[j, ], na.last="keep")		# 1 is low
	}

	## Call C function to generate matching from U, V
	recruitreturn <- recfromUV3urC(Urank, Vrank, reslocal)				#cpp function
	utiladj <- recruitreturn$utiladj

	rownames(utiladj) <- rownames(U)
	colnames(utiladj) <- rownames(V)

	return(utiladj)

}
