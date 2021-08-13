## Generate documentation using roxygen

#' Inference for preference coefficients given observed recruitment chain, covariates, and prior values for preference coefficients.
#'
#' Function that carries out inference step for preference coefficients for RDS.  Given an observed annotated adjacency matrix, covariates, prior values, and MCMC parameters, produces estimates for \eqn{alpha} and \eqn{beta}.
#'
#' @param rdf An \code{rds.data.frame} containing observed recruitment information.  Must also have a column for covariate(s).
#' @param covnames Vector of names of covariates in \code{rdf} to be used in the model.
#' @param ma0 Vector of prior means for \eqn{alpha} and \eqn{ksi}.  Must be equal to \code{covnames+3} in length. The last threes values are for self-match.
#' @param sa0 Matrix of prior variances for \eqn{alpha} and \eqn{ksi}.
#' @param mb0 Vector of prior means for \eqn{beta} and \eqn{zeta}. Must be equal to \code{covnames+1} in length. The last value is for self-match.
#' @param sb0 Matrix of prior variances for \eqn{beta} and \eqn{ksi}.
#' @param Kcalc Desired number of iterations (MCMC parameter). This is total, so Kcalc = burn + interval * noutiter.
#' @param tallowed Number of hours allowed on Hoffman.
#' @param ncoup Number of coupons recruiters had to distribute.
#' @param restrict.local Logical: use additional information about the underlying network? Makes computation faster. If TRUE, assumes that \code{netstage} takes three values: 1 indicates a recruitment tie (and thus a tie in the underlying network); 0 indicates no recruitment tie, but a tie does exist in the underlying network; NA indicates no tie in the underlying network, and thus no possibility of recruitment.
#' @param popadjmat Population adjacency matrix.  Must be supplied if \code{restrict.local=TRUE}.
#' @param updlast Update beta self values using information in last wave? Defaults to true.
#' @param usenetsize Should network size be used as a variable to approximate? If true, prior value must be reflected in ma0, sa0, mb0, and sb0.
#'
#' @return A list containing
#'  \item{alpha}{Value of \eqn{alpha} at each out-iteration.}
#'  \item{beta}{Value of \eqn{beta} at each out-iteration.}
#'  \item{rejcounter}{Number of rejections (due to swapping problem).}
#'  \item{earlystop}{Number of iterations after which Gibbs sampler was forced to abort for time constraints.}

urInfR <- function(rdf,
  covnames,
	ma0,
	sa0,
	mb0,
	sb0,
	Kcalc,
	tallowed,
  ncoup=3,
	restrict.local=FALSE,
	popadjmat=NULL,
  updlast=TRUE,
  usenetsize=FALSE) {

  ###############################################
  ####   INITIAL SET UP
  ####        only done before first iteration
  ###############################################

  n.wav <- max(rdf$wave[rdf$wave!=1000])

  # Get amended adjacency matrix for each wave
  netstage <- list()
  charrec <- list()
  charpeer <- list()
  netrec <- list()
  netpeer <- list()

  for (i in 1:n.wav) {

    temp <- makeNetStageR(rdf,
                          covnames,
                          i-1,
                          restrict.local=restrict.local,
                          popadjmat=popadjmat,
                          usenetsize=usenetsize)

    netstage[[i]] <- temp$netstage
    charrec[[i]] <- temp$chrec
    charpeer[[i]] <- temp$chpeer
    netrec[[i]] <- temp$netrec / max(temp$netrec)
    netpeer[[i]] <- temp$netpeer / max(temp$netpeer)

  }

	## Get network characteristics
  nrec  <- list()
  npeer <- list()
  for (i in 1:n.wav) {
	  nrec[[i]]  <- nrow(netstage[[i]]) - 1
	  npeer[[i]] <- ncol(netstage[[i]]) - 3
  }

	## Draw initial alpha, beta from prior distributions
	alpha0 <- mvrnorm(n = 1,
		mu = ma0,
		Sigma = sa0)
	beta0 <- mvrnorm(n = 1,
		mu = mb0,
		Sigma = sb0)

	## Calculate initial U and V (from linear formula)

	epsilon <- list()
	gamma <- list()
	X <- list()
	Y <- list()
	Ustar <- list()
	Vstar <- list()

	for (i in 1:n.wav) {

	  epsilon[[i]] <- matrix(rnorm(nrec[[i]]*(npeer[[i]]+3)),
		  nrow=nrec[[i]],
		  ncol=npeer[[i]]+3)

  	gamma[[i]] <- matrix(rnorm((nrec[[i]]+1)*npeer[[i]], mean=0, sd=1),
  		nrow=npeer[[i]],
  		ncol=nrec[[i]]+1)

  	# cpp functions
	  X[[i]] <- makeX(charrec[[i]], charpeer[[i]], netstage[[i]])
	  Y[[i]] <- makeY(charrec[[i]], charpeer[[i]], netstage[[i]])

	  if (usenetsize == TRUE) {
	    for (rr in 1:nrec[[i]]) {
	      tmm <- X[[i]][[rr]][1:(length(ma0)-4),]
	      tmme <- X[[i]][[rr]][(length(ma0)-3):(length(ma0)-1),]
	      X[[i]][[rr]] <- rbind(tmm, c(netpeer[[i]], 0, 0, 0), tmme)
	    }
	    for (pp in 1:npeer[[i]]) {
	      tmm <- Y[[i]][[pp]][1:(length(mb0)-2),]
	      tmme <- Y[[i]][[pp]][(length(mb0)-1):(length(mb0)-1),]
	      Y[[i]][[pp]] <- rbind(tmm, c(netrec[[i]], 0), tmme)
	    }
	  }

	  Ustar[[i]] <- calcU(alpha0, X[[i]], epsilon[[i]])
	  Vstar[[i]] <- calcV(beta0, Y[[i]], gamma[[i]])

	}

	## If using additional population information, add NA to Ustar and Vstar AND to X and Y
	reslocal <- FALSE
	if (restrict.local == TRUE) {

	  for (wav in 1:n.wav) {

  	  for (i in 1:nrec[[wav]]) {
  	    for (j in 1:npeer[[wav]]) {
  	      if (is.na(netstage[[wav]][i, j])) {
  	        Ustar[[wav]][i, j] <- NA
  	        Vstar[[wav]][j, i] <- NA
  	        X[[wav]][[i]][j] <- NA
  	        Y[[wav]][[j]][i] <- NA
  	      }
  	    }
  	  }

	  }

	  reslocal <- TRUE

	}

	## Calculate alpha^(1) and beta^(1) using the first wave

	## cpp function
	alphastart <- alphaFromU(X[[1]], Ustar[[1]], ma0, sa0, mvrnorm)
	betastart  <- betaFromV(Y[[1]], Vstar[[1]], mb0, sb0, mvrnorm)

	## Iterate many times
  ## cpp function
	outList <- gibbsC(Kcalc,
		alphastart,
		betastart,
		netstage,
		nrec,
		npeer,
		X,
		Y,
		Ustar,
		Vstar,
		expand.grid,
		rtruncnorm,
		sample,
		ma0,
		sa0,
		mb0,
		sb0,
		tallowed,
		recfromUV3urR,
		reslocal,
		mvrnorm,
		updlast)

	return(outList)

}
