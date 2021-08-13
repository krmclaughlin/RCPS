## Use roxygen to make documentation for package

#' Simulate recruitment chains given parameters.
#'
#' @param nseeds Number of seeds.
#' @param N Population size (including seeds).
#' @param ch Matrix of covariate(s).  Must have \code{N} rows.
#' @param alpha Preference coefficients for recruiters to peers. The last value is for self-match.
#' @param beta Preference coefficients for peers to recruiters. The last value is for self-match.
#' @param ksi Number of recruits parameter.
#' @param zeta Peer self-match parameter.
#' @param waves Number of waves.
#' @param ncoup Maximum number of coupons; set to 3 by default.
#' @param restrict.local Logical: use additional information about the underlying network? Makes computation faster. Must also specify \code{popadjmat}.
#' @param popadjmat Population adjacency matrix.  Must be supplied if \code{restrict.local=TRUE}.
#' @param specifyseeds Logical: use specific nodes as the seeds? Used for simulating inclusion probabilities after inference.  Must also specify \code{seedID}.
#' @param seedID ID of specific seeds to use if \code{specifyseeds=TRUE}.
#' @param targetn Logical: specify target sample size? FALSE by default. If true need to specify value to use n instead of number of waves as stopping point.
#' @param tss Target sample size
#' @param usenetsize Logical: use network size variable in model?
#' @param netsiz Network size of each member of the population. These will then be divided by the max within each set of recruiters and peers.
#' @param netr Preference coefficient recruiters have for large network sizes.
#' @param netp Preference coefficient peers have for large network sizes.
#' @param outcome Outcome measure of interest (e.g., HIV, HCV).
#'
#' @return A list consisting of
#'  \item{net}{network object for the full recruitment chain}
#'  \item{rds}{rds.data.frame object for the full recruitment chain}
#'  \item{rdsfull}{rds.data.frame object including population members not in the sample. For simulation studies. These members of pop not in the sample have recruiter value NA, wave value 1000 (larger than would be encountered).}
#'  \item{seedsrec}{indicator of whether or not the seeds recruited anyone}
#'  \item{alphaorig}{original value of alpha used for simulation}
#'  \item{betaorig}{original value of beta used for simulation}
#'
#' @examples
#' nseeds <- 2
#' N <- 500
#' char <- matrix(sample(0:1, N, replace=TRUE), ncol=1, dimnames=list(1:N,"ch1"))
#' alpha <- 0
#' beta <- 0
#' ksi <- c(1,1,1)
#' zeta <- 0.8
#' waves <- 4
#'
#' simrec(nseeds, N, char, alpha, beta, ksi, zeta, waves)
simrec <- function(nseeds,
	N,
	ch,
	alpha,
	beta,
	ksi,
	zeta,
	waves,
	ncoup=3,
	restrict.local=FALSE,
	popadjmat=NULL,
	specifyseeds=FALSE,
	seedID=NULL,
	targetn=FALSE,
	tss=NULL,
	usenetsize=FALSE,
	netsiz=NULL,
	netr=NULL,
	netp=NULL,
	outcome=NULL) {

	## Wave 1: seeds sample from whole population\seeds

	nrec  <- nseeds
	npeer <- N - nrec

	index <- 1:N
	if (specifyseeds == FALSE) {
	  sindex <- sample(index, size=nrec, replace=FALSE)
	} else {
	  sindex <- seedID
	}
	pindex <- index[-sindex]

	crec  <- ch[sindex, , drop=FALSE]
	cpeer <- ch[pindex, , drop=FALSE]

	if (usenetsize == TRUE) {
	  recnet <- netsiz[sindex]
	  recnet <- recnet / max(recnet)
	  peernet <- netsiz[pindex]
	  peernet <- peernet / max(peernet)
	}

	nchar <- ncol(ch)

	## Perfect recruitment warning

	max.particip <- nseeds * (1 + sum(ncoup ^ c(1:waves)))
	if (max.particip > N) {
		print("WARNING: perfect recruitment would exhaust population")
	}

	## Generate utilities

	epsilon <- matrix(rnorm(nrec * (npeer+3), mean=0, sd=1),
		nrow=nrec,
		ncol=npeer+3)
	gamma   <- matrix(rnorm( npeer*(nrec+1), mean=0, sd=1),
		nrow=npeer,
		ncol=nrec+1)

	U <- makeU(crec, cpeer, alpha, ksi, epsilon)		#cpp function
	rownames(U) <- sindex
	colnames(U) <- c(pindex, "s1", "s2", "s3")

	V <- makeV(crec, cpeer, beta, zeta, gamma)		#cpp function
	rownames(V) <- pindex
	colnames(V) <- c(sindex, "0")

	if (usenetsize == TRUE) {
  	for (i in 1:nrec) {
  	  U[i,1:npeer] <- U[i,1:npeer] + netr*peernet
  	}
  	for (j in 1:npeer) {
  	  V[j,1:nrec] <- V[j,1:nrec] + netp*recnet
  	}
	}

	## If using additional population information
	if (restrict.local == TRUE) {

	  #Make sure dimensionality works out
	  stopifnot(is.null(popadjmat) == FALSE)
	  stopifnot(nrow(popadjmat) == N)

	  for (i in 1:nrec) {
	    ival <- as.integer(rownames(U)[i])
	    for (j in 1:npeer) {
	      jval <- as.integer(rownames(V)[j])
	      if (popadjmat[ival, jval] == 0) {
	        U[i, j] <- NA
	        V[j, i] <- NA
	      }
	    }
	  }

	}


	## Recruitment

	utiladjmat <- recfromUV3urR(U, V, siml=TRUE)			#R function that calls a cpp function
	rownames(utiladjmat) <- sindex
	colnames(utiladjmat) <- pindex

	elbuild <- makeEL(utiladjmat, sindex, pindex, 1)	#cpp function

	seedsrec <- 1
	if (nrow(elbuild) == 0) {
	  print("Error: No seeds recruited. Relax self cost.")
	  seedsrec <- 0
	}

	if (seedsrec == 1) {

  	## Multiple waves of recruitment (use if target n not specified)
  	if (targetn == FALSE) {

    	if (waves > 1) {

    		for (w in 2:waves) {

    			eltemp <- RaddEL(w, N, alpha, beta, ksi, zeta, elbuild, pindex, ch, usenetsize, netsiz, netr, netp)
    			if (nrow(eltemp) > 0) {
    			  elbuild <- rbind(elbuild, eltemp)
    			} else {
    			  print("Break: recruitment stopped before desired number of waves could be attained")
    			  break
    			}

    		}

    	}

  	} else {
  	  # Use n as target sample size

  	  current.n <- nrow(elbuild) + nseeds
  	  w <- 2
  	  hadbreak <- 0
  	  while (current.n < tss) {

  	    eltemp <- RaddEL(w, N, alpha, beta, ksi, zeta, elbuild, pindex, ch, usenetsize, netsiz, netr, netp)
  	    if (nrow(eltemp) > 0) {
  	      elbuild <- rbind(elbuild, eltemp)
  	      w <- w + 1
  	      current.n <- nrow(elbuild) + nseeds
  	    } else {
  	      hadbreak <- 1
  	      print("Break: recruitment stopped before desired number of waves could be attained")
  	      break
  	    }

  	  }

  	  over <- nrow(elbuild) + nseeds - tss
  	  if (over > 0 & hadbreak == 0) {
  	    #randomly sample #over dyads from last wave, remove
  	    nlw <- length(which(elbuild$wave == max(elbuild$wave)))
  	    if(over > nlw) {
  	      print(over)
  	      print(nlw)
  	    }

  	    rms <- sample(c(nrow(elbuild)-nlw+1):nrow(elbuild), over, replace=FALSE)
  	    elbuild <- elbuild[-rms, ]
  	  }

  	}


  	## Make into network object
  	net.t <- network(elbuild)
  	for (i in 1:nchar) {
  	  net.t %v% paste("char", i, sep="") <- ch[, i]
  	}

  	sdf <- data.frame(recruiter = rep(0, nseeds),
  		peer = sindex,
  		wave = rep(0, nseeds))
  	elb <- rbind(sdf, elbuild)
  	if (usenetsize == FALSE) {
  	  elb$net.size <- sample(c(49, 51),
  		dim(elb)[1],
  		replace=TRUE)
  	} else {
  	  elb$net.size <- netsiz[elb$peer]
  	}
  	elb <- cbind(elb, ch[elb$peer, , drop=FALSE])
  	if (is.null(outcome) == FALSE) {
  	  elb <- cbind(elb, outcome[elb$peer, , drop=FALSE])
  	}

  	rdf <- as.rds.data.frame(elb,
  		id = "peer",
  		recruiter.id = "recruiter",
  		network.size = "net.size",
  		population.size= N,
  		max.coupons = ncoup)
  	rdf$peer <- as.numeric(rdf$peer)

  	ns <- nrow(rdf)
  	peerex <- c(1:N)[-c(rdf$peer)]
  	char1ex <- ch[peerex, ]
  	if (usenetsize == FALSE) {
  	  nsex <- sample(c(49, 51),N-ns,replace=TRUE)
  	} else {
  	  nsex <- netsiz[peerex]
  	}
  	extra <- data.frame(recruiter = rep(NA, N - ns),
 	                    peer = peerex,
  	                  wave = rep(1000, N - ns),
 	                    net.size = nsex,
	                    char1ex)
  	if (is.null(outcome) == FALSE) {
  	  outex <- outcome[peerex, ]
  	  extra <- cbind(extra, outex)
  	}

  	names(extra) <- names(rdf)

 	rdffull <- rbind(rdf, extra)


	} else {
	  # if seeds don't recruit anyone

	  net.t <- NULL
	  rdf <- NULL
	  rdffull <- NULL

	}

	aret <- c(alpha, ksi)
	bret <- c(beta, zeta)

	return(list(net = net.t,
		rds = rdf,
		rdsfull = rdffull,
		seedsrec = seedsrec,
		alphaorig = aret,
		betaorig = bret))

}
