#' Reingold Tilford plot with each recruitment chain placed in a separate box.
#'
#' This function also adds specificity for networks generated as part of the preferential recruitment modeling framework
#'
#' @param rdf An object of class \code{rds.data.frame}, with a column named \code{char} that contains node categories.
#' @param char Column name on \code{rdf} used to color the points on the plot.
#' @param title Title of plot.
#'
#' @return A recruitment plot.
#'
#' @examples
#' nseeds <- 3
#' N <- 2000
#' char <- matrix(sample(0:1, N, replace=TRUE), ncol=1, dimnames=list(1:N,"ch1"))
#' alpha <- 0
#' beta <- 0
#' ksi <- c(5,3,0)
#' zeta <- 1
#' waves <- 7
#'
#' a <- simrec(nseeds, N, char, alpha, beta, ksi, zeta, waves)
#' rtpbox(a$rds, char="ch1")

rtpbox <- function(rdf, char="char", title="Recruitment Plot") {

  	# Add more rds information
  	rdf$seed <- get.seed.id(rdf)
  	rdf$wave <- get.wave(rdf)

 	  # Use igraph plotting functionalities
  	tp.el <- matrix(c(get.rid(rdf), get.id(rdf)), ncol=2)

  	g <- graph_from_edgelist(tp.el)

  	lyt <- layout.reingold.tilford(g)
  	lyt <- lyt[-1, ]

  	seed.names <- unique(rdf$seed)
  	n.seeds <- length(seed.names)

  	# Adjust layout to make sure chains non-overlapping
  	for (i in 2:n.seeds) {

  	  sid.p <- seed.names[i-1]
  	  sid.c <- seed.names[i]
  	  cids.p <- which(rdf$seed == sid.p)
  	  cids.c <- which(rdf$seed == sid.c)
  	  max.x.p <- max(lyt[cids.p, 1])
  	  min.x.c <- min(lyt[cids.c, 1])

  	  to.add <- 2 - (min.x.c - max.x.p)

  	  sid.o <- seed.names[i:n.seeds]
  	  cids.o <- which(rdf$seed %in% sid.o)

  	  lyt[cids.o, 1] <- lyt[cids.o, 1] + to.add

  	}

  	# Set up plotting parameters
  	s.mins <- c()
  	s.maxes <- c()
  	seeds <- get.seed.id(rdf)[1:n.seeds]
  	for (i in 2:n.seeds) {

  		ids <- which(rdf$seed == seeds[i])
  			# Add more space between chains
  		lyt[ids, 1] <- lyt[ids, 1] + (i-1)*2

  		s.mins[i-1] <- min(lyt[ids, 1])
  		idsp <- which(rdf$seed == seeds[i-1])
  		s.maxes[i-1] <- max(lyt[idsp, 1])

  	}

  	mps <- (s.mins + s.maxes) / 2

  	n.seg <- dim(tp.el)[1] - n.seeds
  	x0 <- rep(0, n.seg)
  	y0 <- rep(0, n.seg)
  	x1 <- rep(0, n.seg)
  	y1 <- rep(0, n.seg)

  	tp.el.new <- tp.el[-c(1:n.seeds), ]

  	for (i in 1:n.seg) {

  		x0[i] <- lyt[which(get.id(rdf)==tp.el.new[i, 1]), 1]
  		y0[i] <- lyt[which(get.id(rdf)==tp.el.new[i, 1]), 2]
  		x1[i] <- lyt[which(get.id(rdf)==tp.el.new[i, 2]), 1]
  		y1[i] <- lyt[which(get.id(rdf)==tp.el.new[i, 2]), 2]

  	}

  	charcol <- which(names(rdf) == char)
  	if (length(charcol) == 0) {
  	  stop("Error: covariate column name does not exist")
  	}

  	n.cat <- length(unique(rdf[, charcol]))
  	palette(categorical_pal(n.cat))

  	n.wav <- max(rdf$wave)
  	n.rect <- (n.wav+1) %/% 10 + ifelse((n.wav+1) %% 10 > 0, 1, 0)

  	vcex <- ifelse(dim(rdf)[1] < 500, 1.5, 1)

  	par(mar=c(1, 1, 3, 1))

  	plot(lyt,
  		axes=FALSE,
  		main=title,
  		xlab="",
  		ylab="",
  		xlim=c(min(lyt[, 1])-2, max(lyt[, 1])+2),
  		ylim=c(-1, n.wav+1))

  	for (i in 1:n.rect) {
  		rect(xleft=min(lyt[, 1])-2,
  			xright=max(lyt[, 1])+2,
  			ybottom=ifelse(i==n.rect & i>1, -0.5, (n.wav-4.5) - (i-1)*5),
  			ytop=(n.wav+0.5) - (i-1)*10,
  			border=NA,
  			col="gray95")
  	}

  	segments(x0=x0, y0=y0, x1=x1, y1=y1, col="grey70", lwd=1.5)
  	points(lyt, col=as.numeric(rdf[, charcol] + 1), pch=19, cex=vcex)
  	points(lyt, cex=vcex)
  	segments(x0=mps, y0=-1, x1=mps, y1=n.wav+1, col="white", lwd=4)

  	for (i in 0:n.wav) {
  		text(min(lyt[, 1])-1.5, i, n.wav-i)
  		text(max(lyt[, 1])+1.5, i, n.wav-i)
  	}
  	text(min(lyt[,1])-3, n.wav/2, "Wave", srt=90)
  	text(max(lyt[,1])+3, n.wav/2, "Wave", srt=90)

  	par(mar=c(5.1, 4.1, 4.1, 2.1))

}
