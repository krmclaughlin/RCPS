#' Reingold Tilford plot with each recruitment chain placed in a separate box VERTICALLY.
#'
#' This function also adds specificity for networks generated as part of the preferential recruitment modeling framework
#'
#' @param rdf An object of class \code{rds.data.frame}, with a column named \code{char} that contains node categories.
#' @param char Column name on \code{rdf} used to color the points on the plot.
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
#' rtpboxvert(a$rds, char="ch1")

rtpboxvert <- function(rdf, char="char") {

  # Add more rds information
  rdf$seed <- get.seed.id(rdf)
  rdf$wave <- get.wave(rdf)

  maxnw <- max(rdf$wave) + 2
  n.wav <- max(rdf$wave)

  charcol <- which(names(rdf) == char)
  if (length(charcol) == 0) {
    stop("Error: covariate column name does not exist")
  }

  n.cat <- length(unique(rdf[, charcol]))
  palette(categorical_pal(n.cat))

  # Use igraph plotting functionalities
  rdft <- rdf[rdf$seed == unique(rdf$seed)[1], ]
  tp.el <- matrix(c(get.rid(rdft), get.id(rdft)), ncol=2)
  g <- graph_from_edgelist(tp.el)
  lyt <- layout.reingold.tilford(g)
  lyt <- lyt[-1, ]
  n.seg <- dim(tp.el)[1] - 1
  x0 <- rep(0, n.seg)
  y0 <- rep(0, n.seg)
  x1 <- rep(0, n.seg)
  y1 <- rep(0, n.seg)
  tp.el.new <- tp.el[-1, ]
  cc <- as.numeric(rdft[, charcol]) + 1
  for (i in 1:n.seg) {
    x0[i] <- lyt[which(get.id(rdft)==tp.el.new[i, 1]), 1]
    y0[i] <- lyt[which(get.id(rdft)==tp.el.new[i, 1]), 2]
    x1[i] <- lyt[which(get.id(rdft)==tp.el.new[i, 2]), 1]
    y1[i] <- lyt[which(get.id(rdft)==tp.el.new[i, 2]), 2]
  }
  for (us in 2:length(unique(rdf$seed))) {
    # Set up nodes
    rdft <- rdf[rdf$seed == unique(rdf$seed)[us], ]
    tp.el <- matrix(c(get.rid(rdft), get.id(rdft)), ncol=2)
    g <- graph_from_edgelist(tp.el)
    lytemp <- layout.reingold.tilford(g)
    lytemp <- lytemp[-1, ]
    lytemp[, 2] <- lytemp[, 2] + (us-1)*maxnw

    # Set up edges
    n.seg <- dim(tp.el)[1] - 1
    x0temp <- rep(0, n.seg)
    y0temp <- rep(0, n.seg)
    x1temp <- rep(0, n.seg)
    y1temp <- rep(0, n.seg)
    tp.el.new <- tp.el[-1, ]
    for (i in 1:n.seg) {
      x0temp[i] <- lytemp[which(get.id(rdft)==tp.el.new[i, 1]), 1]
      y0temp[i] <- lytemp[which(get.id(rdft)==tp.el.new[i, 1]), 2]
      x1temp[i] <- lytemp[which(get.id(rdft)==tp.el.new[i, 2]), 1]
      y1temp[i] <- lytemp[which(get.id(rdft)==tp.el.new[i, 2]), 2]
    }
    cctemp <- as.numeric(rdft[, charcol]) + 1
    lyt <- rbind(lyt, lytemp)
    x0 <- c(x0, x0temp)
    y0 <- c(y0, y0temp)
    x1 <- c(x1, x1temp)
    y1 <- c(y1, y1temp)
    cc <- c(cc, cctemp)
  }

  vcex <- ifelse(dim(rdf)[1] < 500, 1.5, 1)

  par(mar=c(1, 1, 3, 1))

  plot(lyt,
       axes=FALSE,
       main="Recruitment Plot",
       xlab="",
       ylab="",
       xlim=c(min(lyt[, 1])-2, max(lyt[, 1])+2))

  segments(x0=x0, y0=y0, x1=x1, y1=y1, col="grey70", lwd=1.5)
  points(lyt, col=cc, pch=19, cex=vcex)
  points(lyt, cex=vcex)

  for (nsee in 1:length(unique(rdf$seed))) {
    for (i in 0:n.wav) {
      text(min(lyt[, 1])-1.5, i+(nsee-1)*maxnw, n.wav-i)
      text(max(lyt[, 1])+1.5, i+(nsee-1)*maxnw, n.wav-i)
    }
    text(min(lyt[,1])-4, n.wav/2+(nsee-1)*maxnw, "Wave", srt=90)
    text(max(lyt[,1])+4, n.wav/2+(nsee-1)*maxnw, "Wave", srt=90)
  }

  par(mar=c(5.1, 4.1, 4.1, 2.1))

}
