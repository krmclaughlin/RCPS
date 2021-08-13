#' Methods to produce plots and diagnostics for posterior distribution fits.
#'
#' Allows for the examination of MCMC chains, posterior predictive distributions, and correlation plots between alpha and beta.
#'
#' @param out The outcome of a \code{urInfR} function call.
#' @param rdsout The outcome of a \code{simrec} function call, including the RDS data frame, original alpha and beta, and sample size(s).
#' @param burn Desired burn-in period.
#' @param int Desired thinning interval.
#' @param covnames Vector of covariate names (not including self).
#' @param N Population size.
#' @param type The desired type of outcome plot. Defaults to density. Valid types are \code{dens}, \code{mcmc}, and \code{cor}.
#' @param parallel Logical: were multiple (parallel) simulations and inferences performed?
#' @param nobself Logical: do not allow peers to reject matches. Defaults to false (peers can reject matches). Not currently used.
#'
#' @return The plot specified by \code{type}.
postplot <- function(out,
                     rdsout,
                     burn,
                     int,
                     covnames,
                     N,
                     type=c("dens"),
                     parallel=FALSE,
                     nobself=FALSE) {

  if (parallel == FALSE) {
    nwav <- nrow(out$alpha[[1]])      # number of waves
    noutiter <- 0
    for (i in 1:length(out$alpha)) {
      if (is.null(out$alpha[[i]]) == FALSE) {
        noutiter = noutiter + 1       # number of iterations
      }
    }
    nvar <- ncol(out$alpha[[1]]) - 3      # number of variables (excluding self)

    talpha <- rdsout$alphaorig
    if (nobself == FALSE) {
      tbeta <- rdsout$betaorig
    } else {
      tbeta <- rdsout$betaorig[1:nvar]
    }
    n <- nrow(rdsout$rds)             # sample size
  } else {
    nvar <- ncol(out[[1]][[1]][[1]]) - 1  # number of variables (excluding self)
    nchn <- length(out)               # number of chains

    nwav <- c()
    noutiter <- c()
    n <- c()

    for (chn in 1:nchn) {
      nwav[chn] <- nrow(out[[chn]][[1]][[1]])  # number of waves
      noutiter[chn] <- 0
      for (i in 1:length(out$alpha)) {
        if (is.null(out[[chn]][[1]][[i]]) == FALSE) {
          noutiter[chn] = noutiter[chn] + 1       # number of iterations
        }
      }
      noutiter[chn] <- length(out[[chn]][[1]]) # number of iterations
      n[chn] <- nrow(rdsout[[chn]]$rds)        # sample size
    }

    talpha <- rdsout[[1]]$alphaorig
    if (nobself == FALSE) {
      tbeta <- rdsout[[1]]$betaorig
    } else {
      tbeta <- rdsout[[1]]$betaorig[1:nvar]
    }

  }

  stopifnot(length(covnames) == nvar)

  if (parallel == FALSE) {

    ##################################################
    # Apply specified burn-in and interval
    ##################################################
    aout <- matrix(nrow=noutiter, ncol=nvar + 3)
    if (nobself == FALSE) {
      bout <- matrix(nrow=noutiter, ncol=nvar + 1)
      #for each iteration, sample a random wave to use alpha, beta value
      for (i in 1:noutiter) {
        wav <- sample(1:nwav, 1)
        aout[i, ] <- out$alpha[[i]][wav, ]
        bout[i, ] <- out$beta[[i]][wav, ]
      }
    } else {
      bout <- matrix(nrow=noutiter, ncol=nvar)
      #for each iteration, sample a random wave to use alpha, beta value
      for (i in 1:noutiter) {
        wav <- sample(1:nwav, 1)
        aout[i, ] <- out$alpha[[i]][wav, ]
        bout[i, ] <- out$beta[[i]][wav, 1:nvar]
      }
    }



    aoutbi <- aout[seq(burn+1, noutiter, int), ]
    boutbi <- bout[seq(burn+1, noutiter, int), ]

  } else {

    iterout <- c()

    for (chn in 1:nchn) {

      aout <- matrix(nrow=noutiter[chn], ncol=nvar + 3)
      if (nobself == FALSE) {
        bout <- matrix(nrow=noutiter, ncol=nvar + 1)
        for (i in 1:noutiter[chn]) {
          wav <- sample(1:nwav[chn], 1)
          aout[i, ] <- out[[chn]][[1]][[i]][wav, ]
          bout[i, ] <- out[[chn]][[2]][[i]][wav, ]
        }
      } else {
        bout <- matrix(nrow=noutiter, ncol=nvar)
        for (i in 1:noutiter[chn]) {
          wav <- sample(1:nwav[chn], 1)
          aout[i, ] <- out[[chn]][[1]][[i]][wav, ]
          bout[i, ] <- out[[chn]][[2]][[i]][wav, 1:nvar]
        }
      }



      aoutbitemp <- aout[seq(burn+1, noutiter[chn], int), ]
      boutbitemp <- bout[seq(burn+1, noutiter[chn], int), ]

      iterout[chn] <- nrow(aoutbitemp)

      if (chn == 1) {
        aoutbi <- aoutbitemp
        boutbi <- boutbitemp
      } else {
        aoutbi <- rbind(aoutbi, aoutbitemp)
        boutbi <- rbind(boutbi, boutbitemp)
      }

    }

  }

  iterf <- nrow(aoutbi)

  if (is.null(talpha)) {
    talpha <- colMeans(aoutbi)
    talpha[(nvar+1):(nvar+3)] <- sort(talpha[(nvar+1):(nvar+3)])
  }
  if (is.null(tbeta)) {
    tbeta <- colMeans(boutbi)
  }

  ##################################################
  # Plot of posterior distribution densities
  ##################################################

  if ("dens" %in% type) {

    taord <- c(talpha[1:nvar],
               talpha[length(talpha):(length(talpha)-2)])

    Value <- c(aoutbi, boutbi)
    if (nobself == FALSE) {
      Covariate <- rep(c(covnames, "S3", "S2", "S1", covnames, "S1"),
                     each=iterf)
      Preference <- c(rep("Alpha", iterf*(nvar+3)),
                    rep("Beta", iterf*(nvar+1)))
      Iteration <- rep(1:iterf,
                     times=(nvar+3 + nvar+1))
    } else {
      Covariate <- rep(c(covnames, "S3", "S2", "S1", covnames),
                       each=iterf)
      Preference <- c(rep("Alpha", iterf*(nvar+3)),
                      rep("Beta", iterf*(nvar)))
      Iteration <- rep(1:iterf,
                       times=(nvar+3 + nvar))
    }
    TrueValue <- c(rep(taord, each=iterf),
                   rep(tbeta, each=iterf))
    dfplot <- data.frame(Value, Covariate, Preference, Iteration, TrueValue)

    print(ggplot(dfplot, aes(x = Value), environment = environment()) +
      geom_density(fill="gray") +
      facet_grid(Preference ~ Covariate) +
      ggtitle("Posterior Distributions") +
      theme(plot.title = element_text(face="bold")) +
      ylab("Density") +
      geom_vline(aes(xintercept = TrueValue), colour="red"))

  }

  ##################################################
  # Plot of MCMC draws from posterior distribution
  ##################################################

  if ("mcmc" %in% type) {

    taord <- c(talpha[1:nvar],
               talpha[length(talpha):(length(talpha)-2)])

    Value <- c(aoutbi, boutbi)
    if (nobself == FALSE) {
      Covariate <- rep(c(covnames, "S3", "S2", "S1", covnames, "S1"),
                       each=iterf)
      Preference <- c(rep("Alpha", iterf*(nvar+3)),
                      rep("Beta", iterf*(nvar+1)))
      Iteration <- rep(1:iterf,
                       times=(nvar+3 + nvar+1))
    } else {
      Covariate <- rep(c(covnames, "S3", "S2", "S1", covnames),
                       each=iterf)
      Preference <- c(rep("Alpha", iterf*(nvar+3)),
                      rep("Beta", iterf*(nvar)))
      Iteration <- rep(1:iterf,
                       times=(nvar+3 + nvar))
    }
    TrueValue <- c(rep(taord, each=iterf),
                   rep(tbeta, each=iterf))
    dfplot <- data.frame(Value, Covariate, Preference, Iteration, TrueValue)

    if (parallel == FALSE) {

    print(ggplot(dfplot, aes(x = Iteration, y = Value), environment = environment()) +
      geom_line() +
      facet_grid(Preference ~ Covariate) +
      geom_hline(aes(yintercept = TrueValue), colour="red") +
      ggtitle("MCMC draws from Posterior Distribution") +
      theme(plot.title = element_text(face="bold")))

    } else {

      vlin <- c()
      sm <- 0
      for (ad in 1:(nchn-1)) {
        sm <- iterout[ad] + sm
        vlin[ad] <- sm
      }

      print(ggplot(dfplot, aes(x = Iteration, y = Value), environment = environment()) +
              geom_line() +
              facet_grid(Preference ~ Covariate) +
              geom_hline(aes(yintercept = TrueValue), colour="red") +
              ggtitle("MCMC draws from Posterior Distribution") +
              theme(plot.title = element_text(face="bold")) +
              geom_vline(xintercept = vlin, colour="blue"))


    }

  }
}
