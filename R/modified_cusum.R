library("funtimes")
library("compiler")
library(osc)

Mfun <- function(e, k) {
  T <- length(e)
  # if (length(k) > 1)
  #   k <- unique(sort(k))
  m <- length(k)
  sume <- sum(e)
  m1 <- (sum(e[1:k[1]]) - k[1]*sume/T) / sqrt(k[1])
  mm <- (sum(e[(k[m] + 1):T]) - sume*(T - k[m])/T) / sqrt(T - k[m])
  mi <- numeric()
  if (m > 1) {
    for (i in 2:m) {
      mi[i - 1] <- (sum(e[(k[i - 1] + 1):k[i]]) - sume*(k[i] - k[i - 1])/T) / sqrt(T)
    }
  }
  sum(abs(c(m1, mi, mm)))
}

cMfun <- cmpfun(Mfun)

M1fun <- function(x, e) { # for at most one change point
  T <- length(e)
  me <- .Internal(mean(e)) # 'naked' mean function
  #me <- mean(e)
  m1 <- (sum(e[1:x]) - x*me) / sqrt(x)
  mm <- (sum(e[(x + 1):T]) - me*(T - x)) / sqrt(T - x)
  sum(abs(c(m1, mm)))
}

cM1fun <- cmpfun(M1fun)

MTfun <- function(e, m, k = NULL, x = NULL) { # ignore x, needed for sapply
  if (is.null(k)) { #explore all combinations of at-most-m change points
    T <- length(e)
    M <- K <- as.list(rep(NA, m))
    for (i in 1:m) {
      K[[i]] <- combn(T - 1, i)
      M[[i]] <- sapply(1:ncol(K[[i]]), function(x) cMfun(e, K[[i]][,x]))
    }
  } else {#explore only the pre-defined k's
    #k <- unique(sort(k))
    #m <- length(k)
    M <- K <- as.list(rep(NA, m))
    if (m == 1) {
      if (length(k) == 1) {
        K[[1]] <- k
        M[[1]] <- cM1fun(x = k, e = e)
      } else {
        K[[1]] <- matrix(k, nrow = 1)
        M[[1]] <- sapply(k, cM1fun, e = e)
        #M[[1]] <- sapply(k, function(x) cM1fun(e, x))
      }
    } else {
      K[[1]] <- matrix(k, nrow = 1)
      M[[1]] <- sapply(k, cM1fun, e = e)
      for (i in 2:m) {
        K[[i]] <- combn(k, i)
        M[[i]] <- apply(K[[i]], 2, function(x) cMfun(e, x))
        #M[[i]] <- sapply(1:ncol(K[[i]]), function(x) cMfun(e, K[[i]][,x]))
      }
    }
  }
  mhat <- which.max(sapply(M, max))
  tmp <- which.max(M[[mhat]])
  if (m == 1) {
    # K[[mhat]][,tmp]
    khat <- k[tmp]
    # was: khat <- k
  } else {
    khat <- K[[mhat]][,tmp]
  }
  MT <- M[[mhat]][tmp]
  list(MT = MT, m = mhat, k = khat)
}

cMTfun <- cmpfun(MTfun)

mcusum <- function(e, k, m = length(k),
                         B = 1000, shortboot = TRUE, ksm = FALSE,
                         ksm.arg = list(kernel = "gaussian", bw = "sj"), ...)
{
  DNAME <- deparse(substitute(e))
  T <- length(e)
  e <- e - .Internal(mean(e))
  k <- sort(unique(k))
  phi <- ARest(e, ...)
  MTobs <- cMTfun(e, k = k, m = m)
  if (length(phi) > 0) {
    e <-  as.vector(embed(e, length(phi) + 1L) %*% c(1, -phi))
    e <- e - .Internal(mean(e))
  }
  if (ksm) { #use e from a smoothed distribution of e
    ksm.arg$x <- e #append x to the arguments of the density function
    bw <- do.call(stats::density, ksm.arg) #estimate bandwidth
    bw <- bw$bw
    MTboot <- sapply(1:B, function(b)
      cMTfun(as.vector(arima.sim(T, model = list(order = c(length(phi), 0, 0), ar = phi),
                                 innov = rnorm(T, mean = sample(e, size = T, replace = TRUE), sd = bw))),
             k = k, m = m)$MT
    )
  } else {#use bootstrapped e

    if(shortboot){
      # bootstrap only for a portion of B first:

      B_part <- ceiling(B/4) # use a quarter of B, but could also be less/more
      sig <- ceiling(B/10) # portion of bootstraps that has to be bigger than original for alpha = 0.1
      thr <- sig/B_part # prior threshold to discard time series which won't reach significant threshold anymore

      MTboot <- sapply(1:B_part, function(b)
        cMTfun(as.vector(arima.sim(T, model = list(order = c(length(phi), 0, 0), ar = phi),
                                   innov = sample(e, size = T, replace = TRUE))),
               k = k, m = m)$MT
      )

      if(mean(MTboot >= MTobs$MT) < thr){
        MTboot2 <- sapply(1:(B-B_part), function(b)
          cMTfun(as.vector(arima.sim(T, model = list(order = c(length(phi), 0, 0), ar = phi),
                                     innov = sample(e, size = T, replace = TRUE))),
                 k = k, m = m)$MT
        )
        MTboot <- c(MTboot, MTboot2)
        pval <- mean(MTboot >= MTobs$MT)+1/B
      } else {
        pval <- 999
      }

    } else { # no shortening of bootstrap
      MTboot <- sapply(1:B, function(b)
        cMTfun(as.vector(arima.sim(T, model = list(order = c(length(phi), 0, 0), ar = phi),
                                   innov = sample(e, size = T, replace = TRUE))),
               k = k, m = m)$MT
      )
      pval <- mean(MTboot >= MTobs$MT)+1/B
    }


  }
  STATISTIC <- MTobs$MT
  names(STATISTIC) <- "M_T"
  PARAMETER <- length(MTobs$k)
  names(PARAMETER) <- "mhat"
  ESTIMATE <- list(length(phi), phi, MTobs$k, B)
  names(ESTIMATE) <- c("AR_order", "AR_coefficients", "khat", "B")
  structure(list(method = "Test for at-most-m changes in linear regression model",
                 data.name = DNAME,
                 statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = pval,
                 alternative = paste("at-most-", length(k), " changes exist", sep = ""),
                 estimate = ESTIMATE),
            class = "htest")

}

#' Modified CUSUM test for at-most-m change points
#'
#' This function implements the modified CUSUM by Lyubchich et al. 2020 (https://doi.org/10.1002/env.2591)
#' adjusted with a reduction in combinations of change points compared and with an adaptive
#' bootstrapping that stops if significance can't be reached anymore after a certain
#' amount of bootstraps.
#' @param x time series at a single grid cell
#' @param t a vector containing the first and last possible change point location, all locations within are compared
#' @param m the maximum number of change points allowed
#' @param B number of bootstrap samples
#' @param loc Should the location of the change point be returned?
#' @return p-value of the chosen change point location, if loc = T a vector containing the p-value and the location
#' @export mcusum_function

mcusum_function <- function(x, loc = F){
  library(funtimes)
  library(osc)
  if(any(is.na(x))) return(NA)
  d <- ts(x)
  ehat <- lm(d ~ time(d))$resid # get the residuals
  res <- mcusum(ehat, m = 1, shortboot = T, k = c(4:34), B=500)
  if(loc){
    return(c(res$p.value, res$estimate$khat))
  } else {
    return(res$p.value)
  }
}
