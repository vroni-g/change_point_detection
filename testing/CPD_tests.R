# simulation to check FPR of change point detection methods

trend <- rep(sin(2*pi*c(0:11)/12),20) #+c(rep(0,12*10),rep(0.5,12*10))
n <- 5000

# BFAST ----
#*********************************************************************
library(bfast)
ti <- c()
for (i in 1:n) {
  d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
  res <- bfast(d, max.iter=10)
  if (!is.na(res$Time)) ti <- c(ti, res$Time)
}
length(ti)
if (length(ti) > 5) hist(ti, col="lightblue", br=length(ti)/2)
length(ti)/n # should be 0.05
plot(d) # plot the last t.s. as an example


# strucchange package ----
#*********************************************************************
library(strucchange)

# Generalized fluctuation tests
#************
test_CPD_efp <- function(n, fun, trend, h = NULL){
  p <- c()
  for (i in 1:n) {
    d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
    temp <- efp(d ~ 1, d, type = fun, h = h, rescale = TRUE)
    res <-  sctest(temp)
    p <- c(p, res$p.value)
  }
  fpr <- as.double(sum(p<0.05))/as.double(n)
  cat(paste('False positive rate for efp function ', fun,' is: ', fpr))
  return(p)
}

n <- 5000
# residual based:
rec_cusum <- test_CPD_efp(n, 'Rec-CUSUM', trend)
ols_cusum <- test_CPD_efp(n, 'OLS-CUSUM', trend)
rec_mosum <- test_CPD_efp(n, 'Rec-MOSUM', trend, h = 0.3) # tried different window sizes from 0.01 to 0.3 but always FPR = 0
ols_mosum <- test_CPD_efp(n, 'OLS-MOSUM', trend, h = 0.01)
# estimates based:
recur_estimates <- test_CPD_efp(n, 'RE', trend)
mov_estimates <- test_CPD_efp(n, 'ME', trend, h = 0.05)

# F Statistics (compares model fit for whole series vs. two segments; iterates over whole series as possible break points)
#************
test_CPD_fstat <- function(n){
  p_sup <- c()
  p_exp <- c()
  p_ave <- c()
  for (i in 1:n) {
    d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
    temp <- Fstats(d ~ 1)
    res_sup <- sctest(temp, type = 'supF')
    p_sup <- c(p_sup, res_sup$p.value)
    res_exp <- sctest(temp, type = 'expF')
    p_exp <- c(p_exp, res_exp$p.value)
    res_ave <- sctest(temp, type = 'aveF')
    p_ave <- c(p_ave, res_ave$p.value)
  }
  fpr_sup <- sum(p_sup<0.05)/n
  cat(paste('False positive rate for supF is: ', fpr_sup))
  fpr_exp <- sum(p_exp<0.05)/n
  cat(paste('False positive rate for expF is: ', fpr_exp))
  fpr_ave <- sum(p_ave<0.05)/n
  cat(paste('False positive rate for aveF is: ', fpr_ave))
}

fstat <- test_CPD_fstat(10000)


# AMOC (changepoint package) ----
# at most one change
#*********************************************************************
library(changepoint)

n <-10000

fpr_mean <- c()
for (i in 1:n) {
  d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
  res <- cpt.mean(d, penalty = "Asymptotic", pen.value = 0.05, method = "AMOC",
                     test.stat = "Normal", class = FALSE)
  # in documentation says that it results in NA if no changepoint but it always
  # returns 240 as changepoint location which is the last point...
  if (res[1] != 240) fpr_mean <- c(fpr_mean, res[1])
}
cat(paste('False positive rate for AMOC mean is ', length(fpr_mean)/n))

fpr_var <- c()
for (i in 1:n) {
  d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
  res <- cpt.var(d, penalty = "Asymptotic", pen.value = 0.05, method = "AMOC",
                  test.stat = "Normal", class = FALSE)
  if (res[1] != 240) fpr_var <- c(fpr_var, res[1])
}
cat(paste('False positive rate for AMOC var is ', length(fpr_var)/n))
# detects a handful of change points but only within the first and last two or three time points

fpr_meanvar <- c()
for (i in 1:n) {
  d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
  res <- cpt.meanvar(d, penalty = "Asymptotic", pen.value = 0.05, method = "AMOC",
                  test.stat = "Normal", class = FALSE)
  if (res[1] != 240) fpr_meanvar <- c(fpr_meanvar, res[1])
}
cat(paste('False positive rate for AMOC meanvar is ', length(fpr_meanvar)/n))
fpr_meanvar
# always detects change points at the penultimate and second/third time points...


# Pettitts (trend package) ----
#*********************************************************************
# non parametric pettitt test (also somehow based on wilcoxon rang thing)

library(trend)
#trend <- rep(sin(2*pi*c(0:11)/12),40) # larger time series (adjust the ts() command!)

n <- 10000
p <- c()
for (i in 1:n) {
  d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
  res <- pettitt.test(d)
  p <- c(p, res$p.value)
}
fpr <- sum(p<0.05)/n
cat(paste('False positive rate for pettitts test is: ', fpr))

