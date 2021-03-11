# simulation to check FPR of different change point detection methods:
# Apply test n times to a new simulated random data set without any breakpoints
# and check the number of times it yielded a significant result (i.e. false positive)
# at alpha = 0.0.5
#**************************
library(tidyverse)
n <- 10000

# strucchange package ----
#*********************************************************************
library(strucchange)

# Generalized fluctuation tests
#************
test_CPD_efp <- function(n, fun, h = NULL){
  p <- c()
  for (i in 1:n) {
    d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
    temp <- efp(d ~ 1, d, type = fun, h = h, rescale = TRUE)
    res <-  sctest(temp)
    p <- c(p, res$p.value)
  }
  fpr <- as.double(sum(p<0.05))/as.double(n)
  cat(paste('False positive rate for efp function ', fun,' is: ', fpr))
  return(p)
}

# residual based:
rec_cusum <- test_CPD_efp(n, 'Rec-CUSUM') # fpr: 0.0472, 0.0426, 0.0433, 0.0442
ols_cusum <- test_CPD_efp(n, 'OLS-CUSUM') # fpr: 0.0426, 0.0396, 0.0406
score_cusum <- test_CPD_efp(n, 'Score-CUSUM') #  fpr: 0.0358, 0.0401
# different window sizes (performance of certain window size is probably data specific)
h <- c(0.01, 0.05, 0.1, 0.2, 0.3)

rec_mosum <- function(h){
  rec_mosum <- test_CPD_efp(10000, 'Rec-MOSUM', h = h)
  return(rec_mosum)
}
ols_mosum <- function(h){
  ols_mosum <- test_CPD_efp(10000, 'OLS-MOSUM', h = h)
  return(rec_mosum)
}
rmosum <- purrr::map(h, rec_mosum)
# h = 0.01: 0, h = 0.05: 0.0138, h = 0.1: 0.0289, h = 0.2: 0.0346, h = 0.3: 0.0388
omosum <-  purrr::map(h, ols_mosum)
# h = 0.01: 0, h = 0.05: 0.0262, h = 0.1: 0.0337, h = 0.2: 0.0354, h = 0.3: 0.0398

# estimates based:
recur_estimates <- test_CPD_efp(n, 'RE') # fpr: 0.0386, 0.0419
mov_estimates <- test_CPD_efp(n, 'ME', h = 0.05) # fpr: 0.023, 0.0254

# Generalized M-Fluctuation Tests
#************
p_LM <- c()
p_maxmo <- c()

for (i in 1:n) {
  d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  temp <- gefp(d ~ 1, fit = lm)
  res_LM <- sctest(temp, functional = supLM(0.05))
  p_LM <- c(p_LM, res_LM$p.value)
  res_maxmo <- sctest(temp, functional = maxMOSUM(width = 0.05))
  p_maxmo <- c(p_maxmo, res_maxmo$p.value)
}
fpr_LM <- sum(p_LM<0.05)/n # fpr: 0.0412
cat(paste('False positive rate for supLM is: ', fpr_LM))
fpr_maxmo <- sum(p_maxmo<0.05)/n # fpr: 0.02
cat(paste('False positive rate for maximum MOSUM is: ', fpr_maxmo))


# F Statistics (compares model fits for whole series vs. two segments; iterates over whole series as possible break points)
#************
p_sup <- c()
p_exp <- c()
p_ave <- c()
for (i in 1:n) {
  d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  temp <- Fstats(d ~ 1)
  res_sup <- sctest(temp, type = 'supF')
  p_sup <- c(p_sup, res_sup$p.value)
  res_exp <- sctest(temp, type = 'expF')
  p_exp <- c(p_exp, res_exp$p.value)
  res_ave <- sctest(temp, type = 'aveF')
  p_ave <- c(p_ave, res_ave$p.value)
}
fpr_sup <- sum(p_sup<0.05)/n # fpr: 0.0518, 0.0458
cat(paste('False positive rate for supF is: ', fpr_sup))
fpr_exp <- sum(p_exp<0.05)/n # fpr: 0.0554, 0.0491
cat(paste('False positive rate for expF is: ', fpr_exp))
fpr_ave <- sum(p_ave<0.05)/n # fpr: 0.051, 0.0493
cat(paste('False positive rate for aveF is: ', fpr_ave))


# AMOC (changepoint package) ----
# at most one change
#*********************************************************************
library(changepoint)

fpr_mean <- c()
for (i in 1:n) {
  d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  res <- cpt.mean(d, penalty = "Asymptotic", pen.value = 0.05, method = "AMOC",
                     test.stat = "Normal", class = FALSE)
  # in documentation says that it results in NA if no changepoint but it often
  # returns 240 as changepoint location which is the last point, so only include results unlike 240
  if (res[1] != 240) fpr_mean <- c(fpr_mean, res[1])
}
cat(paste('False positive rate for AMOC mean is ', length(fpr_mean)/n))
# fpr: 0

fpr_var <- c()
for (i in 1:n) {
  d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  res <- cpt.var(d, penalty = "Asymptotic", pen.value = 0.05, method = "AMOC",
                  test.stat = "Normal", class = FALSE)
  if (res[1] != 240) fpr_var <- c(fpr_var, res[1])
}
cat(paste('False positive rate for AMOC var is ', length(fpr_var)/n))
# fpr: 0.012

fpr_meanvar <- c()
for (i in 1:n) {
  d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  res <- cpt.meanvar(d, penalty = "Asymptotic", pen.value = 0.05, method = "AMOC",
                  test.stat = "Normal", class = FALSE)
  if (res[1] != 240) fpr_meanvar <- c(fpr_meanvar, res[1])
}
cat(paste('False positive rate for AMOC meanvar is ', length(fpr_meanvar)/n))
# fpr: 1
# always detects change points at the penultimate and/or second/third time points...


# Pettitts (trend package) ----
#*********************************************************************
# non parametric pettitt test (also somehow based on wilcoxon rang thing)
library(trend)

p <- c()
for (i in 1:n) {
  d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  res <- pettitt.test(d)
  p <- c(p, res$p.value)
}
fpr <- sum(p<0.05)/n
cat(paste('False positive rate for pettitts test is: ', fpr))
# fpr: 0.0428, 0.0403

# cpm package ----
#*********************************************************************
# library(cpm)
#
# test_CPD_cpm <- function(n, fun){
#   bp <- 0
#   for (i in 1:n) {
#     d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
#     res <- detectChangePointBatch(d, cpmType = fun, alpha = 0.05)
#     if (res$changeDetected) bp <- c(bp, res$changePoint)
#   }
#   fpr <- length(bp)/n
#   cat(paste('False positive rate for cpm function ', fun,' is: ', fpr))
#   return(bp)
# }

# stud <- test_CPD_cpm(n, 'Student', trend) # 0.0008
# bart <- test_CPD_cpm(n, 'Bartlett', trend) # 0.03
# glr <- test_CPD_cpm(n, 'GLR', trend) # 0.13
# mw <- test_CPD_cpm(n, 'Mann-Whitney', trend) # 0.003
# mood <- test_CPD_cpm(n, 'Mood', trend) # 0.001
# lep <- test_CPD_cpm(n, 'Lepage', trend) # 0.0227
# ks <- test_CPD_cpm(n, 'Kolmogorov-Smirnov', trend) # 0.0002
# cvm <- test_CPD_cpm(n, 'Cramer-von-Mises', trend) # 0.0117
