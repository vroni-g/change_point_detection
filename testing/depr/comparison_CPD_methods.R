# simulation to check FPR of different change point detection methods
#**************************
library(tidyverse)
library(TSstudio)


trend <- rep(sin(2*pi*c(0:11)/12),200) # trend for annual means

annual_means <- function(trend){
  d <- ts(rnorm(12*200,0,.2)+trend, frequency=12, start=c(1800,1))
  means <- d %>%
    ts_reshape(type='long') %>%
    group_by(year) %>%
    summarise(annual_mean = mean(value))
  m <- ts(means$annual_mean,frequency = 1, start(1800))
  return(m)
}

n <- 10000

# BFAST ----
#*********************************************************************
# we need seasonal data for BFAST
# trend <- rep(sin(2*pi*c(0:11)/12),20)
# library(bfast)
# ti <- c()
# for (i in 1:n) {
#   d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
#   res <- bfast(d, max.iter=10)
#   if (!is.na(res$Time)) ti <- c(ti, res$Time)
# }
# #length(ti)
# #if (length(ti) > 5) hist(ti, col="lightblue", br=length(ti)/2)
# length(ti)/n # is 0.0052


# strucchange package ----
#*********************************************************************
library(strucchange)

# Generalized fluctuation tests
#************
test_CPD_efp <- function(n, fun, trend, h = NULL){
  p <- c()
  for (i in 1:n) {
    #d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
    d <- annual_means(trend)
    #d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
    temp <- efp(d ~ 1, d, type = fun, h = h, rescale = TRUE)
    res <-  sctest(temp)
    p <- c(p, res$p.value)
  }
  fpr <- as.double(sum(p<0.05))/as.double(n)
  cat(paste('False positive rate for efp function ', fun,' is: ', fpr))
  return(p)
}

n <- 10000
# residual based:
rec_cusum <- test_CPD_efp(n, 'Rec-CUSUM', trend) # annual mean: 0.0404, 0.0457; random data (no season): 0.0472, 0.0426
ols_cusum <- test_CPD_efp(n, 'OLS-CUSUM', trend) # annual mean: 0.037, 0.0363; random data (no season): 0.0426, 0.0396
score_cusum <- test_CPD_efp(n, 'Score-CUSUM', trend) # annual mean: 0.0364, random data (no season): 0.0358
# tried different window sizes from 0.01 to 0.3 but always FPR = 0
rec_mosum <- test_CPD_efp(n, 'Rec-MOSUM', trend, h = 0.01) # annual mean: 0, random data (no season): 0
ols_mosum <- test_CPD_efp(n, 'OLS-MOSUM', trend, h = 0.01) # annual mean: 0, random data (no season):0
# estimates based:
recur_estimates <- test_CPD_efp(n, 'RE', trend) # annual mean: 0.0406, random data (no season): 0.0386
mov_estimates <- test_CPD_efp(n, 'ME', trend, h = 0.05) # annual mean: 0.0192, random data (no season): 0.023

# all result in FPR = 0...

# Generalized M-Fluctuation Tests
#************
p_LM <- c()
p_maxmo <- c()
n <- 10000
for (i in 1:n) {
  #d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
  #d <- annual_means(trend)
  d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  temp <- gefp(d ~ 1, fit = lm)
  res_LM <- sctest(temp, functional = supLM(0.05))
  p_LM <- c(p_LM, res_LM$p.value)
  res_maxmo <- sctest(temp, functional = maxMOSUM(width = 0.05))
  p_maxmo <- c(p_maxmo, res_maxmo$p.value)
}
fpr_LM <- sum(p_LM<0.05)/n # annual means: 0.0426, random data (no season): 0.0412
cat(paste('False positive rate for supLM is: ', fpr_LM))
fpr_maxmo <- sum(p_maxmo<0.05)/n # annual means: 0.0236, random data (no season): 0.02
cat(paste('False positive rate for maximum MOSUM is: ', fpr_maxmo))


# F Statistics (compares model fits for whole series vs. two segments; iterates over whole series as possible break points)
#************
p_sup <- c()
p_exp <- c()
p_ave <- c()
for (i in 1:n) {
  #d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
  #d <- annual_means(trend)
  d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  temp <- Fstats(d ~ 1)
  res_sup <- sctest(temp, type = 'supF')
  p_sup <- c(p_sup, res_sup$p.value)
  res_exp <- sctest(temp, type = 'expF')
  p_exp <- c(p_exp, res_exp$p.value)
  res_ave <- sctest(temp, type = 'aveF')
  p_ave <- c(p_ave, res_ave$p.value)
}
fpr_sup <- sum(p_sup<0.05)/n # is 0 # annual mean: 0.0514, random data (no season):0.0518
cat(paste('False positive rate for supF is: ', fpr_sup))
fpr_exp <- sum(p_exp<0.05)/n # is 0 # annual mean: 0.0562, random data (no season):0.0554
cat(paste('False positive rate for expF is: ', fpr_exp))
fpr_ave <- sum(p_ave<0.05)/n # is 0 # annual mean: 0.0544, random data (no season):0.051
cat(paste('False positive rate for aveF is: ', fpr_ave))


# AMOC (changepoint package) ----
# at most one change
#*********************************************************************
library(changepoint)

n <-10000

fpr_mean <- c()
for (i in 1:n) {
  #d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
  d <- annual_means(trend)
  #d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  res <- cpt.mean(d, penalty = "Asymptotic", pen.value = 0.05, method = "AMOC",
                     test.stat = "Normal", class = FALSE)
  # in documentation says that it results in NA if no changepoint but it often
  # returns 240 as changepoint location which is the last point, so only include results unlike 240
  #if (res[1] != 200) fpr_mean <- c(fpr_mean, res[1]) # 200 for annual means
  if (res[1] != 240) fpr_mean <- c(fpr_mean, res[1])
}
cat(paste('False positive rate for AMOC mean is ', length(fpr_mean)/n)) # fpr is 0
# annual mean: 0, random data (no season):

fpr_var <- c()
for (i in 1:n) {
  #d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
  d <- annual_means(trend)
  #d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  res <- cpt.var(d, penalty = "Asymptotic", pen.value = 0.05, method = "AMOC",
                  test.stat = "Normal", class = FALSE)
  #if (res[1] != 200) fpr_var <- c(fpr_var, res[1]) # 200 for annual means
  if (res[1] != 240) fpr_var <- c(fpr_var, res[1])
}
cat(paste('False positive rate for AMOC var is ', length(fpr_var)/n)) # fpr is 0.001
# annual mean: 0.0136, random data (no season):

fpr_meanvar <- c()
for (i in 1:n) {
  #d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
  d <- annual_means(trend)
  #d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  res <- cpt.meanvar(d, penalty = "Asymptotic", pen.value = 0.05, method = "AMOC",
                  test.stat = "Normal", class = FALSE)
  #if (res[1] != 200) fpr_meanvar <- c(fpr_meanvar, res[1]) # 200 for annual means
  if (res[1] != 240) fpr_meanvar <- c(fpr_meanvar, res[1])
}
cat(paste('False positive rate for AMOC meanvar is ', length(fpr_meanvar)/n)) # fpr is 1
fpr_meanvar
# always detects change points at the penultimate and/or second/third time points...
# annual mean: 0.7334, random data (no season):

# Pettitts (trend package) ----
#*********************************************************************
# non parametric pettitt test (also somehow based on wilcoxon rang thing)
library(trend)
#trend <- rep(sin(2*pi*c(0:11)/12),40) # larger time series (adjust the ts() command!)
n <- 10000
p <- c()
for (i in 1:n) {
  #d <- ts(rnorm(12*20,0,.2)+trend, frequency=12, start=c(2000,1))
  #d <- annual_means(trend)
  d <- ts(rnorm(12*20,0,.2), frequency=12, start=c(2000,1))
  res <- pettitt.test(d)
  p <- c(p, res$p.value)
}
fpr <- sum(p<0.05)/n
cat(paste('False positive rate for pettitts test is: ', fpr)) # is 0
# annual mean: 0.0404, random data (no season): 0.0428

# cpm package ----
#*********************************************************************
# library(cpm)
#
# test_CPD_cpm <- function(n, fun, trend){
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
#
# n <- 10000
# stud <- test_CPD_cpm(n, 'Student', trend) # 0.0008
# bart <- test_CPD_cpm(n, 'Bartlett', trend) # 0.03
# glr <- test_CPD_cpm(n, 'GLR', trend) # 0.13
# mw <- test_CPD_cpm(n, 'Mann-Whitney', trend) # 0.003
# mood <- test_CPD_cpm(n, 'Mood', trend) # 0.001
# lep <- test_CPD_cpm(n, 'Lepage', trend) # 0.0227
# ks <- test_CPD_cpm(n, 'Kolmogorov-Smirnov', trend) # 0.0002
# cvm <- test_CPD_cpm(n, 'Cramer-von-Mises', trend) # 0.0117
