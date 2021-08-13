#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# Simulation to check FPR of different change point detection methods:

# Apply test n times to a new simulated random data set without any breakpoints
# and check the number of times it yielded a significant result (i.e. false positive)
# at alpha = 0.0.5
#**************************

library(tidyverse)
n <- 1000


# wbsts package ----
#*********************************************************************
library(wbsts)
d <- ts(rnorm(38,12000,1500), frequency=1, start=1)
res <- wbs.lsw(d)

hits <- 0
for (i in 1:n) {
    d <- ts(rnorm(400,12000,1500), frequency=1, start=1)
    res <- wbs.lsw(d)
    if(!is.null(res$cp.aft)){
      hits <- hits + 1
    }
  }

fpr <- hits/n
cat(paste('False positive rate for wild binary segmentation test is: ', fpr))

# Shows problems for time series short than 300:
# Error in if (cbr[i + 2] == nz) { : missing value where TRUE/FALSE needed


# strucchange package ----
#*********************************************************************
library(strucchange)

# Generalized fluctuation tests
#************
test_CPD_efp <- function(n, fun, h = NULL){
  p <- c()
  for (i in 1:n) {
    d <- ts(rnorm(38,12000,1500), frequency=1, start=1)
    temp <- efp(d ~ time(d), d, type = fun, h = h, rescale = TRUE)
    res <-  sctest(temp)
    p <- c(p, res$p.value)
  }
  fpr <- as.double(sum(p<0.05))/as.double(n)
  cat(paste('False positive rate for efp function ', fun,' is: ', fpr))
  return(p)
}

# residual based:
rec_cusum <- test_CPD_efp(n, 'Rec-CUSUM') # fpr: 0.038, 0.027, 0.023
ols_cusum <- test_CPD_efp(n, 'OLS-CUSUM') # fpr: 0, 0, 0
score_cusum <- test_CPD_efp(n, 'Score-CUSUM') #  fpr: 0.017, 0.012, 0.013
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
library(MASS)
p_LM <- c()
p_maxmo <- c()
p_bb <- c()
p_meanbb <- c()

for (i in 1:n) {
  d <- ts(rnorm(38,12000,1500), frequency=1, start=1)
  temp <- gefp(d ~ time(d), fit = rlm)
  res_LM <- sctest(temp, functional = supLM(0.05))
  p_LM <- c(p_LM, res_LM$p.value)
  res_maxmo <- sctest(temp, functional = maxMOSUM(width = 0.05))
  p_maxmo <- c(p_maxmo, res_maxmo$p.value)
  res_BB <- sctest(temp, functional =rangeBB)
  p_bb <- c(p_bb,res_BB$p.value)
  res_meanBB <- sctest(temp, functional =meanL2BB)
  p_meanbb <- c(p_meanbb,res_meanBB$p.value)

}
fpr_LM <- sum(p_LM<0.05)/n # fpr:  0.168, 0.194, 0.159
cat(paste('False positive rate for supLM is: ', fpr_LM))
fpr_maxmo <- sum(p_maxmo<0.05)/n # fpr: 0.002, 0, 0
cat(paste('False positive rate for maximum MOSUM is: ', fpr_maxmo))
fpr_bb <- sum(p_bb<0.05)/n # fpr: 0.004, 0.001, 0.003; similar for lm, glm an rlm
cat(paste('False positive rate for BB is: ', fpr_bb))
fpr_meanbb <- sum(p_meanbb<0.05)/n # fpr:  0.011, 0.005, 0.013; similar for lm, glm an rlm
cat(paste('False positive rate for meanbb is: ', fpr_meanbb))


# F Statistics (compares model fits for whole series vs. two segments; iterates over whole series as possible break points)
#************
p_sup <- c()
p_exp <- c()
p_ave <- c()
for (i in 1:n) {
  d <- ts(rnorm(38,12000,1500), frequency=1, start=1)
  temp <- Fstats(d ~ time(d))
  res_sup <- sctest(temp, type = 'supF')
  p_sup <- c(p_sup, res_sup$p.value)
  res_exp <- sctest(temp, type = 'expF')
  p_exp <- c(p_exp, res_exp$p.value)
  res_ave <- sctest(temp, type = 'aveF')
  p_ave <- c(p_ave, res_ave$p.value)
}
fpr_sup <- sum(p_sup<0.05)/n # fpr: 0.07, 0.095, 0.086, 0.0808
cat(paste('False positive rate for supF is: ', fpr_sup))
fpr_exp <- sum(p_exp<0.05)/n # fpr: 0.073, 0.108, 0.087, 0.0906
cat(paste('False positive rate for expF is: ', fpr_exp))
fpr_ave <- sum(p_ave<0.05)/n # fpr: 0.056, 0.08, 0.05, 0.0587
cat(paste('False positive rate for aveF is: ', fpr_ave))


# AMOC (changepoint package) ----
# at most one change
#*********************************************************************
library(changepoint)

fpr_mean <- c()
for (i in 1:n) {
  d <- ts(rnorm(38,12000,1500), frequency=1, start=1)
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
  d <- ts(rnorm(38,0,4), frequency=1, start=1)
  res <- pettitt.test(d)
  p <- c(p, res$p.value)
}
fpr <- sum(p<0.05)/n
cat(paste('False positive rate for pettitts test is: ', fpr))
# fpr:  0.0273, 0.0251, 0.0231

# MCUSUM (Lyubchich2020, funtimes package) with our modifications----
#*********************************************************************
devtools::load_all()
p <- c()
for (i in 1:n) {
  d <- ts(rnorm(38,0,4), frequency=1, start=1)
  #d <- 2*c(1:38) + 4 + rnorm(n=38, mean=0, sd=20)
  ehat <- lm(d ~ time(d))[["residuals"]]
  res <- mcusum_function(ehat)
  p <- c(p, res)
}
fpr <- sum(p<0.05)/n
cat(paste('False positive rate for modified cusum test is: ', fpr)) 
# 0 for random data
# 0-0.002 for random data with trend (0.002 with sd = 10, 0 for sd = 20 and sd = 30)



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
