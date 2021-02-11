# Script to check FPR of adjusted cluster detection algorithm
#************************************************************
library("devtools")
load_all()
library(tidyverse)

#results <- readRDS("testing/detrended_temp_data_Wtadjust_nperm_100.rds")
perm_results <- readRDS("detrended_temp_data_Wtadjust_nperm_1000.rds")


sim_length<- 1000
bootstrap_sample<- 100
fpr_length<- 100
alpha<- 0.05

fpr_sim_stcs<- vector(length = sim_length)
fpr_sim_maxt<- vector(length = sim_length)
fpr_sim_wt <- vector(length = sim_length)

for(j in 1:sim_length){
  
  fpr_stcs<- vector(length = fpr_length)
  fpr_maxt<- vector(length = fpr_length)
  fpr_wt<- vector(length = fpr_length)
  
  for (i in 1:fpr_length){
    ind<- sample(x = length(perm_results$maxT), size = bootstrap_sample, replace = TRUE)
    tmp_stcs<- perm_results$stcs[ind]
    tmp_maxt<- perm_results$maxT[ind]
    tmp_wt<- perm_results$wt[ind]
    
    # gets the threshold for the current sample
    q_thr_stcs<- quantile(tmp_stcs, probs = 1-alpha, names = FALSE)
    q_thr_maxt<- quantile(tmp_maxt, probs = 1-alpha, names = FALSE)
    q_thr_wt<- quantile(tmp_wt, probs = 1-alpha, names = FALSE)
    
    # retrieve false positives, i.e. values above current sample based threshold
    fpr_stcs[i]<- tmp_stcs[length(tmp_stcs)] > q_thr_stcs
    fpr_maxt[i]<- tmp_maxt[length(tmp_maxt)] > q_thr_maxt
    fpr_wt[i] <- tmp_wt[length(tmp_wt)] > q_thr_wt
    
    
  }
  fpr_sim_stcs[j]<- sum(fpr_stcs)/length(fpr_stcs)
  fpr_sim_maxt[j]<- sum(fpr_maxt)/length(fpr_maxt)
  fpr_sim_wt[j]<- sum(fpr_wt)/length(fpr_wt)
  
}

par(mfrow = c(1, 3))
hist(fpr_sim_stcs)
hist(fpr_sim_maxt)
hist(fpr_sim_wt)
par(mfrow = c(1, 1))

summary(fpr_sim_stcs)
summary(fpr_sim_maxt)
summary(fpr_sim_wt)
