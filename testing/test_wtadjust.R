#************************************************************
# Check FPR of adjusted cluster detection algorithm with resulsts of detrended data
#************************************************************
library("devtools")
load_all()
library(tidyverse)

#results <- readRDS("testing/detrended_temp_data_Wtadjust_nperm_100.rds")
perm_results <- readRDS("detrend_temp_Wtadjust_nperm_1000.rds")


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

#************************************************************
# Check clustersizes of significant clusters with LAI data
#************************************************************
orig <- readRDS("testing/LAI_Wtadjust_nperm_2.rds")
# in former code the number and sizes of cluster of original data were not returned
# changed this and executed again with only 2 permutations to retrieve cluster sizes
res <- readRDS("testing/LAI_Wtadjust_nperm_100.rds")

clust_size <- orig$original_data$cluster.count
wt_vals <- res$original_wt

res_df <- data.frame(clust_size,wt_vals)
names(res_df) <- c("cluster_size","wt_values")

# values declared significant by tippet
tipp_95 <- res_df %>%
  dplyr::filter(wt_vals >= quantile(res$wt, probs = .95))
tipp_975 <- res_df %>%
  dplyr::filter(wt_vals >= quantile(res$wt, probs = .975))

# values declared significant by stcs
stcs_95 <- res_df %>%
  dplyr::filter(cluster_size >= quantile(res$stcs, probs = .95))
stcs_975 <- res_df %>%
  dplyr::filter(cluster_size >= quantile(res$stcs, probs = .975))

# values declared significant by tippet but not by stcs
tip_small_95 <- res_df %>%
  dplyr::filter(wt_vals >= quantile(res$wt, probs = .95) & cluster_size < quantile(res$stcs, probs = .95))
tip_small_975 <- res_df %>%
  dplyr::filter(wt_vals >= quantile(res$wt, probs = .975) & cluster_size < quantile(res$stcs, probs = .975))


res_df %>%
  ggplot(aes(x = cluster_size, y = wt_vals)) +
  labs(title = "Suprathreshold Cluster of Mann Kendall's Trend Test (alpha = 0.05)") +
  xlab('Suprathreshold Cluster Size') +
  ylab('Tippet combining function value') +
  geom_point(color = "steelblue3", size = 2, alpha = 0.6) +
  theme_bw() +
  geom_hline(yintercept = quantile(res$wt, probs = .95),col = "red") +
  geom_hline(yintercept = quantile(res$wt, probs = .975), col = "goldenrod1") +
  geom_vline(aes(xintercept = quantile(res$stcs, probs = .95), colour = "alpha = 0.05")) +
  geom_vline(aes(xintercept = quantile(res$stcs, probs = .975), colour = "alpha = 0.025")) +
  scale_colour_manual("MaxDistribution Threshold",
                      values = c("alpha = 0.05"="red", "alpha = 0.025"="goldenrod1")
  ) +
  geom_rect(mapping = aes(xmin=22500, xmax=42000, ymin=0.8, ymax=2), fill = NA, color = "black", show.legend = FALSE) +
  annotate("text", x=32000, y=1.79, label= "For maxDistribution alpha = 0.05:", size = 3, fontface = "bold") +
  annotate("text", x=32000, y=1.59, label= paste0("Significant Cluster STCS: ", nrow(stcs_95) ), size = 3) +
  annotate("text", x=32000, y=1.44, label= paste0("Significant Cluster Tippet: ", nrow(tipp_95)), size = 3) +
  annotate("text", x=32000, y=1.29, label= paste0("Significant Cluster Tippet but not STCS: ", nrow(tip_small_95)), size = 3) +
  annotate("text", x=32000, y=1.09, label= paste0("Total number of Cluster: ", nrow(res_df)), size = 3)
