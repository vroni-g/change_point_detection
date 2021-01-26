# SCRIPT FOR TESTING FUNCTIONS
library("devtools")
load_all()

# ADDITIONAL NEEDED LIBRARIES
library(tidyverse) # for pipe: %>%
library(fields) # for nicer basic plotting: image.plot

# GENERAL SETTINGS
str(temp_gistemp)
data=temp_gistemp
fx=sample_mk_function
method="all"
nperm=10
alpha_local=0.05
alpha_global=0.05
null_distribution="normal" # defines if threshold based on alpha level is drawn from normal or t distribution
seed=NULL
block_size=NULL
verbose=TRUE

# FOR TESTING LOOP IN perm_dist.R
i=1

# get results
# nperm = 100

# detrend the data and try! -> check false positive rate
sen0 <- function(y,x){
  zyp.slopediff <- function(i, xx, yy, n) (yy[1:(n - i)] - yy[(i + 1):n])/(xx[1:(n - i)] - xx[(i + 1):n])
  n <- length(y)
  if (missing(x)) x <- c(1:n)
  slopes <- unlist(lapply(1:(n - 1), zyp.slopediff, x, y, n))
  return(median(slopes[is.finite(slopes)], na.rm=TRUE))
}


data_detrend<- data %>% apply(1:2, # apply(1:2,...) will apply function to every cell
                              function(x)
    {
      (x- 1:length(x)*sen0(x))
    }
  )

data_detrend<-  aperm(data_detrend, c(2,3,1)) # transpose it to put lat & long in the first dimensions again

# perm_dist performs obtains the permutation distribution for maxT, STCS and others specified in the respective functions
perm_results<- perm_dist(data=data_detrend, fx=fx, nperm=nperm, alpha_local=alpha_local,
                         alpha_global=alpha_global, null_distribution=null_distribution,
                         seed=seed, block_size=block_size, verbose=verbose)

# test results for nperm=1000, detrended data
# saveRDS(perm_results, file = "testing/detrended_temp_data_nperm_1000.rds")
perm_results<- readRDS("testing/detrended_temp_data_nperm_1000.rds")
str(perm_results)
perm_results$stcs_maxT[!is.finite(perm_results$stcs_maxT)] <- 0

# bootstrap check to check false positive rate
sim_length<- 1000
bootstrap_sample<- 100
fpr_length<- 100
alpha<- 0.05

fpr_sim_stcs<- vector(length = sim_length)
fpr_sim_maxt<- vector(length = sim_length)
fpr_sim_bivariate<- vector(length = sim_length)
for(j in 1:sim_length){

  fpr_stcs<- vector(length = fpr_length)
  fpr_maxt<- vector(length = fpr_length)
  fpr_bivariate<- vector(length = fpr_length)

  for (i in 1:fpr_length){
    ind<- sample(x = length(perm_results$maxT), size = bootstrap_sample, replace = TRUE)
    tmp_stcs<- perm_results$stcs[ind]
    tmp_maxt<- perm_results$maxT[ind]
    tmp_stcs_maxt<- perm_results$stcs_maxT[ind]
    tmp_stcs_maxt[!is.finite(tmp_stcs_maxt)]<- 0

    # gets the threshold for the current sample
    q_thr_stcs<- quantile(tmp_stcs, probs = 1-alpha, names = FALSE)
    q_thr_maxt<- quantile(tmp_maxt, probs = 1-alpha, names = FALSE)
    q_thr_stcs_maxt<- quantile(tmp_stcs_maxt, probs = 1-alpha, names = FALSE)

    # retrieve false positives, i.e. values above current sample based threshold
    fpr_stcs[i]<- tmp_stcs[length(tmp_stcs)] > q_thr_stcs
    fpr_maxt[i]<- tmp_maxt[length(tmp_maxt)] > q_thr_maxt

    # naive bivariate
    # (either STCS or maxT is significant but on 0.025 alpha each)
    fpr_bivariate[i]<- tmp_stcs[length(tmp_stcs)] > quantile(tmp_stcs, probs = 1-alpha/2, names = FALSE) | tmp_stcs_maxt[length(tmp_stcs_maxt)] > quantile(tmp_stcs_maxt, probs = 1-alpha/2, names = FALSE)
    # slightly too liberal - mean and median =.055


    # using bivariate empirical cdf
    # library(bivariate)
    # biv_distr<- ebvcdf(tmp_stcs, tmp_stcs_maxt)
    # fpr_bivariate[i]<- biv_distr(tmp_stcs[length(tmp_stcs)], tmp_stcs_maxt[length(tmp_stcs_maxt)]) > (1-alpha)
    # Too conservative

    # Using empirical quantile contour line
    # contour_line<- get_contour(tmp_stcs, tmp_stcs_maxt, alpha = alpha)
    # fpr_bivariate[i]<- (tmp_stcs[length(tmp_stcs)] > max(quantile(tmp_stcs, probs = 1-alpha/2, names = FALSE), max(contour_line$x))) | (tmp_stcs_maxt[length(tmp_stcs_maxt)] > max(quantile(tmp_stcs_maxt, probs = 1-alpha/2, names = FALSE), max(contour_line$y)))


  }
  fpr_sim_stcs[j]<- sum(fpr_stcs)/length(fpr_stcs)
  fpr_sim_maxt[j]<- sum(fpr_maxt)/length(fpr_maxt)
  fpr_sim_bivariate[j]<- sum(fpr_bivariate)/length(fpr_bivariate)

}
par(mfrow = c(1, 3))
hist(fpr_sim_stcs)
hist(fpr_sim_maxt)
hist(fpr_sim_bivariate)
par(mfrow = c(1, 1))

summary(fpr_sim_stcs)
summary(fpr_sim_maxt)
summary(fpr_sim_bivariate)

# ACCURATE
# bootstrap sample size of 100 gives same results as 300

# following part is based on some previous trials and currently does not work
#*******************
# visualize results for multivariate threshold
# NOTE: outdated df creation
results_df<- do.call(rbind, perm_results$stcs_mvt[[1]]) %>% as.data.frame
results_df$results_length<- lapply(perm_results$stcs_mvt[[1]], function(x) x$results %>% length) %>% do.call(c, .)
results_df$results<- NULL
results_df$id<- 1

nperm = 1000
for(i in 2:nperm){
  tmp<- do.call(rbind, perm_results$stcs_mvt[[i]]) %>% as.data.frame
  tmp$results_length<- lapply(perm_results$stcs_mvt[[i]], function(x) x$results %>% length) %>% do.call(c, .)
  tmp$results<- NULL
  tmp$id<- i
  results_df<- rbind(results_df, tmp)
}
results_df<- unnest(results_df)

library(ggalt)

# all cluster results
results_encircle<- results_df %>% filter(id==3)
results_df %>%
  ggplot(aes(x = results_length, y = abs(maxT), color = (id-25))) +
  geom_point() +
  #xlim(c(-1000,12000)) +
  #ylim(c(-1, 8)) +
  scale_fill_distiller(palette='Spectral') +
  geom_point(aes(x = results_length, y = abs(maxT)),
                data = results_encircle,
                color = "red",
                size = 2)+
  geom_encircle(aes(x = results_length, y = abs(maxT)),
             data = results_encircle,
             color = "red",
             size = 2,
             expand = .01) +
  theme_bw()

# filtered by maximum
results_df %>%
  group_by(id) %>%
  filter(results_length==max(results_length)) %>%
  ggplot(aes(x = results_length, y = abs(maxT), color = (id-25))) +
  geom_point() +
  # xlim(c(0,1000)) +
  scale_fill_distiller(palette='Spectral') +
  geom_point(aes(x = results_length, y = abs(maxT)),
             data = results_encircle %>% filter(results_length==max(results_length)) ,
             color = "red",
             size = 2)+
  geom_encircle(aes(x = results_length, y = abs(maxT)),
                data = results_encircle %>% filter(results_length==max(results_length)),
                color = "red",
                size = 2,
                expand = .08) +
  theme_bw()

# threshold results
out<- threshold_data(perm_results=perm_results, alpha_local=alpha_local,
                     alpha_global=alpha_global, data_dim=dim(data),
                     null_distribution=null_distribution)



####### TEST WITH LAI DATA

load("/home/jose/LAI/data/CHEN_RANGA/AVHRR/yearly_mean/lai_data.RData")

# thresholds for perm distributions

#need to reformat data
tibble_list_to_3d_array<- function(data){
  library(raster)
  out<- matrix(NA, ncol = data$lat %>% unique %>% length, nrow = data$lon %>% unique %>% length)
  data_long<- data %>% dplyr::select(data) %>% ungroup %>% unnest(cols = data)
  n_years<- data_long$t %>% unique %>% length
  data_list<- vector(mode = "list", length = unique(data_long$t) %>% length)
  i<- 1
  for(year in unique(data_long$t)){
    data_raster<- data_long %>% filter(t==year) %>% rename(x = lon, y = lat) %>%
      dplyr::select(x, y, lai) %>% rasterFromXYZ(crs = CRS("+init=epsg:4326"))
    data_list[[i]]<- data_raster$lai %>% as.matrix
    i<- i+1
  }
  data_array<- array(NA, dim= c(dim(data_list[[1]]), length(data_list)))
  for(i in 1:n_years) data_array[,,i]<- data_list[[i]]
  return(data_array)
}

data_lai<- tibble_list_to_3d_array(data)
rm(data)

data_lai_res<- perm_dist(data=data_lai, fx=fx, nperm=100, alpha_local=alpha_local,
          alpha_global=alpha_global, null_distribution=null_distribution,
          seed=seed, block_size=block_size, verbose=verbose)
rm(data)
saveRDS(data_lai_res, file = "testing/lai_res_nperm_1000.rds")



results_df<- do.call(rbind, data_lai_res$stcs_mvt[[1]]) %>% as.data.frame
results_df$results_length<- lapply(data_lai_res$stcs_mvt[[1]], function(x) x$results %>% length) %>% do.call(c, .)
results_df$results<- NULL
results_df$id<- 1

for(i in 2:nperm){
  tmp<- do.call(rbind, data_lai_res$stcs_mvt[[i]]) %>% as.data.frame
  tmp$results_length<- lapply(data_lai_res$stcs_mvt[[i]], function(x) x$results %>% length) %>% do.call(c, .)
  tmp$results<- NULL
  tmp$id<- i
  results_df<- rbind(results_df, tmp)
}
results_df<- unnest(results_df)

library(ggalt)

# all cluster results
results_encircle<- results_df %>% filter(id==nperm)
results_encircle %>% filter(abs(maxT) >= quantile(abs(data_lai_res$maxT), probs = .975) | results_length >= quantile(results_length, probs=.975))

results_df %>%
  ggplot(aes(x = results_length, y = abs(maxT), color = (id-25))) +
  geom_point() +
  #xlim(c(-1000,12000)) +
  #ylim(c(-1, 8)) +
  scale_fill_distiller(palette='Spectral') +
  geom_point(aes(x = results_length, y = abs(maxT)),
             data = results_encircle,
             color = "red",
             size = 2)+
  geom_encircle(aes(x = results_length, y = abs(maxT)),
                data = results_encircle,
                color = "red",
                size = 2,
                expand = .01) +
  theme_bw()

    # filtered by maximum
results_df %>%
  group_by(id) %>%
  filter(results_length==max(results_length)) %>%
  ggplot(aes(x = results_length, y = abs(maxT), color = (id-25))) +
  geom_point() +
  # xlim(c(0,1000)) +
  scale_fill_distiller(palette='Spectral') +
  geom_point(aes(x = results_length, y = abs(maxT)),
             data = results_encircle %>% filter(results_length==max(results_length)) ,
             color = "red",
             size = 2)+
  geom_point(aes(x = results_length, y = abs(maxT)),
                data = results_encircle, #%>% filter(results_length==max(results_length)),
                color = "red",
                size = 2) +
  theme_bw() +
  geom_hline(yintercept = quantile(abs(data_lai_res$maxT), probs = .975)) +
  geom_vline(xintercept = quantile(abs(data_lai_res$stcs), probs = .975)) +
  geom_vline(xintercept = quantile(abs(data_lai_res$stcs), probs = .95), col = "green") +
  geom_vline(xintercept = quantile(abs(data_lai_res$stcs), probs = .9), col = "blue") +
  geom_hline(yintercept = quantile(abs(data_lai_res$maxT), probs = .95), col = "green") +
  geom_hline(yintercept = quantile(abs(data_lai_res$maxT), probs = .9), col = "blue")


