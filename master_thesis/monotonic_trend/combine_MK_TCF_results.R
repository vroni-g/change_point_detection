#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# !! Input Data is not available in GitHub Repository !!


# Script to combine the single results of Mann Kendall permutations
# and to retrieve Tippet combining function (TCF) values
#***************************************************************************
library(tidyverse)

dir <- "/home/veronika/CPD/results/nperm_BU_MK_loc10" # HPC directory with result files
files <-list.files(dir, pattern = "\\.rds$", full.names = T)
files
orig_ind <- 4 # check which file is the last one (nperm) and set it

# retrieve maximum statistics of all permutations
q_results <- lapply(files, function(x) readRDS(x)[[1]]) %>%
  do.call(rbind, .) %>%
  as.data.frame()

# get empirical distribution of peak intensity and stcs
dis_peak_intensity <- ecdf(q_results[,"peak_intensity"])
dis_stcs<- ecdf(q_results[,"stcs"])

get_wt <- function(path, dis_peak_intensity, dis_stcs, nperm, last = FALSE){
   clust_perm <- readRDS(path)[[2]]

    # retrieve p-values for cluster size and cluster maximum for each cluster in the current permutation

    get_p <- function(j, clust_perm,  dis_stcs, dis_peak_intensity){
      p_peak_intensity <- 1 - dis_peak_intensity(clust_perm$cluster.max[j]) + 1/nperm
      if(p_peak_intensity<=0) p_peak_intensity <- 0.000001

      p_stcs <- 1 - dis_stcs(clust_perm$cluster.count[j]) + 1/nperm
      if(p_stcs<=0) p_stcs <- 0.000001

      w <- 1 - min(p_peak_intensity, log(p_stcs))
      if(!is.finite(w)) return(NA)

      return(w)
    }

    js <- seq(1:length(clust_perm$cluster.count))
    # map over each cluster
    w_tmp <- purrr::map(js, get_p, clust_perm = clust_perm,
                       dis_peak_intensity = dis_peak_intensity,
                        dis_stcs = dis_stcs) %>% unlist

    # return maximum tippet value for current permutation
    if(last){
      return(list(max(w_tmp, na.rm = TRUE), w_tmp))
    } else {
      return(max(w_tmp, na.rm = TRUE))
    }
}

# map over all permutations but the last one and retrieve each maximum tippet statistic
or <- files[orig_ind]
files <- files[-orig_ind]

wt <- purrr::map(files, get_wt, dis_peak_intensity = dis_peak_intensity,
                 dis_stcs = dis_stcs, nperm = length(files)+1) %>% unlist

# get values for last permutation and also return original tippet values
l <- get_wt(or, dis_peak_intensity = dis_peak_intensity,
            dis_stcs = dis_stcs, last = TRUE, nperm = length(files)+1)

# append last maximum tippet value to the rest
wt <- c(wt, unlist(l[[1]]))

original_stat <- readRDS(or) %>% .[[3]]
original_cluster <- readRDS(or) %>% .[[2]]

results <- list(maxT = q_results[,"maxT"], stcs = q_results[,"stcs"],
                 peak_intensity = q_results[,"peak_intensity"], wt = wt, original_wt = l[[2]],
                 original_cluster = original_cluster, original_stat = original_stat)

# relative path within GitHub Repository
saveRDS(results, "master_thesis/results/BU_LAI_MK_nperm_1000_al10.rds")
