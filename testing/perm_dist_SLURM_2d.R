#' Permutation distribution for SLURM cluster
#'
#' Obtains the permutation distribution of the maximum statistic and the STCS,
#' named maxT and stcs, respectively. Returns a list containing maxT and stcs:
#' each a vector of length nperm where each entry is the maximum of the nth
#' permutation. Together form the permutation distribution, of which
#' the (1-alpha)th percentile is the threshold of significance adjusted
#' for multiple testing.
#'
#' @param data Data in a 3d array, where the first two dimensions are the
#' physical dimensions (e.g., lon and lat), and the third one is time.
#' @param fx function to be applied at each grid cell. Should be self
#' sufficient - no extra arguments needed. Should return only the test
#' statistic
#' @param nperm number of permutations. Defaults to 1000
#' @param verbose Counter returning when the function is done with 10 function
#' calls
#' @return returns the distribution of maxT, stcs, and stcs_mvt, each a vector in a list
#' @export perm_dist
#'

perm_dist_SLURM_2d<- function(data, fx, nperm=1000,
                              alpha_local, alpha_global, null_distribution,
                              block_size = NULL, seed, verbose = TRUE){

  # convert data to 2d matrix
  data <- array_to_matrix(data)
  
  perm_mat<- perm_matrix(nobs = nrow(data$Y), nperm = nperm, block_size = block_size, seed = seed)
  cat("starting permutations:\n")

  tmp_fn<- function(i, perm_mat, fx, data){
    library(magrittr)
    devtools::load_all("/home/veronika/CPD/change_point_detection/")
    cat("Starting Test for permutation ", i, " at ", date(), "\n")
    tmp1<- apply(data$Y[perm_mat[i,],], 2, fx)
    cat("Test finished for permutation ", i, " at ", date(), "\n")
    if(null_distribution == "p-values"){
      minP <- min(abs(tmp1), na.rm = TRUE)
    } else{
      maxT<- max(abs(tmp1), na.rm = TRUE)
    }
    
    # reinsert NA values
    data_info <- data[2:5]
    tmp <- matrix(NA, ncol = data_info$ncol, nrow = data_info$nrow)
    tmp[data_info$wh.sel]<- tmp1
    
    cat("Starting cluster derivation for permutation ", i, " at ", date(), "\n")
    tmp_stcs<- get_stcs(tmp, alpha_local, null_distribution, tippet = TRUE)
    cat("Clusters derived completely for permutation ", i, " at ", date(), "\n")
    stcs<- tmp_stcs$stcs
    if(null_distribution == "p-values"){
      stcs_minP <- tmp_stcs$stcs_minP
    }else{
      stcs_maxT <- tmp_stcs$stcs_maxT
    }
    
    peak_intensity <- tmp_stcs$peak_intensity
    
    if(null_distribution == "p-values"){
      
      if(i==dim(perm_mat)[[1]]){
        r <- list(c(minP = minP, stcs = stcs, stcs_minP = stcs_minP, peak_intensity=peak_intensity), tmp_stcs$clusters, tmp_stcs$original_stat)
        f <- paste0("/home/veronika/CPD/results/nperm_BU_MOSUM_loc10_h0.325/single_BU_LAI_MOSUM_nperm_", dim(perm_mat)[[1]],"_",i, ".rds")
        saveRDS(r, file = f)
        cat("File saved for permutation ", i)
        return(r)
      } else {
        r <- list(c(minP = minP, stcs = stcs, stcs_minP = stcs_minP, peak_intensity=peak_intensity), tmp_stcs$clusters)
        f <- paste0("/home/veronika/CPD/results/nperm_BU_MOSUM_loc10_h0.325/single_BU_LAI_MOSUM_nperm_", dim(perm_mat)[[1]], "_",i, ".rds")
        saveRDS(r, file = f)
        cat("File saved for permutation ", i)
        return(r)
      }
      
    }else{
      
      if(i==dim(perm_mat)[[1]]){
        r <- list(c(maxT = maxT, stcs = stcs, stcs_maxT = stcs_maxT, peak_intensity=peak_intensity), tmp_stcs$clusters, tmp_stcs$original_stat)
        f <- paste0("/home/veronika/CPD/results/nperm_BU_MOSUM_loc10/single_BU_LAI_MOSUM_nperm_", dim(perm_mat)[[1]],"_",i, ".rds")
        saveRDS(r, file = f)
        cat("File saved for permutation ", i)
        return(r)
      } else {
        r <- list(c(maxT = maxT, stcs = stcs, stcs_maxT = stcs_maxT, peak_intensity=peak_intensity), tmp_stcs$clusters)
        f <- paste0("/home/veronika/CPD/results/nperm_BU_MOSUM_loc10/single_BU_LAI_MOSUM_nperm_", dim(perm_mat)[[1]], "_",i, ".rds")
        saveRDS(r, file = f)
        cat("File saved for permutation ", i)
        return(r)
      }
    }

  }

  library(clustermq)
  results<- Q(tmp_fn,
              i=1:nperm,
              const = list(perm_mat = perm_mat,
                           data = data,
                           fx = fx),
              export = list(alpha_local = alpha_local,
                            null_distribution = null_distribution),
              n_jobs = 75,
              template = list(job_name = "mosum_BU",
                              partition = "all",
                              log_file = paste0("/home/veronika/CPD/logs/BU_MOSUM_",nperm,"n_75j_10000mem.txt"),
                              memory = 10000,
                              n_cpus = 1),
              fail_on_error = FALSE,
              verbose = TRUE)

  # library(tidyverse)
  # library(magrittr)
  # q_results <- lapply(results, function(x) x[[1]]) %>%
  #   do.call(rbind, .)
  # # extract all perm_results and combine in list
  # perm_results <- sapply(results, function(x) x[[2]])
  # # extract original statistic values
  # original_stat <- results[[nperm]][[3]]
  # rm(results)
  # 
  # res <- list(q_results, perm_results, original_stat)
  # return(res)

  
#   # get empirical distribution of maxT_all and stcs
#   dis_maxT_all <- ecdf(q_results[,3])
#   dis_stcs<- ecdf(q_results[,2])
#
#   get_wt <- function(clust_perm, dis_maxT_all, dis_stcs, nperm, last = FALSE){
#     # retrieve p-values for cluster size and cluster maximum for each cluster in the current permutation
#     get_p <- function(j, clust_perm, dis_maxT_all, dis_stcs, nperm){
#       p_maxT_all <- 1 - dis_maxT_all(clust_perm$cluster.max[j]) + 1/nperm
#       if(p_maxT_all<=0) p_maxT_all <- 0.000001
#       p_stcs <- 1 - dis_stcs(clust_perm$cluster.count[j]) + 1/nperm
#       if(p_stcs<=0) p_stcs <- 0.000001
#       # combine in new test statistic
#       w <- 1 - min(log(p_maxT_all), log(p_stcs))
#       if (is.finite(w)){
#         return(w)
#       } else{
#         return(0)
#       }
#     }
#     js <- seq(1:length(clust_perm$cluster.count))
#     # map over each cluster
#     w_tmp <- purrr::map(js, get_p, clust_perm = clust_perm, dis_maxT_all = dis_maxT_all, dis_stcs = dis_stcs, nperm = nperm)
#     # return maximum tippet value for current permutation
#     if(last){
#       return(list(max(unlist(w_tmp), na.rm = TRUE), unlist(w_tmp)))
#     } else {
#       return(max(unlist(w_tmp), na.rm = TRUE))
#     }
#   }
#   # map over all permutations but the last one and retrieve each maximum tippet statistic
#   wt <- purrr::map(perm_results[1:(nperm-1)], get_wt,
#                    dis_maxT_all = dis_maxT_all, dis_stcs = dis_stcs, nperm = nperm) %>%
#     unlist()
#   # get values for last permutation and also return original tippet values
#   l <- get_wt(perm_results[[nperm]], dis_maxT_all = dis_maxT_all,
#               dis_stcs = dis_stcs, last = TRUE, nperm = nperm)
#   # append last maximum tippet value to the rest
#   wt <- c(wt, l[[1]])
#
#   cat("finished!\n\n")
#   return(list(maxT = q_results[,1], stcs = q_results[,2], stcs_maxT = q_results[,3], wt = wt, original_wt = l[[2]],
#               original_cluster = perm_results[[nperm]], original_stat = original_stat))
}



