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

perm_dist_SLURM<- function(data, fx, nperm=1000,
                     alpha_local, alpha_global, null_distribution,
                     block_size = NULL, seed, verbose = TRUE){
  perm_matrix<- perm_matrix(nobs = dim(data)[3], nperm = nperm, block_size = block_size, seed = seed)
  maxT<- vector(length = nperm)
  stcs<- vector(length = nperm)
  stcs_maxT<- vector(length = nperm)
  stcs_maxT_all <- vector(length = nperm)
  perm_results <- vector(length = nperm, mode = 'list') # save all cluster results to derive p-values for cluster
  cat("starting permutations:\n")

  tmp_fn<- function(i, perm_matrix, fx, data){
    library(devtools)
    library(magrittr)
    load_all()
    tmp<- apply(data[,,perm_matrix[i,]], 1:2, fx)
    maxT<- max(abs(as.vector(tmp)), na.rm = TRUE)
    tmp_stcs<- get_stcs(tmp, alpha_local, null_distribution)
    clust_results <- list(tmp_stcs$clusters)
    stcs<- tmp_stcs$stcs
    stcs_maxT<- tmp_stcs$stcs_maxT
    stcs_maxT_all <- tmp_stcs$stcs_maxT_all
    #if(verbose) if((i%%10)==0) cat(i,"\n")
    return(list(c(maxT = maxT, stcs = stcs, stcs_maxT = stcs_maxT, stcs_maxT_all=stcs_maxT_all), clust_results))
  }

  library(clustermq)
  results<- Q(tmp_fn,
              i=1:nperm,
              const = list(perm_matrix = perm_matrix,
                           data = data,
                           fx = fx),
              export = list(alpha_local = alpha_local,
                            null_distribution = null_distribution),
              n_jobs = 100,
              template = list(job_name = "Wt_MTPC",
                              partition = "all",
                              log_file = "test_clustadjWt.txt",
                              memory = 10000,
                              n_cpus = 1),
              fail_on_error = FALSE,
              verbose = TRUE)
  #return(results)
  #results <- readRDS("detrended_temp_data_Wtadjust_nperm_10.rds")
  library(magrittr)
  q_results <- lapply(results, function(x) x[[1]]) %>%
    do.call(rbind, .)
  # extract all perm_results and combine in list
  perm_results <- sapply(results, function(x) x[[2]])

  # get empirical distribution of maxT_all and stcs
  dis_maxT_all <- ecdf(q_results[,4])
  dis_stcs<- ecdf(q_results[,2])
  wt <- vector(length = nperm)

  for(i in 1:nperm){
    clust_perm <- perm_results[[i]]
    w_tmp <- vector(length = length(clust_perm$cluster.count))
    #cat("Number of cluster: ", length(clust_perm$cluster.count))
    for(j in 1:length(clust_perm$cluster.count)){
      # retrieve p-values for each cluster
      p_maxT_all <- 1 - dis_maxT_all(clust_perm$cluster.max[j]) + 1/nperm
      if(p_maxT_all<=0) p_maxT_all <- 0.000001
      p_stcs <- 1 - dis_stcs(clust_perm$cluster.count[j]) + 1/nperm
      if(p_stcs<=0) p_stcs <- 0.000001
      # combine in new test statistic
      w <- 1 - min(log(p_maxT_all), log(p_stcs))
      if (is.finite(w)){
        w_tmp[j] <- w
      } else{
        w_tmp[j] <- 0
      }
    }
    # get max test statistic value for each permutation
    wt[i] <- max(w_tmp, na.rm = TRUE)
  }

  cat("finished!\n\n")
  return(list(maxT = q_results[,1], stcs = q_results[,2], stcs_maxT = q_results[,3], wt = wt, original_wt = w_tmp))
}

