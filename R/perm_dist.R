#' Permutation distribution
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

perm_dist<- function(data, fx, nperm=1000,
                     alpha_local, alpha_global, null_distribution,
                     block_size = NULL, seed, verbose = TRUE){
  perm_matrix<- perm_matrix(nobs = dim(data)[3], nperm = nperm, block_size = block_size, seed = seed)
  #**************testing only**********************
  # remove last/original one to see if it works for a permutation which has clusters by chance
  nperm <- nperm - 1
  perm_matrix <- perm_matrix[1:nperm, 1:dim(perm_matrix)[2]]
  #**************testing only**********************
  perm_matrix

  maxT<- vector(length = nperm)
  stcs<- vector(length = nperm)
  stcs_maxT<- vector(length = nperm)
  stcs_maxT_all <- vector(length = nperm)
  perm_results <- vector(length = nperm, mode = 'list') # save all cluster results to derive p-values for cluster
  cat("starting permutations:\n")

  for(i in 1:nperm){
    cat("Permutation",i,"\n")
    tmp <- apply(data[,,perm_matrix[i,]], 1:2, fx)
    maxT[i]<- max(abs(as.vector(tmp)), na.rm = TRUE)
    tmp_stcs <- get_stcs(tmp, alpha_local, null_distribution)
    stcs[i]<- tmp_stcs$stcs
    perm_results[i] <- list(tmp_stcs$clusters)
    stcs_maxT[i]<- tmp_stcs$stcs_maxT
    stcs_maxT_all[i] <- tmp_stcs$stcs_maxT_all
    if(verbose) if((i%%10)==0) cat(i,"\n")
  }

  # get empirical distribution of maxT_all and stcs
  dis_maxT_all <- ecdf(stcs_maxT_all)
  dis_stcs<- ecdf(stcs)
  wt <- vector(length = nperm)

  for(i in 1:nperm){
    w_tmp <- vector(length = length(perm_results[[i]]$cluster.count))
    clust_perm <- perm_results[[i]]
    for(j in 1:length(clust_perm$cluster.count)){
      # retrieve p-values for each cluster
      p_maxT_all <- 1 - dis_maxT_all(clust_perm$cluster.max[j])
      p_stcs <- 1 - dis_stcs(clust_perm$cluster.count[j])
      # combine in new test statistic
      w <- 1 - min(log(p_maxT_all), log(p_stcs))
      if (is.finite(w)){
        w_tmp <- c(w_tmp, w)
      } else{
        w_tmp <- c(w_tmp, 0)
      }
    }
    # get max test statistic value for each permutation
    wt[i] <- max(w_tmp, na.rm = TRUE)
  }


  cat("finished!\n\n")
  return(list(maxT = maxT, stcs = stcs, wt = wt, stcs_maxT = stcs_maxT,
              original_results = tmp, original_wt = w_tmp))
}
