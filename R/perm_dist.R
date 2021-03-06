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

  maxT<- vector(length = nperm)
  stcs<- vector(length = nperm)
  stcs_maxT <- vector(length = nperm)
  peak_intensity <- vector(length = nperm)
  perm_results <- vector(length = nperm, mode = 'list') # save all cluster results to derive p-values for cluster
  cat("starting permutations:\n")
  
  # convert data to 2d matrix
  data <- array_to_matrix(data)

  for(i in 1:nperm){
    cat("Permutation",i,"\n")
    tmp1 <- apply(data$Y[perm_matrix[i,],], 2, fx)
    maxT[i]<- max(abs(as.vector(tmp1)), na.rm = TRUE)
    
    # reinsert NA values
    data_info <- data[2:5]
    tmp <- matrix(NA, ncol = data_info$ncol, nrow = data_info$nrow)
    tmp[data_info$wh.sel]<- tmp1
    
    tmp_stcs <- get_stcs(tmp, alpha_local, null_distribution)
    perm_results[i] <- list(tmp_stcs$clusters)
    stcs[i]<- tmp_stcs$stcs
    stcs_maxT[i]<- tmp_stcs$stcs_maxT
    peak_intensity[i] <- tmp_stcs$peak_intensity
    rm(tmp_stcs)
    if(verbose) if((i%%10)==0) cat(i,"\n")
  }

  # get empirical distribution of maxT_all and stcs
  dis_peak_intensity <- ecdf(peak_intensity)
  dis_stcs<- ecdf(stcs)
  dis_stcs_maxT<- ecdf(stcs_maxT)

  get_wt <- function(clust_perm, dis_peak_intensity, dis_stcs, dis_stcs_maxT, nperm, last = FALSE){
    # retrieve p-values for cluster size and cluster maximum for each cluster in the current permutation
    get_p <- function(j, clust_perm, dis_peak_intensity, dis_stcs, nperm){
      p_peak_intensity <- 1 - dis_peak_intensity(clust_perm$cluster.max[j]) + 1/nperm
      if(p_peak_intensity<=0) p_peak_intensity <- 0.000001
      p_stcs <- 1 - dis_stcs(clust_perm$cluster.count[j]) + 1/nperm
      if(p_stcs<=0) p_stcs <- 0.000001
      
      p_stcs_maxT <- 1 - dis_stcs_maxT(clust_perm$cluster.max[j]) + 1/nperm
      if(p_stcs_maxT<=0) p_stcs_maxT <- 0.000001
      
      # combine in new test statistic
      w <- 1 - min(log(p_peak_intensity), log(p_stcs), log(p_stcs_maxT))
      if (is.finite(w)){
        return(w)
      } else{
        return(0)
      }
    }
    js <- seq(1:length(clust_perm$cluster.count))
    # map over each cluster
    w_tmp <- purrr::map(js, get_p, clust_perm = clust_perm, dis_peak_intensity = dis_peak_intensity, 
                        dis_stcs = dis_stcs, dis_stcs_maxT = dis_stcs_maxT, nperm = nperm)
    # return maximum tippet value for current permutation
    if(last){
      return(list(max(unlist(w_tmp), na.rm = TRUE), unlist(w_tmp)))
    } else {
      return(max(unlist(w_tmp), na.rm = TRUE))
    }
  }
  # map over all permutations but the last one and retrieve each maximum tippet statistic
  wt <- purrr::map(perm_results[1:(nperm-1)], get_wt,
                   dis_peak_intensity = dis_peak_intensity, dis_stcs = dis_stcs, nperm = nperm) %>%
    unlist()
  # get values for last permutation and also return original tippet values
  l <- get_wt(perm_results[[nperm]], dis_peak_intensity = dis_peak_intensity,
              dis_stcs = dis_stcs, last = TRUE, nperm = nperm)
  # append last maximum tippet value to the rest
  wt <- c(wt, l[[1]])

  cat("finished!\n\n")
  return(list(maxT = maxT, stcs = stcs, wt = wt, peak_intensity = peak_intensity,
              original_results = tmp, original_wt = l[[2]]), original_cluster = perm_results[[nperm]])
}
