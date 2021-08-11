#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# Author of Code: José Cortés
# Modifications: Veronika Grupp


#' Permutation distribution for SLURM cluster
#'
#' Obtains the permutation distribution of the maximum statistic, the STCS and 
#' the peak intensity, named maxT, stcs and peak_intensity, respectively. 
#' Returns a list containing maxT, stcs and peak_intensity:
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
#' @return returns the distribution of maxT, stcs, and peak intensity, each a vector in a list
#' @export perm_dist
#'

perm_dist_SLURM_2d<- function(data, fx, nperm=1000,
                              alpha_local, alpha_global, null_distribution,
                              block_size = NULL, seed, verbose = TRUE, res_dir){

  # convert data to 2d matrix
  data <- array_to_matrix(data)

  perm_mat<- perm_matrix(nobs = nrow(data$Y), nperm = nperm, block_size = block_size, seed = seed)
  cat("starting permutations:\n")

  tmp_fn<- function(i, perm_mat, fx, data){

    library(magrittr)
    #devtools::load_all("/home/veronika/CPD/change_point_detection/")
    devtools::load_all()
    data_info <- data[2:5]
    tmp <- matrix(NA, ncol = data_info$ncol, nrow = data_info$nrow)

    cat("Starting Test for permutation ", i, " at ", date(), "\n")
    tmp1<- apply(data$Y[perm_mat[i,],], 2, fx)
    cat("Test finished for permutation ", i, " at ", date(), "\n")

    if(null_distribution == "p-values"){
      maxT <- min(abs(tmp1), na.rm = TRUE)
    } else{
      maxT<- max(abs(tmp1), na.rm = TRUE)
    }

    # reinsert NA values
    tmp[data_info$wh.sel]<- tmp1


    cat("Starting cluster derivation for permutation ", i, " at ", date(), "\n")
    tmp_stcs<- get_stcs(tmp, alpha_local, null_distribution, tippet = TRUE)
    cat("Clusters derived completely for permutation ", i, " at ", date(), "\n")

    stcs<- tmp_stcs$stcs
    peak_intensity <- tmp_stcs$peak_intensity

    if(i==dim(perm_mat)[[1]]){
      r <- list(c(maxT = maxT, stcs = stcs, peak_intensity=peak_intensity), tmp_stcs$clusters, tmp_stcs$original_stat)
      f <- paste0(res_dir, "single_BU_LAI_MK_nperm_", dim(perm_mat)[[1]],"_",i, ".rds")
      saveRDS(r, file = f)
      cat("File saved for permutation ", i)
      return(r)
    } else {
      r <- list(c(maxT = maxT, stcs = stcs, peak_intensity=peak_intensity), tmp_stcs$clusters)
      f <- paste0(res_dir, "single_BU_LAI_MK_nperm_", dim(perm_mat)[[1]],"_",i, ".rds")
      saveRDS(r, file = f)
      cat("File saved for permutation ", i)
      return(r)
    }
  }

# start parallelization
  library(clustermq)
  results<- Q(tmp_fn,
              i=1:nperm,
              const = list(perm_mat = perm_mat,
                           data = data,
                           fx = fx),
              export = list(alpha_local = alpha_local,
                            null_distribution = null_distribution),
              n_jobs = 4,
              template = list(job_name = "MK_BU",
                              partition = "all",
                              log_file = paste0("/home/veronika/CPD/logs/BU_MK.txt"), # adjust to your path
                              memory = 10000,
                              n_cpus = 1),
              fail_on_error = FALSE,
              verbose = TRUE)
}




