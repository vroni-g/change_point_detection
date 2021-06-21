#' Supra-Threshold Cluster Size
#'
#' From data in a matrix, finds the largest cluster size of grid cells above
#' the specified threshold.
#' @param data Matrix containing data. Should be a statistical image, i.e., the
#' test statistic of each grid cell
#' @param alpa_local Alpha level for derivation of supra-threshold cluster
#' @param data_dim dimensions of original data. Used for calculating df for
#' t distributed test statistics, ignored if the test statistic is normal
#' @return A number, stcs, specifying the size, in number of grid cells, of the largest
#' exceedance cluster and a matrix containing the detected clusters.
#' @export get_stcs
#' @import osc


get_stcs<- function(data, alpha_local, null_distribution, data_dim, tippet=TRUE){
  if(null_distribution == "normal") thr<- qnorm(1-alpha_local/2)
  if(null_distribution == "t") thr<- qt(1-alpha_local/2, df = data_dim[3]-2)
  
  if(null_distribution == "brownian_motion"){  # for Recursive CUSUM
  # hard coded thresholds derived (experimentally) from implementation in 
  # strucchange package to obtain p values for empirical processes with limiting 
  # process brownian motion
  # https://github.com/cran/strucchange/blob/master/R/efp.R line 477
  # for univariate data only
    if(alpha_local == 0.10) thr <- 0.849923734
    if(alpha_local == 0.05) thr <- 0.947898221
    if(alpha_local == 0.025) thr <- 1.03651316
    if(alpha_local == 0.01) thr <- 1.142973511
  }
  # if(null_distribution == "brownian_bridge_increments"){ # for OLS-MOSUM
  #   # if(alpha_local == 0.10) thr <- 1.0369799425 # for h = 0.12
  #   # if(alpha_local == 0.05) thr <- 1.1113399926
  #   # if(alpha_local == 0.025) thr <- 1.1809399861
  #   # if(alpha_local == 0.01) thr <- 1.2639599731
  #   if(alpha_local == 0.10) thr <- 1.338599926 # for h = 0.325
  #   if(alpha_local == 0.05) thr <- 1.4618499877
  }
  if(null_distribution == "p-values"){
    thr <- alpha_local
  }
  
  pixel_sign<- sign(data)
  
  if(null_distribution == "p-values"){
    pixel_significant<- abs(data)<thr
  } else {
    pixel_significant<- abs(data)>thr
  }
  
  pixel_result<- pixel_sign*pixel_significant

  # positive
  pixel_result_pos<- pixel_result
  pixel_result_pos[is.na(pixel_result_pos)]<- -999
  clusters_pos<- osc::cca(pixel_result_pos,count.cells = TRUE, s=1, mode = 2, # only values >0 are included in osc:cca
                          count.max  = length(pixel_sign)) # this returns a matrix of size data that contains integers > 0 that are cluster ids
  stcs_pos<- max(clusters_pos$cluster.count)

  nclust_pos<- length(clusters_pos$cluster.count)
  clusters_sep<- vector(mode = "list", length = 2)

  # negative
  pixel_result_neg<- pixel_result
  pixel_result_neg[pixel_result_neg == -1] = 10 # swap signs so that originally negative values will now be considered in osc:cca
  pixel_result_neg[pixel_result_neg == 1] = -10
  pixel_result_neg[is.na(pixel_result_neg)] = -999
  clusters_neg<- osc::cca(pixel_result_neg,count.cells = TRUE, s=1, mode = 2,
                          count.max = length(pixel_sign))
  stcs_neg<- max(clusters_neg$cluster.count)

  # join positive and negative clusters 
  clusters_neg$clusters[clusters_neg$clusters > 0]<- clusters_neg$clusters[clusters_neg$clusters > 0] + nclust_pos
  if(clusters_neg$cluster.count[1]>0){ # only join negative cluster if they exist to avoid cluster of size 0 and s-Inf values
    clusters_sep[[1]]<- clusters_pos$clusters + clusters_neg$clusters
    clusters_sep[[2]]<- c(clusters_pos$cluster.count, clusters_neg$cluster.count)
  } else{
    clusters_sep[[1]]<- clusters_pos$clusters
    clusters_sep[[2]]<- clusters_pos$cluster.count
  }
  names(clusters_sep)<- c("clusters", "cluster.count")
  # clusters_sep is a matrix same size as data containing cluster ids (ints > 0) for cells of negative and positive clusters

  stcs<- max(clusters_sep$cluster.count, na.rm = TRUE)
  stcs_idx <- which(clusters_sep$cluster.count==stcs)[1]
  if(null_distribution == "p-values"){
    stcs_minP <- min(data[clusters_sep$clusters==stcs_idx], na.rm = TRUE)
  } else {
    stcs_maxT <- max(data[clusters_sep$clusters==stcs_idx], na.rm = TRUE)
  }
  
  #********************************* TIPPET ************************************
  
  if(tippet){
    library(magrittr)
    
    if(null_distribution == "p-values"){
      get_pi <- function(i, clusters_sep){
        clust_min <- data[clusters_sep$clusters==i] %>% abs %>%
          min(.,na.rm = TRUE)
        return(clust_min)
      }
      pis <- lapply(1:length(clusters_sep$cluster.count), get_pi, clusters_sep = clusters_sep) %>% 
        unlist
      clusters_sep$cluster.min <- pis
      peak_intensity <- min(clusters_sep$cluster.min, na.rm = TRUE) # get minimum of all cluster minima regardless the cluster size
      if (!is.finite(peak_intensity)) peak_intensity <- NA
      
      return(list(stcs=stcs, clusters=clusters_sep, peak_intensity=peak_intensity, 
                  stcs_minP=stcs_minP, original_stat = data))
      
    } else {
      get_pi <- function(i, clusters_sep){
        clust_max <- data[clusters_sep$clusters==i] %>% abs %>%
          max(.,na.rm = TRUE)
        return(clust_max)
      }
      pis <- lapply(1:length(clusters_sep$cluster.count), get_pi, clusters_sep = clusters_sep) %>% 
        unlist
      clusters_sep$cluster.max <- pis
      peak_intensity <- max(clusters_sep$cluster.max, na.rm = TRUE) # get maximum of all cluster maxima regardless the cluster size
      if (!is.finite(peak_intensity)) peak_intensity <- NA
      
      return(list(stcs=stcs, clusters=clusters_sep, peak_intensity=peak_intensity,
                  stcs_maxT=stcs_maxT, original_stat = data))
    }
    
    #******************************* TIPPET END ********************************

    
  } else { # no tippet (saves usually quite some time as cluster mins/max don't have to be derived)
    
    if(null_distribution == "p-values"){
      return(list(stcs=stcs, stcs_minP = stcs_minP, clusters=clusters_sep, original_stat = data))
    } else {
      return(list(stcs=stcs, stcs_maxT = stcs_maxT, clusters=clusters_sep, original_stat = data))
    }
  }
  
}
