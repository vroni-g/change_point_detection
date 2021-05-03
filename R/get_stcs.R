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


get_stcs<- function(data, alpha_local, null_distribution, data_dim, data_info=NULL){
  if(null_distribution == "normal") thr<- qnorm(1-alpha_local/2)
  if(null_distribution == "t") thr<- qt(1-alpha_local/2, df = data_dim[3]-2)
  if(null_distribution == "brownian_motion"){ 
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

  if(!is.null(data_info)){
    tmp<- matrix(NA, ncol = data_info$ncol, nrow = data_info$nrow)
    tmp[data_info$wh.sel]<- data
    data<- tmp
  }


  pixel_sign<- sign(data)
  pixel_significant<- abs(data)>thr
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
  if(clusters_neg$cluster.count>0){ # only join negative cluster if they exist to avoid cluster of size 0 and s-Inf values
    clusters_sep[[1]]<- clusters_pos$clusters + clusters_neg$clusters
    clusters_sep[[2]]<- c(clusters_pos$cluster.count, clusters_neg$cluster.count)
  } else{
    clusters_sep[[1]]<- clusters_pos$clusters
    clusters_sep[[2]]<- clusters_pos$cluster.count
  }
  names(clusters_sep)<- c("clusters", "cluster.count")
  # clusters_sep is a matrix same size as data containing cluster ids (ints > 0) for cells of negative and positive clusters

  stcs<- max(clusters_sep$cluster.count, na.rm = TRUE)

  allcluster_max <- c()
  clusters_sep$cluster.max <- vector(length = length(clusters_pos$cluster.count))
  for (i in 1:length(clusters_sep$cluster.count)){ # retrieve maximum of each cluster
    clust_max <- data[clusters_sep$clusters==i] %>% abs %>%
      max(.,na.rm = TRUE)
    clusters_sep$cluster.max[i] <- clust_max # assign each cluster its maximum
    allcluster_max <- c(allcluster_max, clust_max)
  }
  stcs_maxT_all <- max(allcluster_max, na.rm = TRUE) # get maximum of all cluster maxima regardless the cluster size
  if (!is.finite(stcs_maxT_all)) stcs_maxT_all <- 0

  return(list(stcs=stcs, clusters=clusters_sep, stcs_maxT_all=stcs_maxT_all, original_stat = data))#, stcs_maxT = stcs_maxT))
}
