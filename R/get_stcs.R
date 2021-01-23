#' Supra-Threshold Cluster Size
#'
#' From data in a matrix, finds the largest cluster size of grid cells above
#' the specified threshold.
#' @param data Matrix containing data. Should be a statistical image, i.e., the
#' test statistic of each grid cell
#' @param thr Threshold of significance.
#' @return A number, stcs, specifying the size, in number of grid cells, of the largest
#' exceedance cluster and a matrix containing the detected clusters.
#' @param data_dim dimesions of original data. Used for calculating df for
#' t distributed test statistics, ignored if the test statistic is normal
#' @export get_stcs
#' @import osc

get_stcs<- function(data, alpha_local, null_distribution, data_dim){
  if(null_distribution == "normal") thr<- qnorm(1-alpha_local/2)
  if(null_distribution == "t") thr<- qt(1-alpha_local/2, df = data_dim[3]-2)

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

  # join
  clusters_neg$clusters[clusters_neg$clusters > 0]<- clusters_neg$clusters[clusters_neg$clusters > 0] + nclust_pos
  clusters_sep[[1]]<- clusters_pos$clusters + clusters_neg$clusters
  clusters_sep[[2]]<- c(clusters_pos$cluster.count, clusters_neg$cluster.count)
  names(clusters_sep)<- c("clusters", "cluster.count")
  # clusters_sep is a matrix same size as data containing cluster ids (ints > 0) for cells of negative and positive clusters

  stcs<- max(clusters_sep$cluster.count, na.rm = TRUE)
  #stcs_idx<- which(length(clusters_sep$cluster.count)==stcs)
  #!!!!!!!!!!!!!!
  # shouldn't that be
  stcs_idx<- which(clusters_sep$cluster.count==stcs)
  #!!!!!!!!!!!!!!
  stcs_cluster_results<- data[clusters_sep$clusters==stcs_idx] # retrieve all cells (by position in matrix?) that belong to the biggest cluster
  stcs_maxT<- max(stcs_cluster_results, na.rm = TRUE)

  # within cluster properties --- maxT works fine, all others are similar or worse
  # stcs_mvt<- vector(length = length(clusters_sep$cluster.count), mode = "list")
  # for (i in 1:length(clusters_sep$cluster.count)){
  #   #stcs_mvt[[i]]<- vector(length = 11, mode = "list")
  #   #print(i)
  #
  #   #get results for cluster i and save for later
  #   cluster_results<- data[clusters_sep$clusters==i]
  #   stcs_mvt[[i]]$results<- cluster_results
  #
  #   # maxT
  #   stcs_mvt[[i]]$maxT<- max(cluster_results, na.rm = TRUE)
  #
  #   # avgT
  #   stcs_mvt[[i]]$meanT<- mean(cluster_results, na.rm = TRUE)
  #
  #   # medianT
  #   stcs_mvt[[i]]$medianT<- median(cluster_results, na.rm = TRUE)
  #
  #   # quantiles: 0.90, 0.95
  #   stcs_mvt[[i]]$q90T<- unname(quantile(cluster_results, probs = 0.90, na.rm = TRUE))
  #   stcs_mvt[[i]]$q95T<- unname(quantile(cluster_results, probs = 0.95, na.rm = TRUE))
  #
  #   # average of top: 3, 5, 10 grid cells
  #   stcs_mvt[[i]]$meanTop3<- mean(head(sort(cluster_results, decreasing = TRUE), n=3), na.rm = TRUE)
  #   stcs_mvt[[i]]$meanTop5<- mean(head(sort(cluster_results, decreasing = TRUE), n=5), na.rm = TRUE)
  #   stcs_mvt[[i]]$meanTop10<- mean(head(sort(cluster_results, decreasing = TRUE), n=10), na.rm = TRUE)
  #
  #   # average of top: 5% 10%
  #   stcs_mvt[[i]]$meanTop5percent<- mean(head(sort(cluster_results, decreasing = TRUE), n=length(cluster_results)*.05), na.rm = TRUE)
  #   stcs_mvt[[i]]$meanTop10percent<- mean(head(sort(cluster_results, decreasing = TRUE), n=length(cluster_results)*.10), na.rm = TRUE)
  # }

  return(list(stcs=stcs, clusters=clusters_sep, stcs_maxT=stcs_maxT))
}
