# go through the cluster properties / get_stcs step by step:
library("devtools")
load_all()
library(tidyverse)

data=temp_gistemp
# detrend data ----
#***********************
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
data_detrend <-  aperm(data_detrend, c(2,3,1)) # transpose it to put lat & long in the first dimensions again
# set options ----
#***********************
fx=sample_mk_function
method="all"
nperm=4
alpha_local=0.05
alpha_global=0.05
null_distribution="normal" # defines if threshold based on alpha level is drawn from normal or t distribution
seed=NULL
block_size=NULL
verbose=TRUE
data=data_detrend

# perm_dist ----
#***********************
perm_matrix<- perm_matrix(nobs = dim(data)[3], nperm = nperm, block_size = block_size, seed = seed)
nperm <- nperm - 1
perm_matrix <- perm_matrix[1:nperm, 1:dim(perm_matrix)[2]]
maxT<- vector(length = nperm)
stcs<- vector(length = nperm)
stcs_maxT_all<- vector(length = nperm)
#cat("starting permutations:\n")
# for(i in 1:nperm){
#   tmp <- apply(data[,,perm_matrix[i,]], 1:2, fx)
#   perm_results[i] <- tmp
#   maxT[i]<- max(abs(as.vector(tmp)), na.rm = TRUE)
#   tmp_stcs<- get_stcs(tmp, alpha_local, null_distribution)
#   stcs[i]<- tmp_stcs$stcs
#   #stcs_maxT[i]<- tmp_stcs$stcs_maxT
#   stcs_maxT_all <- tmp_stcs$stcs_maxT_all
#   if(verbose) if((i%%10)==0) cat(i,"\n")
# }


# get_stcs step by step ----
#***********************
tmp<- apply(data[,,perm_matrix[3,]], 1:2, fx)
data = tmp
if(null_distribution == "normal") thr<- qnorm(1-alpha_local/2)
if(null_distribution == "t") thr<- qt(1-alpha_local/2, df = data_dim[3]-2)

pixel_sign<- sign(data)
pixel_significant<- abs(data)>thr
pixel_result<- pixel_sign*pixel_significant

# positive
pixel_result_pos<- pixel_result
pixel_result_pos[is.na(pixel_result_pos)]<- -999
clusters_pos<- osc::cca(pixel_result_pos,count.cells = TRUE, s=1, mode = 2, # only values >0 are included in osc:cca
                        count.max  = length(pixel_sign))
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
# clusters_sep is a matrix same size as data

stcs<- max(clusters_sep$cluster.count, na.rm = TRUE)
allcluster_max <- c()
clusters_sep$cluster.max <- vector(length = length(clusters_pos$cluster.count))
for (i in 1:length(clusters_sep$cluster.count)){ # retrieve maximum of each cluster
  clust_max <- data[clusters_sep$clusters==i] %>%
    max(.,na.rm = TRUE)
  clusters_sep$cluster.max[i] <- clust_max # assign each cluster its maximum
  allcluster_max <- c(allcluster_max, clust_max)
}
stcs_maxT_all <- max(allcluster_max, na.rm = TRUE)

# get_stcs modified ----
#***********************
get_stcs_mod<- function(data, alpha_local, null_distribution, data_dim){
  if(null_distribution == "normal") thr<- qnorm(1-alpha_local/2)
  if(null_distribution == "t") thr<- qt(1-alpha_local/2, df = data_dim[3]-2)

  pixel_sign<- sign(data)
  pixel_significant<- abs(data)>thr
  pixel_result<- pixel_sign*pixel_significant

  # positive
  pixel_result_pos<- pixel_result
  pixel_result_pos[is.na(pixel_result_pos)]<- -999
  clusters_pos<- osc::cca(pixel_result_pos,count.cells = TRUE, s=1, mode = 2, # only values >0 are included in osc:cca
                          count.max  = length(pixel_sign))
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

  stcs<- max(clusters_sep$cluster.count, na.rm = TRUE)
  stcs_idx<- which(length(clusters_sep$cluster.count)==stcs)
  stcs_cluster_results<- data[clusters_sep$clusters==stcs_idx] # retrieve all cells (by position in matrix?) that belong to the biggest cluster
  stcs_maxT<- max(stcs_cluster_results, na.rm = TRUE)

  # within cluster properties --- maxT works fine, all others are similar or worse
  stcs_mvt<- vector(length = length(clusters_sep$cluster.count), mode = "list")
  for (i in 1:length(clusters_sep$cluster.count)){
    #stcs_mvt[[i]]<- vector(length = 11, mode = "list")
    #print(i)

    #get results for cluster i and save for later
    cluster_results<- data[clusters_sep$clusters==i]
    stcs_mvt[[i]]$results<- cluster_results

    # maxT
    stcs_mvt[[i]]$maxT<- max(cluster_results, na.rm = TRUE)

    # avgT
    stcs_mvt[[i]]$meanT<- mean(cluster_results, na.rm = TRUE)

    # medianT
    stcs_mvt[[i]]$medianT<- median(cluster_results, na.rm = TRUE)

    # quantiles: 0.90, 0.95
    stcs_mvt[[i]]$q90T<- unname(quantile(cluster_results, probs = 0.90, na.rm = TRUE))
    stcs_mvt[[i]]$q95T<- unname(quantile(cluster_results, probs = 0.95, na.rm = TRUE))

    # average of top: 3, 5, 10 grid cells
    stcs_mvt[[i]]$meanTop3<- mean(head(sort(cluster_results, decreasing = TRUE), n=3), na.rm = TRUE)
    stcs_mvt[[i]]$meanTop5<- mean(head(sort(cluster_results, decreasing = TRUE), n=5), na.rm = TRUE)
    stcs_mvt[[i]]$meanTop10<- mean(head(sort(cluster_results, decreasing = TRUE), n=10), na.rm = TRUE)

    # average of top: 5% 10%
    stcs_mvt[[i]]$meanTop5percent<- mean(head(sort(cluster_results, decreasing = TRUE), n=length(cluster_results)*.05), na.rm = TRUE)
    stcs_mvt[[i]]$meanTop10percent<- mean(head(sort(cluster_results, decreasing = TRUE), n=length(cluster_results)*.10), na.rm = TRUE)
  }

  return(list(stcs=stcs, clusters=clusters_sep, stcs_maxT=stcs_maxT, stcs_mvt=stcs_mvt))
}

