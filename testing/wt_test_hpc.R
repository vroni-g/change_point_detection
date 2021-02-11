# script to send perm_dist to cluster
suppressMessages(library(devtools))
suppressMessages(library(tidyverse))
devtools::load_all()
source("testing/perm_dist_SLURM.R")

#data=temp_gistemp
fx=sample_mk_function
alpha_local=0.05
alpha_global=0.05
null_distribution <- "normal"
seed=NULL
block_size=NULL
verbose=TRUE
nperm = 100

# sen0 <- function(y,x){
#   zyp.slopediff <- function(i, xx, yy, n) (yy[1:(n - i)] - yy[(i + 1):n])/(xx[1:(n - i)] - xx[(i + 1):n])
#   n <- length(y)
#   if (missing(x)) x <- c(1:n)
#   slopes <- unlist(lapply(1:(n - 1), zyp.slopediff, x, y, n))
#   return(median(slopes[is.finite(slopes)], na.rm=TRUE))
# }
#
#
# data_detrend<- data %>% apply(1:2, # apply(1:2,...) will apply function to every cell
#                               function(x)
#                               {
#                                 (x- 1:length(x)*sen0(x))
#                               }
# )
#
# data_detrend <- aperm(data_detrend, c(2,3,1)) # transpose it to put lat & long in the first dimensions again
# data=data_detrend

load("/Users/veronikagrupp/hpc_vroni/home/jose/LAI/data/CHEN_RANGA/AVHRR/yearly_mean/lai_data.RData")
tibble_list_to_3d_array<- function(data){
  library(raster)
  out<- matrix(NA, ncol = data$lat %>% unique %>% length, nrow = data$lon %>% unique %>% length)
  data_long<- data %>% dplyr::select(data) %>% ungroup %>% unnest(cols = data)
  n_years<- data_long$t %>% unique %>% length
  data_list<- vector(mode = "list", length = unique(data_long$t) %>% length)
  i<- 1
  for(year in unique(data_long$t)){
    data_raster<- data_long %>% filter(t==year) %>% rename(x = lon, y = lat) %>%
      dplyr::select(x, y, lai) %>% rasterFromXYZ(crs = CRS("+init=epsg:4326"))
    data_list[[i]]<- data_raster$lai %>% as.matrix
    i<- i+1
  }
  data_array<- array(NA, dim= c(dim(data_list[[1]]), length(data_list)))
  for(i in 1:n_years) data_array[,,i]<- data_list[[i]]
  return(data_array)
}

data_lai<- tibble_list_to_3d_array(data)
rm(data)
# filename <- paste0("/home/veronika/CPD/data/lai_yearlymean_3d.rds")
# saveRDS(data_lai, file = filename)
# data <- readRDS("/home/veronika/CPD/data/lai_yearlymean_3d.rds")
res <- perm_dist_SLURM(data=data_lai, fx=fx, nperm=nperm, alpha_local=alpha_local,
                       alpha_global=alpha_global, null_distribution=null_distribution,
                       seed=NULL, block_size=NULL, verbose=TRUE)
filename <- paste0("testing/LAI_Wtadjust_nperm_", nperm, ".rds")
saveRDS(res, file = filename)
