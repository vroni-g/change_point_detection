# script to send perm_dist to cluster
#setwd("/home/veronika/CPD/change_point_detection")
suppressMessages(library(devtools))
suppressMessages(library(tidyverse))
devtools::load_all("/home/veronika/CPD/change_point_detection/")
source("/home/veronika/CPD/change_point_detection/testing/perm_dist_SLURM.R")
fx=sample_mk_function
alpha_local=0.05
alpha_global=0.05
null_distribution <- "normal"
seed=NULL
block_size=NULL
verbose=TRUE
nperm = 1000
#nperm = 3

#data_lai <- readRDS("/home/veronika/CPD/data/NOAA_LAI/yearly_median/masked/int10000_NOAA_LAI_masked_median_1981_2020.rds")
data_lai <- readRDS("/home/veronika/CPD/data/NOAA_LAI/yearly_median/masked/int10000_NOAA_LAI_masked_median_1981_2019.rds")
# does it matter for CUSUM if values are multiplied by 10000?
# print(object.size(data_lai),units="auto")
# library(Matrix)
# data_lai_sparse<- Matrix(data_lai, sparse = TRUE)
# print(object.size(data_lai_sparse),units="auto")
#

res <- perm_dist_SLURM(data=data_lai,  nperm=nperm, alpha_local=alpha_local, #fx=fx,
                       alpha_global=alpha_global, null_distribution=null_distribution,
                       seed=NULL, block_size=NULL, verbose=TRUE)

filename <- paste0("/home/veronika/CPD/results/NOAA_LAI_tippet_nperm_", nperm, ".rds")
saveRDS(res, file = filename)

################################ 2d approach ####################################
data_lai_2d<- array_to_matrix(data_lai)

res <- perm_dist_SLURM_2d(data=data_lai_2d, fx=fx, nperm=nperm, alpha_local=alpha_local,
                       alpha_global=alpha_global, null_distribution=null_distribution,
                       seed=NULL, block_size=NULL, verbose=TRUE)
################################ 2d approach ####################################




#**************************************
# BU LAI Data
#**************************************

# load("/home/jose/LAI/data/CHEN_RANGA/AVHRR/yearly_mean/lai_data.RData")
# tibble_list_to_3d_array<- function(data){
#   library(raster)
#   out<- matrix(NA, ncol = data$lat %>% unique %>% length, nrow = data$lon %>% unique %>% length)
#   data_long<- data %>% dplyr::select(data) %>% ungroup %>% unnest(cols = data)
#   n_years<- data_long$t %>% unique %>% length
#   data_list<- vector(mode = "list", length = unique(data_long$t) %>% length)
#   i<- 1
#   for(year in unique(data_long$t)){
#     data_raster<- data_long %>% filter(t==year) %>% rename(x = lon, y = lat) %>%
#       dplyr::select(x, y, lai) %>% rasterFromXYZ(crs = CRS("+init=epsg:4326"))
#     data_list[[i]]<- data_raster$lai %>% as.matrix
#     i<- i+1
#   }
#   data_array<- array(NA, dim= c(dim(data_list[[1]]), length(data_list)))
#   for(i in 1:n_years) data_array[,,i]<- data_list[[i]]
#   return(data_array)
# }
# data_lai<- tibble_list_to_3d_array(data)
# rm(data)
