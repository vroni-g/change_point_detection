# script to send perm_dist to cluster

suppressMessages(library(devtools))
suppressMessages(library(tidyverse))
devtools::load_all("/home/veronika/CPD/change_point_detection/")
#source("/home/veronika/CPD/change_point_detection/testing/perm_dist_SLURM.R")
source("/home/veronika/CPD/change_point_detection/testing/perm_dist_SLURM_2d.R")

fx=sample_mk_function
alpha_local=0.05
alpha_global=0.05
null_distribution <- "normal"
seed=NULL
block_size=NULL
verbose=TRUE
nperm = 250
#nperm = 2

#data_lai <- readRDS("/home/veronika/CPD/data/NOAA_LAI/yearly_median/masked/int10000_NOAA_LAI_masked_median_1981_2020.rds")

data_lai <- readRDS("/home/veronika/CPD/data/NOAA_LAI/yearly_median/masked/int10000_NOAA_LAI_masked_median_1981_2019.rds")


# res <- perm_dist_SLURM(data=data_lai,  nperm=nperm, alpha_local=alpha_local, #fx=fx,
#                        alpha_global=alpha_global, null_distribution=null_distribution,
#                        seed=NULL, block_size=NULL, verbose=TRUE)
#
# filename <- paste0("/home/veronika/CPD/results/NOAA_LAI_tippet_nperm_", nperm, ".rds")
# saveRDS(res, file = filename)

################################ 2d approach ####################################
data_lai_2d<- array_to_matrix(data_lai)

res <- perm_dist_SLURM_2d(data=data_lai_2d, fx=fx, nperm=nperm, alpha_local=alpha_local,
                       alpha_global=alpha_global, null_distribution=null_distribution,
                       seed=NULL, block_size=NULL, verbose=TRUE)
################################ 2d approach ####################################

