# script to send perm_dist to cluster

suppressMessages(library(devtools))
suppressMessages(library(tidyverse))
devtools::load_all("/home/veronika/CPD/change_point_detection/")
#source("/home/veronika/CPD/change_point_detection/testing/perm_dist_SLURM.R")
source("/home/veronika/CPD/change_point_detection/testing/perm_dist_SLURM_2d.R")

#fx=cusum_pval_function
#fx=cusum_function
fx=mosum_function
alpha_local=0.1
alpha_global=0.1
#null_distribution <- "p-values"
#null_distribution <- "brownian_motion"
null_distribution <- "brownian_bridge_increments"
seed=9
block_size=NULL
verbose=TRUE
nperm = 1000
#nperm = 4


# ################################ BU data ####################################

data_lai <- readRDS("/home/veronika/CPD/data/10000int_BU_LAI_yearlymean_1981_2018")
#data_lai <- data_lai[1700:2300, 800:1400,]

res <- perm_dist_SLURM_2d(data=data_lai, fx=fx, nperm=nperm, alpha_local=alpha_local,
                       alpha_global=alpha_global, null_distribution=null_distribution,
                       seed=seed, block_size=NULL, verbose=TRUE)
################################ BU data ####################################

################################ LTDR data 2d ####################################
# 
# data_lai <- readRDS("/home/veronika/CPD/data/10000int_LTDR_LAI_yearlymean_1981_2018.rds")
# 
# res <- perm_dist_SLURM_2d(data=data_lai, fx=fx, nperm=nperm, alpha_local=alpha_local,
#                           alpha_global=alpha_global, null_distribution=null_distribution,
#                           seed=NULL, block_size=NULL, verbose=TRUE)
################################ LTDR data 2d #####################################

################################ NOAA data 2d ####################################

# data_lai <- readRDS("/home/veronika/CPD/data/10000int_NOAA_LAI_yearlymean_1981_2018.rds")
# 
# res <- perm_dist_SLURM_2d(data=data_lai, fx=fx, nperm=nperm, alpha_local=alpha_local,
#                           alpha_global=alpha_global, null_distribution=null_distribution,
#                           seed=NULL, block_size=NULL, verbose=TRUE)
################################ NOAA data 2d ####################################


