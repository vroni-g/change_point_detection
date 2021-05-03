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
nperm = 1000
#nperm = 2


# ################################ BU data 2d ####################################
# 
# data_lai <- readRDS("/home/veronika/CPD/data/10000int_BU_LAI_yearlymean_1981_2018")
# data_lai_2d<- array_to_matrix(data_lai)
# 
# res <- perm_dist_SLURM_2d(data=data_lai_2d, fx=fx, nperm=nperm, alpha_local=alpha_local,
#                        alpha_global=alpha_global, null_distribution=null_distribution,
#                        seed=NULL, block_size=NULL, verbose=TRUE)
################################ BU data 2d ####################################

################################ LTDR data 2d ####################################

data_lai <- readRDS("/home/veronika/CPD/data/10000int_LTDR_LAI_yearlymean_1981_2018.rds")
data_lai_2d<- array_to_matrix(data_lai)

res <- perm_dist_SLURM_2d(data=data_lai_2d, fx=fx, nperm=nperm, alpha_local=alpha_local,
                          alpha_global=alpha_global, null_distribution=null_distribution,
                          seed=NULL, block_size=NULL, verbose=TRUE)
################################ LTDR data 2d ####################################

