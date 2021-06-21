# script to send perm_dist to cluster

suppressMessages(library(devtools))
suppressMessages(library(tidyverse))
devtools::load_all("/home/veronika/CPD/change_point_detection/")
source("/home/veronika/CPD/change_point_detection/testing/perm_dist_SLURM_2d.R")

fx=mcusum_function
alpha_local=0.1
alpha_global=0.1
null_distribution <- "p-values"
seed=9
block_size=NULL
verbose=TRUE
#nperm = 500
nperm = 4


# ################################ BU data ####################################

data_lai <- readRDS("/home/veronika/CPD/data/10000int_BU_LAI_yearlymean_1982_2018")
data <- data_lai[3580:3609, 380:409,]

res <- perm_dist_SLURM_2d(data=data_lai, fx=fx, nperm=nperm, alpha_local=alpha_local,
                       alpha_global=alpha_global, null_distribution=null_distribution,
                       seed=seed, block_size=NULL, verbose=TRUE)
################################ BU data ####################################
