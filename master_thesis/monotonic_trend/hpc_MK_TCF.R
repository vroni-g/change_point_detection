#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# !! Input Data is not available in GitHub Repository !!


# Script to send Mann Kendall trend test with permutation based multiple testing
# correction to cluster
#************************************

suppressMessages(library(devtools))
suppressMessages(library(tidyverse))
devtools::load_all()
source("master_thesis/monotonic_trend/parallel_MK_TCF_2d.R")

fx=sample_mk_function
alpha_local=0.1
alpha_global=0.1
null_distribution <- "normal"
seed=9
block_size=NULL
verbose=TRUE
nperm = 1000
res_dir = "home/veronika/CPD/results/"



data_lai <- readRDS("/home/veronika/CPD/data/10000int_BU_LAI_yearlymean_1981_2018") # data from HPC

res <- perm_dist_SLURM_2d(data=data_lai, fx=fx, nperm=nperm, alpha_local=alpha_local,
                       alpha_global=alpha_global, null_distribution=null_distribution,
                       seed=seed, block_size=NULL, verbose=TRUE, res_dir = res_dir)



