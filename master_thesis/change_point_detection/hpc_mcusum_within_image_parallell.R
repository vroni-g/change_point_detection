#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# !! Input Data is not available in GitHub Repository !!

# Within Image Parallelization for MMCUSUM change point detection
#***************************************************************

#devtools::load_all("/home/veronika/CPD/change_point_detection/")
devtools::load_all()
library(osc)
library(clustermq)

data_lai <- readRDS("/home/veronika/CPD/data/10000int_BU_LAI_yearlymean_1982_2018")
data <- array_to_matrix(data_lai)
rm(data_lai)

seed=9
block_size=NULL
verbose=TRUE

# pre-specified permutation matrix with 100 perms
perm_mat<- readRDS("master_thesis/data_preprocessing/perm_matrix.rds")

tmp_fn <- function(d){
  #source("/home/veronika/CPD/change_point_detection/R/modified_cusum.R")
  source("R/modified_cusum.R")
  res <- mcusum_function(d)
  cat("computed modified cusum")
  return(res)
}

for(i in 1:30){ # specify the indices for permutation derived from the matrix
  img <- data$Y[perm_mat[i,],]
  img_list <- lapply(seq_len(ncol(img)), function(i) img[,i])

  t <- system.time({
    results<- Q(tmp_fn,
                d=img_list,
                n_jobs = 200,
                export = list(mcusum_function = mcusum_function),
                template = list(job_name = "MC_pl_V",
                                partition = "all",
                                log_file = paste0("/home/veronika/CPD/logs/Parallel_MCUSUM_200j_900mem_%a.txt"),
                                memory = 900,
                                n_cpus = 1),
                fail_on_error = FALSE,
                verbose = TRUE)

    img2 <- matrix(unlist(results), ncol = length(img_list))
    data_info <- data[2:5]
    fin <- matrix(NA, ncol = data_info$ncol, nrow = data_info$nrow)
    fin[data_info$wh.sel]<- img2
  })

  saveRDS(list(fin,t),paste0("master_thesis/results/BU_MCUSUM_permmat_", i,".rds"))

}
