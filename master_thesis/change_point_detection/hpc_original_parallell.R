# try parallelization of modified cusum on unpermuted image

#devtools::load_all("/home/veronika/CPD/change_point_detection/")
devtools::load_all()
library(osc)

data_lai <- readRDS("/home/veronika/CPD/data/10000int_BU_LAI_yearlymean_1982_2018")
data <- array_to_matrix(data_lai)
rm(data_lai)

img <- data$Y
# conversion to list
img_list <- lapply(seq_len(ncol(img)), function(i) img[,i])

tmp_fn <- function(d){
  #source("/home/veronika/CPD/change_point_detection/R/modified_cusum.R")
  source("master_thesis/R/modified_cusum.R")
  res <- mcusum_function(d, loc = T) # loc = T returns also timing and ar order estimates
  cat("computed modified cusum")
  return(res)
}

library(clustermq)
t <- system.time({
results<- Q(tmp_fn,
            d=img_list,
            n_jobs = 200,
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

saveRDS(list(fin,t),"/results/MCUSUM_BU_orig_pvals_ar_locs.rds")
