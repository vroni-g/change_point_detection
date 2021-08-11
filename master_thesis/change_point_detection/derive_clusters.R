#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# !! Input Data is not available in GitHub Repository !!


# derive clusters for all permuted images MCUSUM results
#******************************

get_res <- function(i){
  devtools::load_all("/home/veronika/CPD/change_point_detection/")
  res <- readRDS(paste0("/home/veronika/CPD/MCUSUM/results/BU_MCUSUM_permmat_", i,".rds"))
  pmat <- res[[1]]
  pmat <- matrix(as.numeric(unlist(pmat)),nrow=nrow(pmat))
  minP <- min(as.vector(pmat),na.rm = T)
  tmp_stcs <- get_stcs(data=pmat, alpha_local=0.1, null_distribution = "p-values",tippet=T)
  stcs<- tmp_stcs$stcs
  stcs_minP <- tmp_stcs$stcs_minP
  peak_intensity <- tmp_stcs$peak_intensity
  r <- list(c(minP = minP, stcs = stcs, stcs_minP = stcs_minP, peak_intensity=peak_intensity), tmp_stcs$clusters)
  saveRDS(r, paste0("/home/veronika/CPD/MCUSUM/results/cluster/BU_MCUSUM_cluster_", i,".rds"))
  cat("finished file ", i, "\n")
}

library(clustermq)
results<- Q(get_res,
            i=1:30,
            n_jobs = 30,
            template = list(job_name = "cluster_Mcusum",
                            partition = "all",
                            log_file = paste0("/home/veronika/CPD/logs/BU_MCUSUM_getcluster.txt"),
                            memory = 10000,
                            n_cpus = 1),
            fail_on_error = FALSE,
            verbose = TRUE)

