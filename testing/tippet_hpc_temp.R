# script to send perm_dist to cluster
suppressMessages(library(devtools))
suppressMessages(library(tidyverse))
devtools::load_all()
source("testing/perm_dist_SLURM.R")

data=temp_gistemp
fx=sample_mk_function
alpha_local=0.05
alpha_global=0.05
null_distribution <- "normal"
seed=NULL
block_size=NULL
verbose=TRUE
nperm = 1000

sen0 <- function(y,x){
  zyp.slopediff <- function(i, xx, yy, n) (yy[1:(n - i)] - yy[(i + 1):n])/(xx[1:(n - i)] - xx[(i + 1):n])
  n <- length(y)
  if (missing(x)) x <- c(1:n)
  slopes <- unlist(lapply(1:(n - 1), zyp.slopediff, x, y, n))
  return(median(slopes[is.finite(slopes)], na.rm=TRUE))
}


data_detrend<- data %>% apply(1:2, # apply(1:2,...) will apply function to every cell
                              function(x)
                              {
                                (x- 1:length(x)*sen0(x))
                              }
)

data_detrend <- aperm(data_detrend, c(2,3,1)) # transpose it to put lat & long in the first dimensions again
rm(data)
res <- perm_dist_SLURM(data=data_detrend, fx=fx, nperm=nperm, alpha_local=alpha_local,
                       alpha_global=alpha_global, null_distribution=null_distribution,
                       seed=NULL, block_size=NULL, verbose=TRUE)
filename <- paste0("testing/detrend_temp_tippet_nperm_", nperm, ".rds")
saveRDS(res, file = filename)

