#' CUSUSM structural change test with recursive residuals
#'
#' Computes the recursive residuals for the given time series and returns the test statistic
#'
#' @param x time series of a single grid cell
#' @return Recursive residuals CUSUM test statistic
#' @export cusum_stat

cusum_stat <- function(x){
  if(any(is.na(x))) x<- x[!is.na(x)]
  if(length(x)<8) return(NA)

  res <- ts(x, frequency=1, start=c(1980)) %>%
    strucchange::efp(. ~ 1, ., type = 'Rec-CUSUM') %>%
    strucchange::sctest()
  return(res$statistic)
}

x <- rnorm(2*20,0,.2) %>%
  ts(., frequency=1, start=c(1980))

setwd("/home/veronika/CPD/change_point_detection")
devtools::load_all()

# create random matrix
mat <- function(nrow){
  m <- matrix(rnorm(1200, 0,.2),nrow=nrow, ncol = 30)
  return(m)
}
v <- rep(40, 10)
mat_list <- lapply(v, mat)
data <- abind::abind( matlist, along=3 )

#data <- readRDS("/home/veronika/CPD/data/NOAA_LAI/yearly_median/masked/int10000_NOAA_LAI_masked_median_1981_2020.rds")

perm_matrix<- perm_matrix(nobs = dim(data)[3], nperm = 2, block_size = NULL, seed = NULL)

tmp <- apply(data[,,perm_matrix[2,]], 1:2, cusum_stat)
