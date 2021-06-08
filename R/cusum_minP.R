#' Recursive CUSUSM structual change test function with pre-whitening
#' 
#' Computes the recursive residuals for the given time series and returns the test statistic
#'
#' @param x time series of a single grid cell
#' @return Recursive residuals CUSUM test statistic
#' @export cusum_stat

cusum_pval <- function(x){
  res <- ts(x, frequency=1, start=c(1980)) %>%
    strucchange::efp(. ~ time(.), ., type = 'Rec-CUSUM') %>%
    strucchange::sctest()
  if(rlang::is_empty(res$p.value)) return(NA)
  if(is.na(res$p.value)) return(NA)
  return(res$p.value)
}

#' This function is an example of a valid function as an input to the multiple
#' testing correction function. Autocorrelation lag 1 , r, is estimated and
#' removed from each grid cell time series ,x_t, so thatso that the new time
#' series is y_t = x_t - rx_t-1
#' @param x time series at a single grid cell
#' @return recursive CUSUM test statistic, corrected for temporal autocorrelation
#' @export cusum_function

cusum_pval_function<- function(x){
  if(any(is.na(x))) x<- x[!is.na(x)]
  #if(length(x)<8) return(NA)

  x <- x[is.finite(x)]
  xn <- (x[-1] - (x[-length(x)] * rk_fn(x)))
  p <- cusum_pval(xn)

  if(is.na(p)) return(NA)
  #if(rlang::is_empty(s)) return(NA)

  return(p)
}
