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
