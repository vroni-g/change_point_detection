#' CUSUSM structural change test with recursive residuals
#'
#' Computes the recursive residuals for the given time series and returns the test statistic
#'
#' @param x time series of a single grid cell
#' @return Recursive residuals CUSUM test statistic
#' @export cusum_stat

cusum_stat <- function(x){
  res <- ts(x, frequency=1, start=c(1980)) %>%
    strucchange::efp(. ~ time(.), ., type = 'Rec-CUSUM') %>%
    strucchange::sctest()
  if(rlang::is_empty(res$statistic)) return(NA)
  if(is.na(res$statistic)) return(NA)
  return(res$statistic)
}
