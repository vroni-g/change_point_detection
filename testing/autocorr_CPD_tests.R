# check false positive rate of CUSUM method applied to autocorrelated noise data
#***********************
library(strucchange)
library(tidyverse)

coeffs <- seq(0, 0.9, by=0.1)
fprs <- vector("list", 11)
i <- 1

for (c in coeffs){
  p <- c()
  n <- 10000
  for (i in 1:n){
    d <- arima.sim(model = list(order = c(1, 0, 0), ar = c), n = 250)
    res <-  efp(d ~ 1, d, type = 'Rec-CUSUM') %>%
      sctest()
    p <- c(p, res$p.value)
  }
  fpr <- sum(p<0.05)/n
  cat(paste('False positive rate for Rec-CUSUM with autocorrelation ', c,' is: ', fpr, '\n'))
  fprs[i] <- c(c, fpr)
  i <- i + 1
}

# for n = 10000:
# False positive rate for Rec-CUSUM with autocorrelation  0  is:  0.0474
# False positive rate for Rec-CUSUM with autocorrelation  0.1  is:  0.082
# False positive rate for Rec-CUSUM with autocorrelation  0.2  is:  0.1218
# False positive rate for Rec-CUSUM with autocorrelation  0.3  is:  0.2006
# False positive rate for Rec-CUSUM with autocorrelation  0.4  is:  0.2862
# False positive rate for Rec-CUSUM with autocorrelation  0.5  is:  0.3883
# False positive rate for Rec-CUSUM with autocorrelation  0.6  is:  0.529
# False positive rate for Rec-CUSUM with autocorrelation  0.7  is:  0.6743
# False positive rate for Rec-CUSUM with autocorrelation  0.8  is:  0.8249
# False positive rate for Rec-CUSUM with autocorrelation  0.9  is:  0.9544
