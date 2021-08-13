#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************


# replicate plot of Beaulieu 2018 p.9521 (https://doi.org/10.1175/JCLI-D-17-0863.1)
# to show possible pitfalls in change point detection in presence of spatial autocorrelation
#************************************************************

library(tidyverse)
# naming convention: true_fitted
pdf("misuse_cases.pdf", width =13 , height = 11)
par(mfrow=c(3,2))
set.seed(1234)
# mean_trend 
#****************************
m1 <- 5
m2 <- 15
n1 <- 30
n2 <- 30
l <- n1+n2
sd = 3
meanchange <- c(rnorm(n=n1, mean=m1, sd=sd), rnorm(n=n2, mean=m2, sd=sd))
plot(meanchange, type = "l", lwd = 1.5, axes = FALSE, xlab = "", ylab = "", main = "a)",adj = 0, cex.main = 2)
rect(xleft = 0, ybottom = m1-sd, xright = n1, ytop = m1+sd, border = NA, col = rgb(0.5,0.5,0.5,1/4))
rect(xleft = n1, ybottom = m2-sd, xright = l, ytop = m2+sd, border = NA, col = rgb(0.5,0.5,0.5,1/4))
segments(x0 = 0, y0 = m1, x1 = n1, y1 = m1, col = "black", lty="longdash", lwd = 2.5)
segments(x0 = n2, y0 = m2, x1 = n1+n2, y1 = m2, col = "black", lty="longdash", lwd = 2.5)
abline(lm(meanchange ~ c(1:l)), col = "blue", lwd = 3.5)
#***************************

# trend_mean
#****************************
slope <- 4
i <- 3
sdt = 70
trend <- slope*c(1:l) + i + rnorm(n=l, mean=0, sd=sdt)
plot(trend, type = "l", lwd = 1.5,axes = FALSE, xlab = "", ylab = "", main = "b)",adj = 0, cex.main = 2)
x <- c(0,0,l,l)
y <- c(i-sdt, i+sdt, slope*l+i+sdt, slope*l+i-sdt)
polygon(x, y, border = NA, col = rgb(0.5,0.5,0.5,1/4))
abline(lm(trend ~ c(1:l)), col = "black", lty="longdash", lwd = 2.5)
segments(x0 = 0, y0 = 20*i, x1 = 1.7*l/3, y1 = 20*i, col = "blue", lwd = 3.5)
segments(x0 = 1.7*l/3, y0 = 65*i, x1 = l, y1 = 65*i, col = "blue", lwd = 3.5)

#****************************


# AR(1)_trend
#****************************
ar1 <- arima.sim(n = l, model = list(order = c(1, 0, 0), ar = 0.9), sd = 2)
plot(ar1, type = "l", lwd = 1.5,axes = FALSE, xlab = "", ylab = "", col = "red", main = "c)",adj = 0, cex.main = 2)
rect(xleft = 0, ybottom = 0, xright = l, ytop = 9, border = NA, col = rgb(1,0,0,1/4))
abline(h=mean(ar1), col = "black", lty="longdash", lwd = 2.5)
abline(lm(ar1 ~ c(1:l)),col = "blue", lwd = 3.5)

#****************************


# AR(1)_meanshift
#****************************
ar1_m <- arima.sim(n = l, model = list(order = c(1, 0, 0), ar = 0.9), sd = 2)
plot(ar1_m, type = "l", lwd = 1.5,axes = FALSE, xlab = "", ylab = "", col = "red",  main = "d)",adj = 0, cex.main = 2)
rect(xleft = 0, ybottom = mean(ar1_m)-6, xright = l, ytop = mean(ar1_m)+6, border = NA, col = rgb(1,0,0,1/4))
abline(h=mean(ar1_m), col = "black", lty="longdash", lwd = 2.5)
segments(x0 = 0, y0 = -7.5, x1 = l/3, y1 = -7.5, col = "blue", lwd = 3.5)
segments(x0 = l/3, y0 = mean(ar1_m), x1 = 2*l/3, y1 = mean(ar1_m), col = "blue", lwd = 3.5)
segments(x0 = 2*l/3, y0 = 2, x1 = l, y1 = 2, col = "blue", lwd = 3.5)

#****************************


# meanshift_AR(1)
#****************************
fit <- forecast::Arima(meanchange, order = c(1, 0, 0))
plot(meanchange, type = "l", lwd = 1.5, axes = FALSE, xlab = "", ylab = "", main = "e)",adj = 0, cex.main = 2)
rect(xleft = 0, ybottom = m1-sd, xright = n1, ytop = m1+sd, border = NA, col = rgb(0.5,0.5,0.5,1/4))
rect(xleft = n1, ybottom = m2-sd, xright = l, ytop = m2+sd, border = NA, col = rgb(0.5,0.5,0.5,1/4))
segments(x0 = 0, y0 = m1, x1 = n1, y1 = m1, col = "black", lty="longdash", lwd = 2.5)
segments(x0 = n2, y0 = m2, x1 = n1+n2, y1 = m2, col = "black", lty="longdash", lwd = 2.5)
lines(fitted(fit),col="blue", lwd = 3.5)
#****************************
cx <- 2
plot(0,type ="n", axes = FALSE, xlab = "", ylab = "")
segments(x0 = 0.78, y0 = 0.7, x1 = 0.87, y1 = 0.7, col = "blue", lwd = 3.5)
text(0.9, 0.65, " fit", cex = cx, adj = c(0,0))
segments(x0 = 0.78, y0 = 0.45, x1 = 0.87, y1 = 0.45, col = "black", lty="longdash", lwd = 2.5)
text(0.9, 0.4, " true signal", cex = cx, adj = c(0,0))
segments(x0 = 0.78, y0 = 0.2, x1 = 0.87, y1 = 0.2, col = "black", lwd = 1.5)
text(0.9, 0.15, " white noise + signal time series", cex = cx, adj = c(0,0))
rect(0.78, -0.005, 0.87, -0.12, border = NA, col = rgb(0.5,0.5,0.5,1/4))
text(0.9, -0.1, " white noise variability", cex = cx, adj = c(0,0))
segments(x0 = 0.78, y0 = -0.3, x1 = 0.87, y1 = -0.3, col = "red", lwd = 1.5)
text(0.9, -0.35, " AR(1) time series", cex = cx, adj = c(0,0))
rect(0.78, -0.505, 0.87, -0.62, border = NA, col = rgb(1,0,0,1/4))
text(0.9, -0.6, " AR(1) variability", cex = cx, adj = c(0,0))

dev.off()

