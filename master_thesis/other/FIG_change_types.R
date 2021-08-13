#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************


# generate plots for different types of changes like mean, variance, slope etc.
#*******************************************************
set.seed(1234)

plot_lines <- function(n1, n2, sd1, sd2, mean1, mean2){
  segments(x0 = 0, y0 = mean1, x1 = n1-0.5, y1 = mean1, col = "black")
  segments(x0 = 0, y0 = mean1 + sd1, x1 = n1-0.5, y1 = mean1 + sd1, col = "black", lty="dotted")
  segments(x0 = 0, y0 = mean1 - sd1, x1 = n1-0.5, y1 = mean1 - sd1, col = "black", lty="dotted")
  segments(x0 = n1+0.5, y0 = mean2, x1 = n1+n2, y1 = mean2, col = "black")
  segments(x0 = n1+0.5, y0 = mean2+sd2, x1 = n1+n2, y1 = mean2+sd2, col = "black", lty="dotted")
  segments(x0 = n1+0.5, y0 = mean2-sd2, x1 = n1+n2, y1 = mean2-sd2, col = "black", lty="dotted")
}

change_meanvar <- function(n1, n2, sd1, sd2, mean1, mean2, title, omitx=FALSE) {
  change <- c(rnorm(n=n1, mean=mean1, sd=sd1), rnorm(n=n2, mean=mean2, sd=sd2))
  if(omitx){
    plot(change, type ="l", xlab = "", ylab = "y", main = title)
  } else {
    plot(change, type ="l", xlab = "Time", ylab = "y", main = title)
  }
  plot_lines(n1, n2, sd1, sd2, mean1, mean2)
  abline(v = n1, col="red", lty=2)
}

plot_slopes <- function(n1, n2, i1, i2, slope1, slope2){
  segments(x0 = 0, y0 = i1, x1 = n1-0.5, y1 = n1*slope1+i1, col = "black")
  segments(x0 = 0, y0 = i1 + 100, x1 = n1-0.5, y1 = n1*slope1+i1+100, col = "black", lty="dotted")
  segments(x0 = 0, y0 = i1 - 100, x1 = n1-0.5, y1 = n1*slope1+i1-100, col = "black", lty="dotted")
  segments(x0 = n1+0.5, y0 = i2, x1 = n1+n2, y1 = (n2)*slope2+i2, col = "black")
  segments(x0 = n1+0.5, y0 = i2+100, x1 = n1+n2, y1 = (n2)*slope2+i2+100, col = "black", lty="dotted")
  segments(x0 = n1+0.5, y0 = i2-100, x1 = n1+n2, y1 = (n2)*slope2+i2-100, col = "black", lty="dotted")
}

change_intertrend <- function(n1, n2, i1, i2, slope1, slope2, title, omitx=FALSE){
  change <- c((slope1*c(1:n1) + i1 + rnorm(n=n1, mean=0, sd=100)), (slope2*c(1:n2) + i2 + rnorm(n=n2, mean=0, sd=100)))
  if(omitx){
    plot(change, type ="l", xlab = "", ylab = "y", main = title)
  } else {
    plot(change, type ="l", xlab = "Time", ylab = "y", main = title)
  }
  plot_slopes(n1, n2, i1, i2, slope1, slope2)
  abline(v = n1, col="red", lty=2)
}

change_pos_intertrend <- function(n1, n2, c){
  #par(mfrow = c(3, 1))
  change_intertrend(n1, n2, 5, 200, 3, 3, "Change in intercept", omitx = TRUE)
  change_intertrend(n1, n2, 5, -20, -1, 1.5, "Change in intercept and trend",omitx = TRUE, omit_title = TRUE)
  autocorr(n1, n2, c)
  #par(mfrow = c(1,1))
}

# Plot
#***************************
pdf_fn <- "change_types.pdf"
pdf(pdf_fn, height=7.6, width=10)
par(mfrow = c(3, 2))
n1 <- 30
n2 <- 30
change_meanvar(n1, n2, 1.5, 1.5, -1, 3, "Change in Mean", omitx = TRUE)
change_intertrend(n1, n2, 5, 290, 3, 3, "Change in Intercept", omitx = TRUE)
change_meanvar(n1, n2, 1.5, 4, 0, 0, "Change in Variance")
change_intertrend(n1, n2, 5, 130, 4, -2, "Change in Slope", omitx = TRUE)
change_meanvar(n1, n2, 2, 1, 2, 0, "Change in Mean and Variance", omitx = TRUE)
change_intertrend(n1, n2, 5, 20, -3, 4, "Change in Intercept and Slope")
dev.off()
