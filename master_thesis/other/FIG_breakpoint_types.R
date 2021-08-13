#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# generate plots for different type of breakpoints 
#*******************************************************
set.seed(1234)

plot_slopes <- function(n1, n2, i1, i2, slope1, slope2){
  segments(x0 = 0, y0 = i1, x1 = n1, y1 = n1*slope1+i1, col = "gray57")#, lty="dotted", lwd = 3)
  #segments(x0 = 0, y0 = i1 + 50, x1 = n1, y1 = n1*slope1+i1+50, col = "black", lty="dotted")
  #segments(x0 = 0, y0 = i1 - 50, x1 = n1, y1 = n1*slope1+i1-50, col = "black", lty="dotted")
  segments(x0 = n1+1, y0 = i2, x1 = n1+n2, y1 = (n2)*slope2+i2, col = "gray57")#, lty="dotted", lwd = 3)
  #segments(x0 = n1+1, y0 = i2+50, x1 = n1+n2, y1 = (n2)*slope2+i2+50, col = "black", lty="dotted")
  #segments(x0 = n1+1, y0 = i2-50, x1 = n1+n2, y1 = (n2)*slope2+i2-50, col = "black", lty="dotted")
}

change_intertrend <- function(n1, n2, i1, i2, slope1, slope2, title, clr,omitx=FALSE){
  change <- c((slope1*c(1:n1) + i1 + rnorm(n=n1, mean=0, sd=100)), (slope2*c(1:n2) + i2 + rnorm(n=n2, mean=0, sd=100)))
  if(omitx){
    plot(change, type ="l", xlab = "", ylab = "y", main = title)
  } else {
    plot(change, type ="l", xlab = "Time", ylab = "y", main = title)
  }
  plot_slopes(n1, n2, i1, i2, slope1, slope2)
  abline(v = n1+0.5, col="red", lty=2)
  rect(n1+n2-10, par("usr")[[4]]-40, n1+n2+2, par("usr")[[4]]-10, col = clr)
}


# Plot
#***************************
n1 <- 30
n2 <- 30
pdf_fn <- "breakpoint_types.pdf"
pdf(pdf_fn, height=12, width=10)
par(mfrow = c(4, 2))

# browning to greening (turquiose): #0FEAD2
change_intertrend(40, 40, 10, -180, -5, 2.5, "Browning to Greening", omitx = TRUE, clr = "#0FEAD2")
# greening to browning (yellow or orange): #F7F121
change_intertrend(50, 30, 5, 230, 5, -3, "Greening to Browning", omitx = TRUE, clr = "#F7F121")
# continued greening (dark green): #145607
change_intertrend(25, 55, 5, 120, 6, 3, "Continued Greening", omitx = TRUE, clr = "#145607")
# continued browning (dark blue):#181C67
change_intertrend(55, 25, 5, -160, -3, -6, "Continued Browning", omitx = TRUE, clr = "#181C67")
# greening onset (middle green): #138211
change_intertrend(30, 50, 5, 10, 0, 3.5, "Greening Onset", omitx = TRUE, clr = "#138211")
# browning onset (middle blue): #1724FB
change_intertrend(40, 40, 5, 5, 0, -3, "Browning Onset", omitx = TRUE, clr = "#1724FB")
# stalled greening (light green): #63CE47
change_intertrend(55, 25, 5, 210, 4, 0, "Stalled Greening", clr = "#63CE47")
# stalled browning (light blue): #8DB5F1
change_intertrend(40, 40, 5, -150, -3, 0, "Stalled Browning",clr = "#8DB5F1")
# non non (gray): gray77

dev.off()

