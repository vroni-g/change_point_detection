# OPTION 1
# Emprical bivariate cdf
library("bivariate")

stcs<- perm_results$stcs
maxT<- perm_results$maxT
biv<- ebvcdf(stcs, maxT)
plot(biv)
plot(biv, FALSE)
plot(biv, FALSE, fb = c(.90,.95, .99))
abline(v = quantile(stcs, c(.90, .95, .99), names = FALSE), col = "red")
biv(1000, 5)

mf.Fh(biv)

plot(median(stcs), mean(maxT))
?plot.bv


bvmat(biv)

fbiv2<- nbvpdf(stcs, maxT)
plot(biv2)

# OPTION 2
# .95 quantile contour line


library(ggplot2)


get_contour<- function(x, y, alpha){
  library(ks)
  d <- data.frame(x=x,y=y)
  kd <- ks::kde(d, compute.cont=TRUE)
  contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                      z=estimate, levels=cont[paste0(alpha*100, "%")])[[1]])
  data.frame(contour_95)[,-1]
}
set.seed(1001)
d <- data.frame(x=tmp_stcs, y=tmp_maxt)

kd <- ks::kde(d, H=hscv1, compute.cont=TRUE)
contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["5%"])[[1]])
contour_95 <- data.frame(contour_95)

data_frame(stcs=tmp_stcs, maxT=tmp_maxt) %>%
  ggplot(aes(stcs, maxT)) +
  geom_point() +
  geom_path(aes(x, y), data=contour_95) +
  theme_bw() +
  geom_vline(xintercept = median(stcs), col = "red") +
  geom_hline(yintercept = median(maxT), col = "red") +
  geom_vline(xintercept = mean(stcs), col = "blue") +
  geom_hline(yintercept = mean(maxT), col = "blue") +
  geom_vline(xintercept = quantile(stcs, .95), col = "green") +
  geom_hline(yintercept = quantile(maxT, .95), col = "green")
  # geom_point(aes(x=tmp_stcs[length(tmp_stcs)],
  #                y=tmp_maxt[length(tmp_maxt)]),
  #            col="red")


#check if point is in the contour line
library(rgeos)
library(sp)


contour_polygon<- contour_95 %>% select(x, y) %>% Polygon
contour_polygon<- Polygons(list(contour_polygon), 1)
contour_polygon<- SpatialPolygons(list(contour_polygon))
proj4string(contour_polygon) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

plot(contour_polygon)

data = data.frame(f=99.9)
spdf = SpatialPolygonsDataFrame(contour_polygon,data)
summary(spdf)
spplot(spdf)
gContains(contour_95, c(2,3))




# one sim and its plot
sim_and_plot<- function(aaa=1){
  ind<- sample(x = length(perm_results$maxT), size = bootstrap_sample, replace = FALSE)
  tmp_stcs<- perm_results$stcs[ind]
  tmp_maxt<- perm_results$maxT[ind]

  fpr_stcs[i]<- tmp_stcs[length(tmp_stcs)] > quantile(tmp_stcs, probs = 1-alpha, names = FALSE)
  fpr_maxt[i]<- tmp_maxt[length(tmp_maxt)] > quantile(tmp_maxt, probs = 1-alpha, names = FALSE)

  # naive bivariate
  fpr_bivariate[i]<- fpr_stcs[i] | fpr_maxt[i] # | too liberal, & too conservative

  # using bivariate empirical cdf
  # biv_distr<- ebvcdf(tmp_stcs, tmp_maxt)
  # fpr_bivariate[i]<- biv_distr(tmp_stcs[length(tmp_stcs)], tmp_maxt[length(tmp_maxt)]) > (1-alpha) # Too conservative

  # Using empirical quantile contour line
  contour_line<- get_contour(tmp_stcs, tmp_maxt, alpha)
  fpr_bivariate[i]<- tmp_stcs[length(tmp_stcs)] > max(contour_line$x) | tmp_maxt[length(tmp_maxt)] > max(contour_line$y)

  data_frame(stcs=tmp_stcs, maxT=tmp_maxt) %>%
    ggplot(aes(stcs, maxT)) +
    geom_point() +
    geom_path(aes(x, y), data=contour_line) +
    theme_bw() +
    geom_vline(xintercept = median(stcs), col = "red") +
    geom_hline(yintercept = median(maxT), col = "red") +
    geom_vline(xintercept = mean(stcs), col = "blue") +
    geom_hline(yintercept = mean(maxT), col = "blue") +
    geom_vline(xintercept = quantile(stcs, .95), col = "green") +
    geom_hline(yintercept = quantile(maxT, .95), col = "green") +
    geom_point(aes(x=tmp_stcs[length(tmp_stcs)],
                   y=tmp_maxt[length(tmp_maxt)]),
               col="red") +
    geom_vline(xintercept = max(contour_line$x), col = "orange") +
    geom_hline(yintercept = max(contour_line$y), col = "orange")


}
sim_and_plot()



x <- 10*1:nrow(volcano)
y <- 10*1:ncol(volcano)
cl <- contourLines(x, y, volcano)
## summarize the sizes of each the contour lines :
cbind(lev = vapply(cl, `[[`, .5, "level"),
      n  = vapply(cl, function(l) length(l$x), 1))

z <- outer(-9:25, -9:25)
pretty(range(z), 10) # -300 -200 ... 600 700
utils::str(c2 <- contourLines(z))
# no segments for {-300, 700};
#  2 segments for {-200, -100, 0}
#  1 segment  for  100:600




#### kde copula package

# cannot load

### https://stackoverflow.com/questions/23437000/how-to-plot-a-contour-line-showing-where-95-of-values-fall-within-in-r-and-in
###  second answer

kk <- with(d,MASS::kde2d(x,y))
library(reshape2)
dimnames(kk$z) <- list(kk$x,kk$y)
dc <- melt(kk$z)

getLevel <- function(x,y,prob=0.95) {
  kk <- MASS::kde2d(x,y)
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}
L95 <- getLevel(d$x,d$y)
library(ggplot2); theme_set(theme_bw())
ggplot(d,aes(x,y)) +
  stat_density2d(geom="tile", aes(fill = ..density..),
                 contour = FALSE)+
  stat_density2d(colour="red",breaks=L95)

ggplot(dc,aes(x=Var1,y=Var2))+
   geom_tile(aes(fill=value))+
   geom_contour(aes(z=value),breaks=L95,colour="red")


# next answer


library(ggplot2)
library(MASS)
library(reshape2)
# create data:
set.seed(8675309)
Sigma <- matrix(c(0.1,0.3,0.3,4),2,2)
mv <- data.frame(mvrnorm(4000,c(1.5,16),Sigma))

# get the kde2d information:
mv.kde <- kde2d(d$x, d$y)
dx <- diff(mv.kde$x[1:2])  # lifted from emdbook::HPDregionplot()
dy <- diff(mv.kde$y[1:2])
sz <- sort(mv.kde$z)
c1 <- cumsum(sz) * dx * dy

# specify desired contour levels:
prob <- c(0.975,0.95,0.5, .6, .7, .8, .9)

# plot:
dimnames(mv.kde$z) <- list(mv.kde$x,mv.kde$y)
dc <- melt(mv.kde$z)
dc$prob <- approx(sz,1-c1,dc$value)$y
p <- ggplot(dc,aes(x=Var1,y=Var2))+
  geom_contour(aes(z=prob,color=..level..),breaks=prob)+
  geom_point(aes(x=x,y=y),data=d,alpha=0.1,size=1)
print(p)

### another
set.seed(1001)
d <- data.frame(x=rnorm(1000),y=rnorm(1000))
getLevel <- function(x,y,prob=0.95) {
  kk <- MASS::kde2d(x,y)
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}
L95 <- getLevel(d$x,d$y)
library(ggplot2); theme_set(theme_bw())
ggplot(d,aes(x,y)) +
  stat_density2d(geom="tile", aes(fill = ..density..),
                 contour = FALSE)+
  stat_density2d(colour="red",breaks=L95)

## using ks::kde
hscv1 <- Hscv(d)
fhat <- ks::kde(d, H=hscv1, compute.cont=TRUE)

dimnames(fhat[['estimate']]) <- list(fhat[["eval.points"]][[1]],
                                     fhat[["eval.points"]][[2]])
library(reshape2)
aa <- melt(fhat[['estimate']])

ggplot(aa, aes(x=Var1, y=Var2)) +
  geom_tile(aes(fill=value)) +
  geom_contour(aes(z=value), breaks=fhat[["cont"]]["50%"], color="red") +
  geom_contour(aes(z=value), breaks=fhat[["cont"]]["5%"], color="purple") +
  geom_point(data = d, aes(x, y), col = "red")

# wikipedia

library(ks)
data(faithful)
H <- Hpi(x=faithful)
fhat <- kde(x=faithful, H=H)

# my own ata
plot(fhat, display="filled.contour2")
points(d, cex=0.5, pch=16)
