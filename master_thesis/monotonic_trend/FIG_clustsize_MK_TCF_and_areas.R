#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************



# Check clustersizes of significant clusters of MK Trend Test with TCF correction
#************************************************************
library(tidyverse)

#res <- readRDS("master_thesis/results/BU_LAI_MK_nperm_1000_al5.rds")
res <- readRDS("master_thesis/results/BU_LAI_MK_nperm_1000_al10.rds")

dataset <- "BU"
alpha_local <- 0.1
alpha_global <- 0.1

clust_size <- res$original_cluster$cluster.count
length(res$original_cluster$cluster.count)
wt_vals <- res$original_wt
res_df <- data.frame(clust_size, wt_vals)
names(res_df) <- c("cluster_size","wt_values")
sum(res_df$cluster_size)

# values declared significant by tippet
tipp_9 <- res_df %>%
  dplyr::filter(wt_vals >= quantile(res$wt, probs = .9))
sum(tipp_9$cluster_size)
# values declared significant by stcs
stcs_9 <- res_df %>%
  dplyr::filter(cluster_size >= quantile(res$stcs, probs = .9))
sum(stcs_9$cluster_size)

# values declared significant by tippet but not by stcs
tip_small_9 <- res_df %>%
  dplyr::filter(wt_vals >= quantile(res$wt, probs = .9) & cluster_size < quantile(res$stcs, probs = .95))


# ********* PLOT WITH BASE R ******************
# ********************************************
#**** Alpha global 10 % *****

pdf_fn <- paste0("TCF_MK_plot_",dataset,
                 "_MK_1981-2018_al",alpha_local*100,"_ag", alpha_global*100, ".pdf")
pdf(pdf_fn, height=7, width=8)
options(scipen=999)


plot(res_df$cluster_size, res_df$wt_values, main="Suprathreshold Cluster of MK Trend Test (local alpha = 0.1)",
     xlab = 'Suprathreshold Cluster Size',
     ylab = 'Tippet combining function value')
abline(h=quantile(res$wt, probs = .95),col = "red", lty = 3)
abline(h=quantile(res$wt, probs = .9),col = "red")
abline(v=quantile(res$stcs, probs = .95),col = "red", lty = 3)
abline(v=quantile(res$stcs, probs = .9),col = "red")
text(x=42000,y=quantile(res$wt, probs = .95)+0.15,pos=4,label = "global alpha = 0.05", col = "red",cex = 0.77)
text(x=42000,y=quantile(res$wt, probs = .9)+0.15,pos=4,label = "global alpha = 0.1", col = "red",cex = 0.77)
text(x = 35000, y=2.9, pos=4, label= "For global alpha = 0.1 :",cex = 0.7, font = 2)
text(x=35000, y=2.3, label= paste0("Significant Cluster STCS:  ", nrow(stcs_9),
                                   "\nSignificant Cluster TCF:  ", nrow(tipp_9),
                                   "\nSignificant Cluster TCF but not STCS: ", nrow(tip_small_9)),
     cex = 0.7, pos=4)
text(x=35000, y=1.75, label= paste0("Total number of Cluster: ", nrow(res_df)),
     cex = 0.7, font = 3, pos=4)
rect(34500, 1.5, 55500, 3.2, border = "black")
dev.off()
#***********************************



# AREAS of significant cluster for different correction methods:
#*********************************************************************
library(raster)
cluster <- res$original_cluster$clusters
cluster[cluster == 0] <- NA
cluster[!is.na(cluster)] <- 1

ras <- raster(t(cluster), xmn=-180, xmx=180, ymn=-90, ymx=90,
              crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
plot(ras)

# Area of uncorrected results
#****************
area_uncorr <- tapply(area(ras,na.rm = T), ras[], sum)
#****************

inds_wt_non <- which(res$original_wt <= quantile(res$wt, probs = 1-alpha_global))
inds_stcs_non <- which(res$original_cluster$cluster.count <= quantile(res$stcs, probs = 1-alpha_global))

cluster <- res$original_cluster$clusters
cluster[cluster == 0] <- NA
cluster[cluster %in% inds_stcs_non] <- NA
cluster[!is.na(cluster)] <- 1
ras <- raster(t(cluster), xmn=-180, xmx=180, ymn=-90, ymx=90,
              crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
plot(ras)

# Area of STCS results
#****************
area_stcs <- tapply(area(ras,na.rm = T), ras[], sum)
#****************

cluster <- res$original_cluster$clusters
cluster[cluster == 0] <- NA
cluster[cluster %in% inds_wt_non] <- NA
cluster[!is.na(cluster)] <- 1
ras <- raster(t(cluster), xmn=-180, xmx=180, ymn=-90, ymx=90,
              crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
plot(ras)

# Area of TCF results
#****************
area_tip <- tapply(area(ras,na.rm = T), ras[], sum)
#****************

# Greening and Browning Cluster
#*********************************

gb <- res[["original_stat"]]*cluster
gb[gb >= 0] <- 10
gb <- dplyr::na_if(gb, 0)
gb[!is.na(gb)] <- 1
table(gb)

inds_wt <- which(res$original_wt > quantile(res$wt, probs = 1-alpha_global))

inds_brown <- c()
for(i in inds_wt){
  cluster <- res$original_cluster$clusters
  cluster[cluster != i] <- 0
  cluster <- dplyr::na_if(cluster, 0)
  stat <- res[["original_stat"]]*cluster
  m <- mean(as.vector(stat), na.rm = T)
  cat("cluster ",i, " has mean ", m, "\n")
  if(m<0) inds_brown <- c(inds_brown, i)
}

sum(res$original_cluster$cluster.count[inds_brown])

cluster <- res$original_cluster$clusters
cluster[cluster %in% inds_brown] <- -100
cluster[cluster != -100] <- 0
cluster <- dplyr::na_if(cluster, 0)
table(cluster)

ras <- raster(t(cluster), xmn=-180, xmx=180, ymn=-90, ymx=90,
              crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
plot(ras)

# Area of TCF browning cluster
#****************
area_brown <- tapply(area(ras,na.rm = T), ras[], sum)



