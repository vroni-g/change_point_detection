#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************



# Maps and Frequency Barplots of breakpoint types
#**************************************************************************


#****************************************************************************
# Maps
#****************************************************************************

slopes_df <- readRDS("master_thesis/results/sig_mq_df.rds")

# make histogram of change point timing estimates
pdf_fn <- "MCUSUM_bp_timing_hist.pdf"
pdf(pdf_fn, height=8, width=11)
hist(slopes_df$bpl_year,breaks = 35, xlab = "Frequency", ylab = "Year of Change Point Estimate",
     main = "Global Frequencies of Change point Timing Estimates", cex.main = 2,
     cex.axis = 1.5, cex.lab = 1.5)
dev.off()

mat_bpt <- matrix(NA, nrow = 4320, ncol = 2160)
for(i in 1:nrow(slopes_df)){
  mat_bpt[slopes_df[i,"row"], slopes_df[i,"col"]] <- slopes_df[i,"changetype_int"]
}

library(stars)
library(tidyverse)
library(sp)
library(raster)

PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
sp <- st_as_stars(raster(matrix(1, nrow = 2160, ncol = 4320), xmn=-180, xmx=180, ymn=-90, ymx=90,
                         crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

sp_typ <- transmute(sp, cluster = as.vector(mat_bpt))
plot(sp_typ)

poly <- st_as_sf(x = sp_typ, as_points = FALSE, merge = TRUE, na.rm = TRUE) %>%
  st_make_valid
st_crs(poly) <- 4326
poly <- st_transform(poly, crs = st_crs(PROJ))
poly$cluster <- as.factor(poly$cluster)

# Additional data for plotting
#*****************************
masked <- readRDS("master_thesis/data_preprocessing/barren_land_ice_poly.RDS") %>%
  st_as_sf %>%
  st_transform(crs = st_crs(PROJ))
#plot(st_geometry(masked))

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = st_crs(PROJ))
load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
NE_box_rob <- spTransform(NE_box, CRSobj = PROJ)
rm(NE_box,NE_graticules, NE_places, NE_countries)

lons <- seq(-180, 180, by = 30)
lats <- seq(-60, 90, by = 15)
#*****************************

# WITHOUT majority vote
#***************************************************************************

# colors for change type maps and graphs:
# greening to browning (yellow or orange): #F7F121
# browning to greening (turquiose): #0FEAD2
# continued greening (dark green): #145607
# greening onset (middle green): #138211
# stalled greening (light green): #63CE47
# continued browning (dark blue):#181C67
# browning onset (middle blue): #1724FB
# stalled browning (light blue): #8DB5F1
# non non (gray): gray77

pdf_fn <- "MCUSUM_bp_types_map.pdf"
pdf(pdf_fn, height=8, width=11)
print(
  ggplot(data = world, size = 0.05) +
    geom_polygon(data=NE_box_rob, aes(x=long, y=lat), fill="aliceblue") +
    geom_sf(fill= "antiquewhite", size = 0.1) +
    xlab("") + ylab("") +
    labs(title = "Modified CUSUM Breakpoint Types") +
    geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
    theme(plot.title = element_text(hjust = 0.5, size=18),
          panel.background = element_rect(fill = "white"),
          plot.caption = element_text(size = 13)) +
    geom_sf(data = masked, fill = "white", color = NA) +
    geom_sf(data = poly, aes(fill = cluster), color = NA) +
    scale_fill_manual(values = c("#F7F121", "#0FEAD2","#145607", "#138211", "#63CE47",
                                 "#181C67","#1724FB", "#8DB5F1", "gray77"),
                      labels = c("Greening to Browning", "Browning to Greening",
                                 "Continued Greening", "Greening Onset", "Stalled Greening",
                                 "Continued Browning", "Browning Onset", "Stalled Browning",
                                 "Both Non-Significant")) +
    labs(fill='Breakpoint Types'))
dev.off()


# WITH majority vote
#***************************************************************************
res <- readRDS("master_thesis/results/BU_MCUSUM_cluster_original.rds")
clus <- res[[2]]

inds <- readRDS("master_thesis/results/sig_nperm30_combined_stcsmq.rds")[[2]]
sig_mat <- readRDS("master_thesis/results/sig_nperm30_combined_stcsmq.rds")[[1]]

get_types <- function(i,mat,types){
  mat[mat %in% i] <- -10
  mat[mat != -10] <- 0
  mat <- dplyr::na_if(mat, 0)
  mat[!is.na(mat)] <- TRUE

  l <- types*mat
  l <- as.vector(l)
  l <- l[!is.na(l)]
  return(as.vector(l))
}

clust_mat <- clus$clusters
sig_clust <- clust_mat*sig_mat
types_clust <- lapply(inds, get_types, mat = sig_clust, types = mat_bpt)

types_maj <- lapply(types_clust, function(x) names(sort(-table(x)))[1]) %>%
  unlist
table(types_maj)

sig_maj <- sig_clust

j = 1
for(i in inds){
  sig_maj[sig_maj  == i] <- types_maj[j]
  j = j+1
}
sig_maj <- type.convert(sig_maj, as.is = FALSE)

maj_types_df <- c()
for(i in 1:nrow(slopes_df)){
  maj_types_df <- c(maj_types_df, sig_maj[slopes_df[i, "row"], slopes_df[i, "col"]])
}

slopes_df$ctype_major <- maj_types_df

sp_typ2 <- transmute(sp, cluster = as.vector(sig_maj))
plot(sp_typ2)

poly2 <- st_as_sf(x = sp_typ2, as_points = FALSE, merge = TRUE, na.rm = TRUE) %>%
  st_make_valid
st_crs(poly2) <- 4326
poly2 <- st_transform(poly2, crs = st_crs(PROJ))
poly2$cluster <- as.factor(poly2$cluster)


pdf_fn <- "/home/veronika/CPD/MCUSUM/MCUSUM_bp_types_majority.pdf"
pdf(pdf_fn, height=8, width=11)
print(
  ggplot(data = world, size = 0.05) +
    geom_polygon(data=NE_box_rob, aes(x=long, y=lat), fill="aliceblue") +
    geom_sf(fill= "antiquewhite", size = 0.1) +
    xlab("") + ylab("") +
    labs(title = "Modified CUSUM Breakpoint Types") +
    geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
    theme(plot.title = element_text(hjust = 0.5, size=18),
          panel.background = element_rect(fill = "white"),
          plot.caption = element_text(size = 13)) +
    geom_sf(data = masked, fill = "white", color = NA) +
    geom_sf(data = poly2, aes(fill = cluster), color = NA) +
    scale_fill_manual(values = c("#F7F121", "#0FEAD2","#145607", "#138211", "#63CE47",
                                 #"#181C67",
                                 "#1724FB", "#8DB5F1", "gray77"),
                      labels = c("Greening to Browning", "Browning to Greening",
                                 "Continued Greening", "Greening Onset", "Stalled Greening",
                                 #"Continued Browning",
                                 "Browning Onset", "Stalled Browning",
                                 "Both Non-Significant")) +
    labs(fill='Breakpoint Types'))
dev.off()


#****************************************************************************
# Frequency Barplots
#****************************************************************************

# without majority vote
#************************************************
types <- slopes_df %>% count(changetype)

crls <- c("#145607", "#138211","#63CE47", "#F7F121", "#181C67", "#1724FB",
          "#8DB5F1", "#0FEAD2","gray77")

nam <- c("Continued Greening", "Greening Onset", "Stalled Greening", "Greening to Browning",
         "Continued Browning", "Browning Onset", "Stalled Browning", "Browning to Greening",
         "Both Non-Significant")

type_order <- c("green_green", "non_green", "green_non", "green_brown", "brown_brown",
                "non_brown", "brown_non", "brown_green", "non_non")
frac <- lapply(type_order, function(x) dplyr::filter(types, changetype == x)[[2]]/sum(types$n)) %>% unlist

# with majority vote
#************************************************
types2 <- slopes_df %>% count(ctype_major)
type_order2 <- c(3,4,5,1,7,8,2,9)
frac2 <- lapply(type_order2, function(x) dplyr::filter(types2, ctype_major == x)[[2]]/sum(types2$n)) %>% unlist

crls2 <- c("#145607", "#138211","#63CE47", "#F7F121", "#1724FB",
           "#8DB5F1", "#0FEAD2","gray77")

nam2 <- c("Continued Greening", "Greening Onset", "Stalled Greening", "Greening to Browning",
          "Browning Onset", "Stalled Browning", "Browning to Greening",
          "Both Non-Significant")


rotate_x <- function(data, labels_vec, rot_angle, tit, corls) {
  plt <- barplot(data, ylab = "% of pixels", xaxt="n", main = tit, col = corls, cex.lab = 1.2, cex.main = 1.7)
  text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE)#, cex=0.9)
}


pdf_fn <- "MCUSUM_bp_types_hists.pdf"
pdf(pdf_fn, height=11, width=12)
par(mfrow = c(2,1))
rotate_x(frac, nam, 30, "Change Point Types without Cluster Majority Vote",corls = crls)
rotate_x(frac2, nam2, 30, "Change Point Types with Cluster Majority Vote",corls = crls2)
dev.off()


