#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************


# overlay / compare change point detection with MK Trend Test
#****************************************************************
#*
library(tidyverse)
alpha_global <- 0.1

sig_df <- readRDS("master_thesis/results/sig_mq_df.rds")

res <- readRDS("master_thesis/results/BU_LAI_tippet_nperm_1000_al10.rds")
inds_wt <- which(res$original_wt >= quantile(res$wt, probs = 1-alpha_global))

cluster <- res$original_cluster$clusters
cluster[cluster %in% inds_wt] <- -10
cluster[cluster != -10] <- 0
cluster <- dplyr::na_if(cluster, 0)
cluster[cluster== -10] <- 1

inds_over <- rep(FALSE, nrow(sig_df))

for(i in 1:nrow(sig_df)){
  val <- cluster[[sig_df[i, "row"], sig_df[i, "col"]]]
  if(is.na(val)) val <- -999
  if(val == 1) inds_over[i] <- TRUE
}
sum(inds_over)
sum(inds_over)/nrow(sig_df)
sig_df$MK_sig <- inds_over

saveRDS(sig_df, "master_thesis/results/sig_mq_df.rds")
sig_df <- readRDS("master_thesis/results/sig_mq_df.rds")

# BARPLOTS
#*************************************************************************
types_MK <- filter(sig_df, MK_sig) %>% count(changetype)

crls <- c("#145607", "#138211","#63CE47", "#F7F121", "#1724FB",
          "#8DB5F1", "#0FEAD2","gray77")

nam <- c("Continued Greening", "Greening Onset", "Stalled Greening", "Greening to Browning",
          "Browning Onset", "Stalled Browning", "Browning to Greening",
         "Both Non-Significant")

type_order <- c("green_green", "non_green", "green_non", "green_brown",
                "non_brown", "brown_non", "brown_green", "non_non")
freq <- lapply(type_order, function(x) dplyr::filter(types_MK, changetype == x)[[2]]) %>% unlist



types_MK_non <- filter(sig_df, !MK_sig) %>% count(changetype)
crls2 <- c("#145607", "#138211","#63CE47", "#F7F121", "#181C67", "#1724FB",
          "#8DB5F1", "#0FEAD2","gray77")

nam2 <- c("Continued Greening", "Greening Onset", "Stalled Greening", "Greening to Browning",
         "Continued Browning", "Browning Onset", "Stalled Browning", "Browning to Greening",
         "Both Non-Significant")

type_order2 <- c("green_green", "non_green", "green_non", "green_brown", "brown_brown",
                "non_brown", "brown_non", "brown_green", "non_non")
freq2 <- lapply(type_order2, function(x) dplyr::filter(types_MK_non, changetype == x)[[2]]) %>% unlist


rotate_x <- function(data, labels_vec, rot_angle, tit, corls) {
  plt <- barplot(data, ylab = "Count of pixels", xaxt="n", main = tit, col = corls, cex.lab = 1.2, cex.main = 1.7)
  text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE)#, cex=0.9)
}

pdf_fn <- "MCUSUM_bp_types_MK_overlay.pdf"
pdf(pdf_fn, height=11, width=12)
par(mfrow = c(2,1))
rotate_x(frac, nam, 30, "Change Point Types of Pixels with significant MK Test",corls = crls)
rotate_x(frac2, nam2, 30, "Change Point Types of Pixels with non-significant MK Test",corls = crls2)
dev.off()
#*************************************************************************



# Map of MCUSUM pixel in two colors for non/significant MK Tests
#****************************************************************

mat <- matrix(NA, nrow = dim(res$original_stat)[1], ncol = dim(res$original_stat)[2])
for(i in 1:nrow(sig_df)){
  if(sig_df[i,"MK_sig"]){
    mat[sig_df[i,"row"], sig_df[i,"col"]] <- 1
  } else {
    mat[sig_df[i,"row"], sig_df[i,"col"]] <- 2
  }
}
library(stars)
library(sp)

PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
sp <- read_stars("/home/veronika/CPD/data/spatial_BU.tif")

sp_typ <- transmute(sp, mk = as.vector(mat))
plot(sp_typ)

poly <- st_as_sf(x = sp_typ, as_points = FALSE, merge = TRUE, na.rm = TRUE) %>%
  st_make_valid
st_crs(poly) <- 4326
poly <- st_transform(poly, crs = st_crs(PROJ))
poly$mk <- as.factor(poly$mk)

# Additional data for plotting
#*****************************
masked <- readRDS("/home/veronika/CPD/data/barren_land_ice_poly.RDS") %>%
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


pdf_fn <- "/home/veronika/CPD/MCUSUM/MCUSUM_bp_typ_map_MK.pdf"
pdf(pdf_fn, height=8, width=11)
print(
  ggplot(data = world, size = 0.05) +
    geom_polygon(data=NE_box_rob, aes(x=long, y=lat), fill="aliceblue") +
    geom_sf(fill= "antiquewhite", size = 0.1) +
    xlab("") + ylab("") +
    labs(title = "Change Point Detection and Monotonic Trend") +
    geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
    theme(plot.title = element_text(hjust = 0.5, size=18),
          panel.background = element_rect(fill = "white"),
          plot.caption = element_text(size = 13)) +
    geom_sf(data = masked, fill = "white", color = NA) +
    geom_sf(data = poly, aes(fill = mk), color = NA) +
    scale_fill_manual(values = c("#056F1F", "#1E88E5"),
                      labels = c("Significant Monotonic Trend", "No Significant Monotonic Trend")) +
    labs(fill=''))
dev.off()
