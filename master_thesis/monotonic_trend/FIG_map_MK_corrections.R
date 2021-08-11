#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************


# Maps for tippet combining function Mann Kendall results
#****************************************************
library(stars)
library(raster)
library(tidyverse)
library(sp)

dataset <- "BU"
alpha_local <- 0.1
alpha_global <- 0.1

res <- readRDS("master_thesis/results/BU_LAI_MK_nperm_1000_al10.rds")
cluster_mat<- res$original_cluster[["clusters"]] %>% 
  as.matrix
cluster_mat[cluster_mat== 0] <- NA
stat_mat <- res$original_stat

inds_wt_non <- which(res$original_wt <= quantile(res$wt, probs = 1-alpha_global))
inds_stcs_non <- which(res$original_cluster$cluster.count <= quantile(res$stcs, probs = 1-alpha_global))

sp <- st_as_stars(raster(matrix(1, nrow = 2160, ncol = 4320), xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))  

PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Uncorrected results
#****************
uc <- cluster_mat
uc[!is.na(uc)] = TRUE
uc_sig <- sign(stat_mat)*uc
sp_uc <- transmute(sp, cluster = as.vector(uc_sig))
#plot(sp_uc)

uc_poly <- st_as_sf(x = sp_uc, as_points = FALSE, merge = TRUE, na.rm = TRUE) %>%
  st_make_valid 
st_crs(uc_poly) <- 4326
uc_poly$cluster <- as.factor(uc_poly$cluster)
uc_poly <- st_transform(uc_poly, crs = st_crs(PROJ))

# TCF results
#****************
tip <- cluster_mat
tip[tip %in% inds_wt_non] <- NA
tip[!is.na(tip)] = TRUE
tip_sig <- sign(stat_mat)*tip
sp_tip <- transmute(sp, cluster = as.vector(tip_sig))
#plot(sp_tip)

tip_poly <- st_as_sf(x = sp_tip, as_points = FALSE, merge = TRUE, na.rm = TRUE) %>%
  st_make_valid 
st_crs(tip_poly) <- 4326
tip_poly$cluster <- as.factor(tip_poly$cluster)
tip_poly <- st_transform(tip_poly, crs = st_crs(PROJ))

# STCS results
#****************
stcs <- cluster_mat
stcs[stcs %in% inds_stcs_non] <- NA
stcs[!is.na(stcs)] = TRUE

stcs_sig <- sign(stat_mat)*stcs
sp_stcs <- transmute(sp, cluster = as.vector(stcs_sig))
#plot(sp_stcs)

stcs_poly <- st_as_sf(x = sp_stcs, as_points = FALSE, merge = TRUE, na.rm = TRUE) %>%
  st_make_valid 
st_crs(stcs_poly) <-4326
stcs_poly$cluster <- as.factor(stcs_poly$cluster)
stcs_poly <- st_transform(stcs_poly,crs = st_crs(PROJ))


# Additional data for plotting
#**********************************************
masked <- readRDS("master_thesis/data_preprocessing/barren_land_ice_poly.RDS") %>%
  st_as_sf %>%
  st_transform(crs = st_crs(PROJ))

library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = st_crs(PROJ))
load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
NE_box_rob <- spTransform(NE_box, CRSobj = PROJ)
rm(NE_box,NE_graticules, NE_places, NE_countries)

lons <- seq(-180, 180, by = 30)
lats <- seq(-60, 90, by = 15)


# PLOT MULTIPLE
#************************************
stcs_plot <- ggplot(data = world, size = 0.05) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), fill="aliceblue") +
  geom_sf(fill= "antiquewhite", size = 0.1) +
  xlab("") + ylab("") +
  labs(subtitle = "Suprathreshold Cluster Size (STCS)",
       title = "Mann Kendall Trend Test") +
  theme(plot.title = element_text(hjust = 0.5,size=13),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        panel.background = element_rect(fill = "white")) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = masked, fill = "white", color = NA) +
  geom_sf(data = stcs_poly, aes(fill = cluster), colour = NA) +
  scale_fill_manual(values=c("-1" = "#5F4220","1" = "#47AB0F"), guide = F)

tip_plot <- ggplot(data = world, size = 0.05) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), fill="aliceblue") +
  geom_sf(fill= "antiquewhite", size = 0.1) +
  xlab("") + ylab("") +
  labs(title = " Tippet Combining Function (TCF)") +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = masked, fill = "white", color = NA) +
  geom_sf(data = tip_poly, aes(fill = cluster), colour = NA) +
  scale_fill_manual(values=c("-1" = "#5F4220","1" = "#47AB0F"), guide = F)

uc_plot <- ggplot(data = world, size = 0.05) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), fill="aliceblue") +
  geom_sf(fill= "antiquewhite", size = 0.1) +
  xlab("") + ylab("") +
  labs(title = " Uncorrected Results",
       caption = paste0("Alpha local = ",alpha_local, ", Alpha global = ", alpha_global )) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        panel.background = element_rect(fill = "white"),
        plot.caption = element_text(size = 6)) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
  geom_sf(data = masked, fill = "white", color = NA) +
  geom_sf(data = uc_poly, aes(fill = cluster), colour = NA) +
  scale_fill_manual(values=c("-1" = "#5F4220","1" = "#47AB0F"), guide = F)

pdf_fn <- paste0("TCF_MK_map_",dataset,
                    "_MK_1981-2018_al",alpha_local*100,"_ag", alpha_global*100, ".pdf")
pdf(pdf_fn, height=7, width=4)
print(ggpubr::ggarrange(stcs_plot, tip_plot, uc_plot, nrow = 3))
dev.off()


