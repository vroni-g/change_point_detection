#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# !! Input Data is not available in GitHub Repository !!

# Script to create a map of quality issues of BU LAI data
#************************************************************

library(dplyr)
library(stars)
library(sp)

qf <- readRDS("/home/veronika/CPD/data/BU_1981_2018/shareQF_complSeries_BU.rds")
sp <- read_stars("/home/veronika/CPD/data/spatial_BU.tif")

sp_qf <- transmute(sp, cluster = as.vector(qf))
plot(sp_qf)

qf_poly <- st_as_sf(x = sp_qf, as_points = FALSE, merge = TRUE, na.rm = TRUE) %>%
  st_make_valid
st_crs(qf_poly) <- 4326

PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
qf_poly_rob <- st_transform(qf_poly, crs = st_crs(PROJ))
rm(qf_poly)


# Additional data for plotting
#******************************
library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = st_crs(PROJ))
load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
NE_box_rob <- spTransform(NE_box, CRSobj = PROJ)
rm(NE_box,NE_graticules, NE_places, NE_countries)

lons <- seq(-180, 180, by = 30)
lats <- seq(-60, 90, by = 15)
#******************************

# save as pdf
pdf_fn <- paste0("QF_map.pdf")
pdf(pdf_fn, height=5, width=8)

print(
qf_plot <-  ggplot(data = world, size = 0.05) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), fill="aliceblue") +
  geom_sf(fill= "grey86", size = 0.1) +
  xlab("") + ylab("") +
  labs(title = "% of Time Points with Quality Issues") +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = qf_poly_rob, aes(fill = cluster), colour = NA) +
  scale_fill_gradientn(colors = c("#f8f2f1","#f8e2e0","#F1B399","#E6806F","#DA5333","#B72310","#80150C"), na.value = "transparent") +
  guides(fill = guide_colourbar(title = "Percentage %")) +
  theme(plot.title = element_text(size=10))
)

dev.off()
