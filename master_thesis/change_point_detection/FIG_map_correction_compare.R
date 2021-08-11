#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************



# MAP of uncorrected, stcs and manual combination correction
#***************************************************************************

library(stars)
library(sp)
library(raster)

PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
clus <- readRDS("master_thesis/results/BU_MCUSUM_cluster_original.rds")[[2]]

sp <- st_as_stars(raster(matrix(1, nrow = 2160, ncol = 4320), xmn=-180, xmx=180, ymn=-90, ymx=90, 
                         crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))  
clust_mat <- clus[[1]]

# uncorrected
#*****************************************
uc <- clust_mat
uc <- dplyr::na_if(uc, 0)
uc[!is.na(uc)] <- 1
table(uc)
sum(clus$cluster.count)

# get area
#****************
ras <- raster(t(uc), xmn=-180, xmx=180, ymn=-90, ymx=90,
              crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
plot(ras)
area_uc <- tapply(area(ras,na.rm = T), ras[], sum)
#****************

sp_uc <- transmute(sp, cluster = as.vector(uc))
uc_poly <- st_as_sf(x = sp_uc, as_points = FALSE, merge = TRUE, na.rm = TRUE) %>%
  st_make_valid
st_crs(uc_poly) <- 4326
uc_poly$cluster <- as.factor(uc_poly$cluster)
uc_poly <- st_transform(uc_poly, crs = st_crs(PROJ))

# stcs corrected
#********************************************

st <- clust_mat
st[st %in% inds_stcs] <- -10
st[st!=-10] <- 0
st <- dplyr::na_if(st, 0)
st[!is.na(st)] <- 1
table(st)
sum(clus$cluster.count[inds_stcs])

# get area
#****************
ras <- raster(t(st), xmn=-180, xmx=180, ymn=-90, ymx=90,
              crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
plot(ras)
area_st <- tapply(area(ras,na.rm = T), ras[], sum)
#****************

sp_st <- transmute(sp, cluster = as.vector(st))
st_poly <- st_as_sf(x = sp_st, as_points = FALSE, merge = TRUE, na.rm = TRUE) %>%
  st_make_valid
st_crs(st_poly) <- 4326
st_poly$cluster <- as.factor(st_poly$cluster)
st_poly <- st_transform(st_poly, crs = st_crs(PROJ))

# corrected with stcs and manual quantiles
#**********************************************************

mq <- readRDS("master_thesis/results/sig_nperm30_combined_stcsmq.rds")[[1]]

# get area
#****************
ras <- raster(t(mq), xmn=-180, xmx=180, ymn=-90, ymx=90,
              crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
plot(ras)
area_mq <- tapply(area(ras,na.rm = T), ras[], sum)
#****************

sp_mq <- transmute(sp, cluster = as.vector(mq))
mq_poly <- st_as_sf(x = sp_mq, as_points = FALSE, merge = TRUE, na.rm = TRUE) %>%
  st_make_valid
st_crs(mq_poly) <- 4326
mq_poly$cluster <- as.factor(mq_poly$cluster)
mq_poly <- st_transform(mq_poly, crs = st_crs(PROJ))


#********************************* PLOT *************************************

# Additional data for plotting
#******************************
masked <- readRDS("master_thesis/data_preprocessing/data/barren_land_ice_poly.RDS") %>%
  st_as_sf %>%
  st_transform(crs = st_crs(PROJ))

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = st_crs(PROJ))
load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
NE_box_rob <- spTransform(NE_box, CRSobj = PROJ)
rm(NE_box,NE_graticules, NE_places, NE_countries)

lons <- seq(-180, 180, by = 30)
lats <- seq(-60, 90, by = 15)

# PLOT MULTIPLE
#************************************
alpha_global <- 0.1
alpha_local <- 0.1

stcs_plot <- ggplot(data = world, size = 0.05) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), fill="aliceblue") +
  geom_sf(fill= "antiquewhite", size = 0.1) +
  xlab("") + ylab("") +
  labs(subtitle = "Suprathreshold Cluster Size (STCS)",
       title = "Modified CUSUM") +
  theme(plot.title = element_text(hjust = 0.5,size=13),
        plot.subtitle = element_text(hjust = 0.5, size=8),
        panel.background = element_rect(fill = "white")) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
  geom_sf(data = masked, fill = "white", color = NA) +
  geom_sf(data = uc_poly, fill = "#1E88E5", colour = NA) +
  geom_sf(data = st_poly, fill = "#FFC107", colour = NA)

mq_plot <- ggplot(data = world, size = 0.05) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), fill="aliceblue") +
  geom_sf(fill= "antiquewhite", size = 0.1) +
  xlab("") + ylab("") +
  labs(title = "Manual Combination",
       caption = paste0("Alpha local = ",alpha_local, ", Alpha global = ", alpha_global )) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
  theme(panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size=8)) +
  geom_sf(data = masked, fill = "white", color = NA) +
  geom_sf(data = uc_poly, fill = "#1E88E5", colour = NA) +
  geom_sf(data = mq_poly, fill = "#FFC107", colour = NA)

pdf_fn <- paste0("MCUSUM_BU_map_stcs_mq.pdf")
pdf(pdf_fn, height=6, width=5.2)
print(ggpubr::ggarrange(stcs_plot, mq_plot, nrow = 2))
dev.off()





# cluster sizes of "discarded clusters" vs. significant cluster
#*****************************************************************
inds_adj <-  readRDS("master_thesis/results/sig_nperm30_combined_stcsmq.rds")[[2]]
inds_rem <- c(1:length(or[[2]]$cluster.count))[-inds_adj]
rem_counts <- or[[2]]$cluster.count[inds_rem]
sig_counts <- or[[2]]$cluster.count[inds_adj]
sig_counts <- sig_counts[sig_counts < thr]
hist(rem_counts, breaks = 40)
hist(sig_counts, breaks = 40)

bins <- c(1,2,6,11,21,51,101, 201, 301, 401, 501, 601, 701, 801, 909, 80000)
tags <- c("1", "2-5", "6-10", "11-20", "21-50", "51-100", "101-200", "201-300",
          "301-400", "401-500", "501-600", "601-700", "701-800", "801-900", ">900")

df_sig<- tibble::tibble(sig_counts)
group_tags_s <- cut(df_sig$sig_counts,
                    breaks=bins,
                    include.lowest=TRUE,
                    right=FALSE,
                    labels=tags)
df_sig$bins <- group_tags_s

df_rem<- tibble::tibble(rem_counts)
group_tags_r <- cut(df_rem$rem_counts,
                    breaks=bins,
                    include.lowest=TRUE,
                    right=FALSE,
                    labels=tags)
df_rem$bins <- group_tags_r

cnts_sig <- df_sig %>% count(bins)
cnts_rem <- df_rem %>% count(bins)

pdf_fn <- paste0("/home/veronika/CPD/MCUSUM/MCUSUM_clustsize_compare_sigrem.pdf")
pdf(pdf_fn, height=14, width=13)
par(mfrow=c(2,1))
barplot(height = cnts_sig$n, names = cnts_sig$bins, xlab = "Cluster size groups",
        ylab = "Count of significant cluster", ylim = c(0,120),
        main = "Size distribution of corrected significant clusters",
        cex.lab=1.5, cex.main=1.8)
barplot(height = cnts_rem$n, names = cnts_rem$bins, xlab = "Cluster size groups",
        ylab = "Count of non-significant cluster", ylim = c(0,140),
        main = "Size distribution of uncorrected clusters",
        cex.lab=1.5, cex.main=1.8)
dev.off()
