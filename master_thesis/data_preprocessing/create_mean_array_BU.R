#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# !! Input Data is not available in GitHub Repository !!

# Create array of all yearly means of BU GIMMS data
#************************

library(tidyverse)
# directory to store final array
dir <- "/home/veronika/CPD/data/"
# yearly mean files produced by José Cortés, stored on HPC facilities of the 
# GIScience group, Institute of Geography, FSU Jena
dir <- "/home/jose/LAI/data/CHEN_RANGA/AVHRR/yearly_mean"
files <-list.files(dir, pattern = "\\.nc$", full.names = T)

arr <- array(NA,dim = c(4320, 2160,length(files)))

for(i in 1:length(files)){
  load(files[i])
  data[is.nan(data)] <- NA
  arr[,,i] <- data
}

saveRDS(arr, paste0(dir,"BU_LAI_yearlymean_1981_2018"))
arr82 <- arr[-1]
saveRDS(arr82, paste0(dir,"BU_LAI_yearlymean_1982_2018"))

arr <- arr*10000
mode(arr) <- "integer"
saveRDS(arr,  paste0(dir,"10000int_BU_LAI_yearlymean_1981_2018"))
arr82 <- arr[-1]
saveRDS(arr,  paste0(dir,"10000int_BU_LAI_yearlymean_1982_2018"))










