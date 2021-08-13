#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************


# Script to derive percentages of annual NA values per pixel
#***************************************

na_perc <- function(x){
  sum(is.na(x))/length(x)
}

data_lai <- readRDS("/home/veronika/CPD/data/10000int_BU_LAI_yearlymean_1981_2018")

perc <- apply(data_lai, 1:2, na_perc)
tp <- table(perc)
(sum(perc>=0.9)-sum(perc==1))/sum(perc >= 0.9) # amount of data between 90% and 100% data availability
(sum(perc>0)-sum(perc>=0.9))/sum(perc>0) # amount of data between >0% and 90% data availability in relation to all data with at least one data point

# the vast majority of those points (min. 90% data) were a complete record, 
# only 38 out of 7578132 had more than 90% and less then 100% values
# less than 1% of the data points with any values had less than 90% values
