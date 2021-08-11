#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# !! Input Data is not available in GitHub Repository !!


# Script to combine quality flag values for all years
#********************************************************


# YEARLY Quality flag (QF) MATRICES
#***************************************
# read binary files from Chi Chen and Ranga, and create yearly rds files for quality flags
#***************************************

create_yearly_nc_files_chen_ranga<- function(data_path, yearly_aggregating_fn, ..., res_dir){
  # data_path: location of binary files
  # yearly_aggregating_fn: how to summarize yearly data
  # ... extra arguments for yearly_aggregating_fn
  # res_dir: directory within data_path where results are saved. 
  #         Created automatically if it does not exist.
  
  NCORES<- 5
  library(future.apply)
  
  plan(multisession, workers = NCORES)
  options(future.globals.maxSize= +Inf)
  on.exit(future:::ClusterRegistry("stop"))
  
  
  files<- list.files(data_path, "\\.bin$")
  
  # ... are extra inputs for aggregating fn, e.g., the quantile
  # na.rm = TRUE is assumed
  
  # years and dim for each data product
  product<- tail(unlist(strsplit(data_path, split = "/")), n = 1)
  if(product == "AVHRR") {
    years<- 1981:2018
    data_dim<- c(4320, 2160)
    data_na<- -32768
    files<- files[grep("Percentile", files)] # the other option is percentile
  } else if(product == "MODIS") {
    years<- 2000:2019
    data_dim<- c(3600, 7200)
    data_na<- -32768
  }
  
  print(paste(product, range(years)[1], range(years)[2], data_dim[1], 
              data_dim[2], sep = ", "))
  
  create_dir<- !dir.exists(res_dir)
  if(create_dir) dir.create(res_dir)
  
  create_parallel<- function(year, ...){
    files_in_year<- files[grep(year, files)]
    data<- array(NA, dim = c(length(files_in_year), data_dim[1]*data_dim[2]))
    cat("Files for year ", year, ": ", files_in_year, "\n")# check years are good
    
    for(i in 1:length(files_in_year)){
      tmp<- readBin(paste0(data_path, "/", files_in_year[i]), what = "integer", 
                    n = (data_dim[1]*data_dim[2]), size = 2, endian = "big") 
      tmp[tmp == data_na]<- NA
      tmp <- floor(tmp / 2e3)
      tmp[tmp == 2] <- 1
      data[i,]<- tmp
    }
    
    data<- apply(data, 2, yearly_aggregating_fn, ...)
    data<- matrix(data, nrow = data_dim[1], ncol = data_dim[2])
    
    saveRDS(data, file = paste0(res_dir, "/SumQF.LAI.", product, ".", data_dim[1], ".", 
                             data_dim[2], ".", year, ".rds"))
    return()
  }
  future_lapply(years[1]:range(years)[2], create_parallel, ...)
  future:::ClusterRegistry("stop")
}


sum_fn<- function(x) sum(x)
# path were quality flag data is stored
dpath <- "/home/jose/LAI/data/CHEN_RANGA/AVHRR"
# target directory to save yearly values
res_dir <- "/home/veronika/CPD/data/BU_1981_2018/sumQF_yearly"
create_yearly_nc_files_chen_ranga(dpath, 
                                  yearly_aggregating_fn = sum_fn, # values are summed up for each year
                                  res_dir = res_dir)


# COMPLETE QF MATRIX
#***************************************
# read yearly nc files for quality flags and sum up
#***************************************

files <- list.files(dpath,"\\.bin$")
count <- files[grep("Percentile", files)]

files <- list.files(res_dir, "\\.rds$", full.names = T)

arr <- array(NA,dim = c(4320, 2160,length(files)))

for(i in 1:length(files)){
  arr[,,i] <- readRDS(files[i])
}

sum(arr[1200, 500,])
sums <- apply(arr, 1:2, sum)
perc <- sums/length(count) # get share of images with quality issues
perc[perc == 0] <- NA
saveRDS(perc, "/home/veronika/CPD/data/BU_1981_2018/shareQF_complSeries_BU.rds")

# Quick Visualization to check:
#******************************
# library(stars)
# library(tidyverse)
# sp <- read_stars("/home/veronika/CPD/data/spatial_BU.tif")
# sp_bu <- transmute(sp, lai = as.vector(perc))
# plot(sp_bu)

