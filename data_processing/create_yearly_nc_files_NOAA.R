# read binary files from Chi Chen and Ranga, and create yearly nc files, 
# should have same style as MPI files from Uli Weber
# FILES CREATED ARE NOT IN ORDER (WITHIN EACH YEAR)
# since we take yearly averages and quantiles, this doesnt really matter
# DATA HAS FACTOR 1000, this saves memory # IGNORE
# (for true values divide by 1000)

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
  
  # from readme with bin files: 
  # Filepath = ‘Please_Give_FilePath.bin’;
  # fid = fopen(Filepath,'r');
  # data = fread(fid,[2160,4320],'int16',’ ieee-be’);
  # fclose(fid);
  # below we read the data as per the specifications in readme, i.e., 
  # big endian
  # integer
  # dim size
  # size = 2 bytes (16 bits)
  
  # years and dim for each data product
  product<- tail(unlist(strsplit(data_path, split = "/")), n = 1)
  if(product == "AVHRR") {
    years<- 1981:2018
    data_dim<- c(4320, 2160)
    data_na<- -32768
    files<- files[grep("BULAI", files)] # the other option is percentile
  } else if(product == "MODIS") {
    years<- 2000:2019
    data_dim<- c(3600, 7200)
    data_na<- -32768
  }
  
  print(paste(product, range(years)[1], range(years)[2], data_dim[1], 
               data_dim[2], sep = ", "))
  
  create_dir<- !dir.exists(paste0(data_path,"/", res_dir))
  if(create_dir) dir.create(paste0(data_path, "/", res_dir))
  
  create_parallel<- function(year, ...){
    files_in_year<- files[grep(year, files)]
    data<- array(NA, dim = c(length(files_in_year), data_dim[1]*data_dim[2]))
    print(files_in_year)# check years are good
    for(i in 1:length(files_in_year)){
      tmp<- readBin(paste0(data_path, "/", files_in_year[i]), what = "integer", 
                   n = (data_dim[1]*data_dim[2]), size = 2, endian = "big") 
      tmp[tmp == data_na]<- NA
      data[i,]<- tmp
    }
    data<- apply(data, 2, yearly_aggregating_fn, ...)
    data<- data/1000 # scaling it back to original value range
    data<- matrix(data, nrow = data_dim[1], ncol = data_dim[2])
    # image(data) # looks good
    
    save(data, file = paste0(data_path,"/", res_dir, "/LAI.", product, ".", data_dim[1], ".", 
    data_dim[2], ".", year, ".nc"))
    # not an nc file but this way it integrates with already existing code 
    # to process later 
    return()
  }
  future_lapply(years[1]:range(years)[2], create_parallel, ...)
  future:::ClusterRegistry("stop")
}

data_files<- c("data/CHEN_RANGA/AVHRR", "data/CHEN_RANGA/MODIS")

# data_path<- data_files[2] # for testing...
mean_fn<- function(x) mean(x, na.rm = TRUE)
quantile_fn<- function(x, probs) quantile(x, probs = probs, na.rm = TRUE)
create_yearly_nc_files_chen_ranga(data_files[1], yearly_aggregating_fn = mean_fn, res_dir = "yearly_mean")
create_yearly_nc_files_chen_ranga(data_files[1], yearly_aggregating_fn = quantile_fn, probs = .1, res_dir = "yearly_q_10")
create_yearly_nc_files_chen_ranga(data_files[1], yearly_aggregating_fn = quantile_fn, probs = .25, res_dir = "yearly_q_25")
create_yearly_nc_files_chen_ranga(data_files[1], yearly_aggregating_fn = quantile_fn, probs = .5, res_dir = "yearly_q_50")
create_yearly_nc_files_chen_ranga(data_files[1], yearly_aggregating_fn = quantile_fn, probs = .75, res_dir = "yearly_q_75")
create_yearly_nc_files_chen_ranga(data_files[1], yearly_aggregating_fn = quantile_fn, probs = .9, res_dir = "yearly_q_90")

create_yearly_nc_files_chen_ranga(data_files[2], yearly_aggregating_fn = mean_fn, res_dir = "yearly_mean")
create_yearly_nc_files_chen_ranga(data_files[2], yearly_aggregating_fn = quantile_fn, probs = .1, res_dir = "yearly_q_10")
create_yearly_nc_files_chen_ranga(data_files[2], yearly_aggregating_fn = quantile_fn, probs = .25, res_dir = "yearly_q_25")
create_yearly_nc_files_chen_ranga(data_files[2], yearly_aggregating_fn = quantile_fn, probs = .5, res_dir = "yearly_q_50")
create_yearly_nc_files_chen_ranga(data_files[2], yearly_aggregating_fn = quantile_fn, probs = .75, res_dir = "yearly_q_75")
create_yearly_nc_files_chen_ranga(data_files[2], yearly_aggregating_fn = quantile_fn, probs = .9, res_dir = "yearly_q_90")

# interdecile range
interdecile_fn<- function(x) quantile(x, probs = c(0.90), na.rm = TRUE) - quantile(x, probs = c(0.10), na.rm = TRUE)
create_yearly_nc_files_chen_ranga(data_files[1], yearly_aggregating_fn = interdecile_fn, res_dir = "yearly_interdecile_range")


