# create yearly median aggregates for nc data NOAA CDR
#****************************************

create_yearly_aggregate<- function(files, aggregating_fn, res_dir, i){
  curr_file = files[i]
  if(i>=2) prev_file<- files[(i-1)]
  if(i == 1) prev_file<- NA
  # note aggregating_fn must also have ... as input  - required by clusterR fn
  # beginCluster(n = NCORES)
  # on.exit(endCluster())

  create_dir<- !dir.exists(paste0(data_path,"/", res_dir))
  if(create_dir) dir.create(paste0(data_path, "/", res_dir))

  data_info<- unlist(base::strsplit(curr_file, split = "[.]"))
  # here, ALWAYS element [1] is var_name and element [length-1] is year
  year<- data_info[length(data_info) - 1]

  monthly_data<- brick(paste0(data_path, "/", MONTHLY_RES_DIR, "/", curr_file))
  if(i>=2) prev_monthly_data<- brick(paste0(data_path, "/", MONTHLY_RES_DIR, "/", prev_file))
  if(i==1) prev_monthly_data<- 0 # subtracting NA from  a number gives NA so we put 0 instead
  # X_PREV<- prev_monthly_data %>% mean(na.rm =TRUE)
  tmp<- monthly_data-prev_monthly_data
  rm(monthly_data)
  rm(prev_monthly_data)
  tmp<- tmp %>% mean(na.rm = TRUE)
  # tmp2<- tmp2 %>% mean(na.rm = TRUE)

  writeRaster(tmp, filename = paste0(res_dir, "/", curr_file))
  return(paste0("file saved as: ", res_dir, "/", curr_file))
}


clustermq_fn<- function(data_path){
  source("CPD/lib/source_all.R")
  tempfile(tmpdir="/tmp")
  NCORES<- 5
  files<- list.files(data_path,pattern = "\\.nc$")
  VAR_NAME<- unlist(base::strsplit(files[1], split = "[.]"))[1]
  MONTHLY_RES_DIR<- "monthly_mean"

  # for(i in 1:length(files)) create_yearly_aggregate(data_path = data_path,
  #                                                    aggregating_fn = deseasonalized_yearly_average, res_dir = "deseasonalized_yearly_average_2", i)
  mclapply_fn<- function(i, data_path){
    create_yearly_aggregate(data_path = data_path, aggregating_fn = deseasonalized_yearly_average, res_dir = "deseasonalized_yearly_average_2", i)
  }


  mclapply(1:length(files), mclapply_fn, data_path = data_path, mc.cores = NCORES)
}

source("lib/source_all.R")
NCORES<- 5
data_dirs<- c("data/gimms_lai_3gv1", "data/LAI_FAPAR_v2019", "data/LTDR/v5/FAPAR",
              "data/LTDR/v5/LAI", "data/LTDR/v5/NDVI", "data/MOD15A6_006", "data/PROBA_V")
# data_dirs<- data_dirs[c(2,3,4,5)]
# data_path<- data_dirs[2]

library(clustermq)
tempfile(tmpdir="/tmp")

Q(clustermq_fn,
  data_path = data_dirs,
  n_jobs = length(data_dirs),
  template = list(n_cpus = NCORES,
                  partition = "all",
                  job_name = "processing",
                  log_file = "/home/jose/LAI/lai_processing%a.txt"))

# clustermq_fn(data_dirs[1])
# clustermq_fn(data_dirs[2])
# clustermq_fn(data_dirs[3])
# clustermq_fn(data_dirs[4])
# clustermq_fn(data_dirs[5])
# clustermq_fn(data_dirs[6])
# clustermq_fn(data_dirs[7])

# read the creeated file
# a<- brick(paste0(data_path, "/", res_dir, "/", files[9]))

# # mclapply version for when servers are busy
# tempfile(tmpdir="/tmp")
#
# source("lib/source_all.R")
#
# data_dirs<- c("data/gimms_lai_3gv1", "data/LAI_FAPAR_v2019", "data/LTDR/v5/FAPAR",
#               "data/LTDR/v5/LAI", "data/LTDR/v5/NDVI", "data/MOD15A6_006", "data/PROBA_V")
#
#
# mclapply_fn<- function(i, data_path){
#   create_yearly_aggregate(data_path = data_path, aggregating_fn = deseasonalized_yearly_average, res_dir = "deseasonalized_yearly_average", i)
#   print(i)
# }
# data_path<- data_dirs[5]
# files<- list.files(data_path,pattern = "\\.nc$")
# VAR_NAME<- unlist(base::strsplit(files[1], split = "[.]"))[1]
#
# MONTHLY_RES_DIR<- "monthly_mean"
#
#
# mclapply(1:length(files), mclapply_fn, data_path = data_path, mc.cores = 19)




