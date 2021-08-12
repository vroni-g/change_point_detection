#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************


# Simulations to compare power for different time points of change and different
# magnitudes of change for MCUSUM
#****************************************************************************


# SIMULATION
#***********************************************************

library(magrittr)
#devtools::load_all("/home/veronika/CPD/change_point_detection/")
devtools::load_all()

# tp is timepoint of change
# n is length of time series
# mag is magnitude of change in parameter

cusum_power <- function(mag, tp, l, changetype,  nsim, alpha){

  if(changetype == "mean"){
    get_data <-  function(tp, l, mag){
      me <- sample(seq(0.2,2,0.2),size=1) # mean and sd values based on distribution of lai data
      sd <- sample(seq(0.02,0.4,0.02),size=1)
      d <- c(rnorm(n=tp, mean=me, sd=sd), rnorm(n=l-tp, mean=me*mag, sd=sd)) %>%
        ts(., frequency=1, start=c(0,1))
    }
  } else if (changetype == "variance"){
    get_data <-  function(tp, l, mag){
      me <- sample(seq(0.2,2,0.2),size=1)
      sd <- sample(seq(0.02,0.4,0.02),size=1)
      d <- c(rnorm(n=tp, mean=me, sd=sd), rnorm(n=l-tp, mean=me, sd=sd*mag)) %>%
        ts(., frequency=1, start=c(0,1))
    }
  } else if (changetype == "slope"){
    get_data <-  function(tp, l, mag){
      mag <- mag*sample(c(-1,1), 1)
      sl <- sample(seq(0.02,0.4,0.02),size=1)
      i1 <- sample(seq(0,2,0.2),size=1)
      i2 <- sl*tp + i1 + rnorm(n=1, mean=1, sd=1)
      d <- c((sl*c(1:tp) + i1 + rnorm(n=tp, mean=0.5, sd=2)),
             ((sl*mag)*c(1:(l-tp)) + i2 + rnorm(n=l-tp, mean=0.5, sd=2))) %>%
        ts(., frequency=1, start=1)
    }
  } else if (changetype == "intercept"){
    get_data <-  function(tp, l, mag){
      sl <- 0.12
      i <- sample(seq(0,2,0.2),size=1)
      if(mag > 1){
        i2 <- (sl*tp + i + rnorm(n=1, mean=1, sd=1))*mag
      } else {
        i2 <- i*mag
      }
      d <- c((sl*c(1:tp) + i + rnorm(n=tp, mean=0.1, sd=0.5)),
             (sl*c(1:(l-tp)) + i2 + rnorm(n=l-tp, mean=0.1, sd=0.5))) %>%
        ts(., frequency=1, start=c(0,1))
    }
  }

  get_result <- function(tp, l, mag, alpha){
    d <- get_data(tp, l, mag)
    res <- mcusum_function(d)
    return(res)
  }

  tps <- rep(tp, nsim)

  res <- purrr::map(tps, purrr::possibly(get_result, otherwise = NA), l=l, mag=mag, alpha=alpha) %>%
    unlist()
  res <- res[res!=999]
  pow <- sum((res <=alpha))/length(res)
  return(pow)

}


#t <- cusum_power(tp = 19, l = 38, changetype = "slope", mag = 2, nsim = 10, alpha=0.1)
#t

# change magnitude
#***********
# for every change repeat with new random data 50 times for a change of a certain
# magnitude at a fixed time step
# for each change magnitude record the number of times the change was detected

mags <- seq(0.1, 3, by=0.1)
mags <- mags[-10]

cat("Start: magnitude mean")
pow_mean_mag <- purrr::map(mags, cusum_power, tp = 17, l=37, changetype = "mean", nsim = 50, alpha = 0.1)%>%
  do.call(rbind, .) %>%
  as.data.frame()
pow_mean_mag$change_magnitude <- mags
saveRDS(pow_mean_mag, file = "mag_power_mcusum_meanchange_a10p_50each.rds")
cat("Saved magnitude  mean")

cat("Start: magnitude  variance")
pow_var_mag <- purrr::map(mags, cusum_power, tp = 17, l=37, changetype = "variance", nsim = 50, alpha = 0.1)%>%
  do.call(rbind, .) %>%
  as.data.frame()
pow_var_mag$change_magnitude <- mags
saveRDS(pow_var_mag, file = "mag_power_mcusum_varchange_a10p_50each.rds")
cat("Saved magnitude  variance")

cat("Start: magnitude  slope")
pow_slope_mag <- purrr::map(mags, cusum_power, tp = 17, l=37, changetype = "slope", nsim = 50, alpha = 0.1)%>%
  do.call(rbind, .) %>%
  as.data.frame()
pow_slope_mag$change_magnitude <- mags
saveRDS(pow_slope_mag, file = "mag_power_mcusum_slopechange_a10p_50each.rds")
cat("Saved magnitude  slope")

cat("Start: magnitude  inter")
pow_inter_mag <- purrr::map(mags, cusum_power, tp = 17, l=37, changetype = "intercept", nsim = 50, alpha = 0.1)%>%
  do.call(rbind, .) %>%
  as.data.frame()
pow_inter_mag$change_magnitude <- mags
saveRDS(pow_inter_mag, file = "mag_power_mcusum_interchange_a10p_50each.rds")
cat("Saved magnitude  inter")

# time point
#***********
# for every change type, repeat with new random data 1000 (?) times for a change at each time step
# for each time step record the number of times the change was detected

tps <-2:36
cat("Start: timepoints mean ", date())
pow_mean <- purrr::map(tps, cusum_power, l=37, changetype = "mean", mag = 0.5, nsim = 50, alpha = 0.1) %>%
  do.call(rbind, .) %>%
  as.data.frame()
saveRDS(pow_mean, file = paste0("time_power_cusum_meanchange_a10p_50each.rds"))
cat("Saved timepoints mean", date())

cat("Start: timepoints variance ", date())
pow_var <- purrr::map(tps, cusum_power, l=37, changetype = "variance", mag = 0.5, nsim = 50, alpha = 0.1) %>%
  do.call(rbind, .) %>%
  as.data.frame()
  saveRDS(pow_var, file = paste0("time_power_cusum_varchange_a10p_50each.rds"))
cat("Saved timepoints variance", date())

cat("Start: timepoints slope ", date())
pow_slope <- purrr::map(tps, cusum_power, l=37, changetype = "slope", mag = 2.5, nsim = 50, alpha = 0.1) %>%
do.call(rbind, .) %>%
as.data.frame()
saveRDS(pow_slope, file = paste0("time_power_cusum_slopechange_a10p_50each.rds"))
cat("Saved timepoints slope", date())

cat("Start: timepoints inter ", date())
pow_inter <- purrr::map(tps, cusum_power, l=37, changetype = "intercept", mag = 0.5, nsim = 50, alpha = 0.1) %>%
do.call(rbind, .) %>%
as.data.frame()
saveRDS(pow_inter, file = paste0("time_power_cusum_interchange_a10p_50each.rds"))
cat("Saved timepoints inter", date())

