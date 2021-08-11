#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# !! Not all Input Data is available in GitHub Repository !!

# Get Subtrends and assign breakpoint types
#******************************************************
library(trend)
library(tidyverse)

sig_mat <- readRDS("master_thesis/results/sig_nperm30_combined_stcsmq.rds")[[1]]
table(sig_mat)
sig <- which(sig_mat == 1, arr.ind = TRUE)
locs <- readRDS("master_thesis/results/MCUSUM_BU_orig_pvals_ar_locs.rds")[[3]]
data_lai <- readRDS("/home/veronika/CPD/data/10000int_BU_LAI_yearlymean_1982_2018")
data_lai <- data_lai/10000

# compute p values of MK trend test
get_subtrends <- function(i, inds, data, locations){
  x <- data[inds[[i,1]],inds[[i,2]],]
  bpl <- locations[inds[[i,1]],inds[[i,2]]]
  x1 <- x[1:bpl-1]
  x2 <- x[bpl:length(x)]
  mk1 <- trend::mk.test(x1)
  mk2 <- trend::mk.test(x2)
  res <- c(p1 = mk1$p.value, p2 = mk2$p.value)
  return(list(res = res, x1 = x1, x2 = x2, loc = bpl))
}

# for significant MK trend results compute linear slope
get_slope <- function(i, inds, data, locations){
  p <- get_subtrends(i, inds, data, locations)

  p1 <- p[["res"]][["p1"]]
  if(!is.numeric(p1) | is.na(p1)) {
    res1 <- c(NA, NA)
  } else if(p[["res"]][["p1"]] <= 0.1){
    d <- p[["x1"]]
    m <- lm(d ~ c(1:length(d)))
    res1 <- c(pval1 = p[["res"]][["p1"]], slope1 = m[["coefficients"]][["c(1:length(d))"]])
  } else {
    res1 <- c(pval1 = p[["res"]][["p1"]], slope1 = NA)
  }

  p2 <- p[["res"]][["p2"]]
  if(!is.numeric(p2) | is.na(p2)){
    res2 <- c(NA,NA)
  } else if(p[["res"]][["p2"]] <= 0.1){
    d <- p[["x2"]]
    m <- lm(d ~ c(1:length(d)))
    res2 <- c(pval2 = p[["res"]][["p2"]], slope2 = m[["coefficients"]][["c(1:length(d))"]])
  } else {
    res2 <- c(pval2 = p[["res"]][["p2"]], slope2 = NA)
  }
  return(c(row = inds[[i,1]], col=inds[[i,2]],bpl = p[["loc"]],
           bpl_year = 1981+p[["loc"]], res1, res2))
}


# iterate over significant pixels and create a dataframe with columns: array indices, slope results
slopes_df <- lapply(1:(length(sig)/2), get_slope, inds = sig, data = data_lai, locations = locs)  %>%
  do.call(rbind, .) %>%
  as.data.frame()

# change types:
# green_green = greening to greening
# brown_brown = browning to browning
# green_brown = greening to browning
# brown_green = browning to greening
# non_green = insignificant to greening
# green_non = greening to insignificant
# non_brown = insignificant to browning
# brown_non = browning to insignificant
# non_non = insignificant to insignificant

slopes_df <- slopes_df %>%
  mutate(., changetype = case_when(
    is.na(slope1) & is.na(slope2) ~ "non_non", # type 1
    is.na(slope1) & slope2>0 ~ "non_green", # type 2
    is.na(slope2) & slope1>0 ~ "green_non", # type 3
    is.na(slope1) & slope2<0 ~ "non_brown", # type 4
    is.na(slope2) & slope1<0 ~ "brown_non", # type 5
    slope1>0 & slope2>0 ~ "green_green", # type 6
    slope1>0 & slope2<0 ~ "green_brown", # type 7
    slope1<0 & slope2<0 ~ "brown_brown", # type 8
    slope1<0 & slope2>0 ~ "brown_green", # type 9
  ))
slopes_df <- slopes_df %>%
  mutate(., changetype_int = case_when(
    changetype == "green_brown" ~ 1, # type 1
    changetype == "brown_green" ~ 2, # type 2
    changetype == "green_green"~ 3, # type 3
    changetype == "non_green" ~ 4, # type 4
    changetype == "green_non" ~ 5, # type 5
    changetype == "brown_brown" ~ 6, # type 6
    changetype == "non_brown" ~ 7, # type 7
    changetype == "brown_non" ~ 8, # type 8
    changetype == "non_non" ~ 9, # type 9
  ))

slopes_df$changetype_int <- as.factor(slopes_df$changetype_int)
str(slopes_df)
saveRDS(slopes_df, "master_thesis/results/sig_mq_df.rds")

# GN <- dplyr::filter(slopes_df, changetype == "green_non")
# table(GN$bpl_year)
# sum(table(GN$bpl_year)[13:17])/sum(table(GN$bpl_year))
# 
# NB <- dplyr::filter(slopes_df, changetype == "non_brown")
# table(NB$bpl_year)
# sum(table(NB$bpl_year)[1:5])/sum(table(NB$bpl_year))


