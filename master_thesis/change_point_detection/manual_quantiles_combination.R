#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# !! Input Data is not available in GitHub Repository !!


# Manual Combination procedure to correct for multiple testing
#*************************************************************

# Bin cluster sizes into classes and derive quantiles for each class
# "manual quantile thresholds"

library(tidyverse)

nperm <- 30

size <- c()
minp <- c()
stcs <- c()

for(i in 1:nperm){
  res <- readRDS(paste0("/home/veronika/CPD/MCUSUM/results/cluster/BU_MCUSUM_cluster_",i,".rds"))
  stcs <- c(stcs, res[[1]][[2]])
  clust <- res[[2]]
  size <- c(size, clust$cluster.count)
  minp <- c(minp, clust$cluster.min)
}
rm(clust, i,res)
or <- readRDS("/home/veronika/CPD/MCUSUM/results/cluster/BU_MCUSUM_cluster_original.rds")
stcs <- c(stcs, or[[1]][[2]])

# derive stcs threshold
thr <- sort(stcs, decreasing = T)[ceiling(length(stcs)*0.1)]

# get all clusters that are above stcs threshold
inds <- which(size >= thr)
if(length(inds)>0){
  size <- size[-inds]
  minp <- minp[-inds]
}

# define groups of cluster sizes:
hist(size, breaks = 50)
table(size)
bins <- c(1,2,6,11,21,51,101, 201, 301, 401, 501, 601, 701, 801, 909, 80000)
tags <- c("1", "2-5", "6-10", "11-20", "21-50", "51-100", "101-200", "201-300",
          "301-400", "401-500", "501-600", "601-700", "701-800", "801-900", ">900")

df <- tibble::tibble(size, minp)
group_tags <- cut(df$size,
                  breaks=bins,
                  include.lowest=TRUE,
                  right=FALSE,
                  labels=tags)
df$bins <- group_tags

quantiles <- df %>%
  group_by(.,bins) %>%
  summarize(quant10 = quantile(minp, probs = 0.1))
quantiles <- tibble::deframe(quantiles)
quantiles <- c(quantiles, "701-800"=0.002, "801-900"=0.002, ">900"=0.002)
quantiles

saveRDS(list(bins=bins, tags=tags, quantiles = quantiles), "/home/veronika/CPD/MCUSUM/manual_quantiles/quantiles_nperm30.rds")

# using all data of all permutations, threshold differs from 0.002 in cluster sizes 1, 2-5, 6-10, 301-400, 601-700
# if I don't use all data than for some cluster sizes the threshold might be missing..


# check clusters in original data
#********************************************************************************
clus <- or[[2]]
sum(clus$cluster.count)

size_or <- or[[2]]$cluster.count
minp_or <- or[[2]]$cluster.min

df_or <- tibble::tibble(size =size_or, minp=minp_or)
group_tags <- cut(df_or$size,
                  breaks=bins,
                  include.lowest=TRUE,
                  right=FALSE,
                  labels=tags)
df_or$bins <- group_tags

# check overall number of cluster and pixel
#*******************************************

# only size and 0.002:
inds_minP <- which(clus$cluster.min==0.002)
inds_stcs <- which(clus$cluster.count>=thr)
inds_simple <- unique(c(inds_minP, inds_stcs)) # 193 cluster
sum(clus$cluster.count[inds_simple]) # 36026 pixel
hist(clus$cluster.count[inds_simple])

# size and quantile threshold
df_or <- mutate(df_or, thres = quantiles[bins], ind = row_number())
inds_adj <- df_or[df_or$minp<=df_or$thres, ]$ind

hist(clus$cluster.count[inds_adj], breaks = 50)
table(clus$cluster.count[inds_adj])
smalls <- clus$cluster.count[inds_adj]
sum(smalls<10)/length(smalls)

inds_adj <- unique(c(inds_adj, inds_stcs)) # 352 cluster
inds_adj_df <- filter(df_or, ind %in% inds_adj)
cnts <-inds_adj_df %>% count(bins)

clust_mat <- clus[[1]]
mq <- clust_mat
mq[mq %in% inds_adj] <- -10
mq[mq!=-10] <- 0
mq <- dplyr::na_if(mq, 0)
mq[!is.na(mq)] <- 1
table(mq)
saveRDS(list(mat = mq, inds = inds_adj), "master_thesis/results/sig_nperm30_combined_stcsmq.rds")


# barplot of significant MCUSUM cluster and their size distribution
#******************************************************************
pdf_fn <- paste0("MCUSUM_clustersize_barplot.pdf")
pdf(pdf_fn, height=9.3, width=14)
barplot(height = cnts$n, names = cnts$bins, xlab = "Cluster size groups",
        ylab = "Count of significant cluster", ylim = c(0,120),
        main = "Sizes of significant MCUSUM cluster", cex.lab=1.5, cex.main=1.8)
dev.off()
#******************************************************************


all(inds_simple %in% inds_adj)
hist(clus$cluster.count[inds_adj])
sum(clus$cluster.count[inds_adj], na.rm = T) # 37507 pixel, only 1500 pixel more..

sum(clus$cluster.count[inds_stcs], na.rm = T)

# clusters that are additional in the adjusted version
inds_new <- setdiff(inds_adj, inds_simple)
hist(clus$cluster.count[inds_new])
table(clus$cluster.count[inds_new])


# # check classes with different thresholds
# #*******************************************
# six_hun <- filter(df_or, bins == "601-700")
# nrow(filter(six_hun, minp <= 0.002))
# nrow(filter(six_hun, minp <= quantiles["601-700"]))
# 
# three_hun <- filter(df_or, bins == "301-400")
# nrow(filter(three_hun, minp <= 0.002))
# nrow(filter(three_hun, minp <= quantiles["301-400"]))
# 
# six_ten <- filter(df_or, bins == "6-10")
# nrow(filter(six_ten, minp <= 0.002))
# nrow(filter(six_ten, minp <= quantiles["6-10"]))
# 
# two_five <- filter(df_or, bins == "2-5")
# nrow(filter(two_five, minp <= 0.002))
# nrow(filter(two_five, minp <= quantiles["2-5"]))
# 
# one <- filter(df_or, bins == "1")
# nrow(filter(one, minp <= 0.002))
# nrow(filter(one, minp <= quantiles["1"]))

