#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# !! Input Data is not available in GitHub Repository !!


# ********* BU QUALITY FLAGS ******************
# ********************************************
library(tidyverse)
#res5 <- readRDS("/home/veronika/CPD/results/BU_LAI_tippet_nperm_1000_al5.rds")
res10 <- readRDS("/home/veronika/CPD/results/BU_LAI_tippet_nperm_1000_al10.rds") 
  
dataset <- "BU"
alpha_local <- 0.05
alpha_global <- 0.05

# for every cluster declared significant by Tippet retrieve QF values
# inds_tip5 <- which(res5$original_wt > quantile(res5$wt, probs = 0.95))
# cluster_mat5<- res5$original_cluster[["clusters"]] %>% 
#   as.matrix
# cluster_mat5[cluster_mat5== 0] <- NA
# cluster_mat5[!(cluster_mat5 %in% inds_tip5)] <- NA

inds_tip10 <- which(res10$original_wt > quantile(res10$wt, probs = 0.9))
cluster_mat10<- res10$original_cluster[["clusters"]] %>% 
  as.matrix
cluster_mat10[cluster_mat10== 0] <- NA
cluster_mat10[!(cluster_mat10 %in% inds_tip10)] <- NA

qf <- readRDS("/home/veronika/CPD/data/BU_1981_2018/shareQF_complSeries_BU.rds")

# numerical results per significant cluster
#************************************************************
# extract qf values for a given cluster index
median_qf <- function(ind, cluster_mat, qf){
  cluster_mat[cluster_mat != ind] <- NA
  cluster_mat[!is.na(cluster_mat)] <- TRUE
  vals <- cluster_mat*qf
  vals <- vals[!is.na(vals)]
  med <- median(vals)
  iq <- quantile(vals, probs = 0.75) - quantile(vals, probs = 0.25)
  return(c(med, iq))
}

# med_sig5 <- lapply(inds_tip5, median_qf, cluster_mat = cluster_mat5, qf = qf) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame()
# colnames(med_sig5) <- c("median","IQR")
# 
# clust_size5 <- res5$original_cluster$cluster.count[inds_tip5]
# med_sig5$size <- clust_size5
# 
# nrow(dplyr::filter(med_sig5, median > 0.4))/nrow(med_sig5)

med_sig10 <- lapply(inds_tip10, median_qf, cluster_mat = cluster_mat10, qf = qf) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()
colnames(med_sig10) <- c("median","IQR")

clust_size10 <- res10$original_cluster$cluster.count[inds_tip10]
med_sig10$size <- clust_size10

nrow(dplyr::filter(med_sig10, median > 0.4))/nrow(med_sig10)
#plot(med_sig$size, med_sig$median)
#plot(med_sig$size, med_sig$IQR)

library(ggplot2)
# Alpha = 0.05
#**********************
# pdf_fn <- paste0("TCF_MK_QF_clustsize_",dataset, 
#                  "_MK_1981-2018_al5_ag5.pdf")
# pdf(pdf_fn, height=5, width=6)
# print(
#   ggplot() + 
#   geom_point(data = med_sig5, aes(x = size, y = median*100, size = IQR*100), alpha = 0.7) +
#   theme(panel.grid = element_blank(),
#         panel.background = element_blank(),   
#         panel.border = element_rect(fill = NA)) +
#   xlab("Size of Significant Cluster") +
#   ylab("Median of Quality Issues %") +
#   labs(size = "IQR") +
#   annotate(geom="text", x=33000, y=70, label="Local alpha = 0.05", size = 3.2)+
#   annotate(geom="text", x=33000, y=66, label="Global alpha = 0.05", size = 3.2)
# )
# dev.off()

# Alpha = 0.1
#**********************
pdf_fn <- paste0("TCF_MK_QF_clustsize_",dataset, 
                 "_MK_1981-2018_al10_ag10.pdf")
pdf(pdf_fn, height=5, width=6)
print(
  ggplot() + 
  geom_point(data = med_sig10, aes(x = size, y = median*100, size = IQR*100), alpha = 0.7) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),   
        panel.border = element_rect(fill = NA)) +
  xlab("Size of Significant Cluster") +
  ylab("Median of Quality Issues %") +
  labs(size = "IQR") +
    annotate(geom="text", x=45000, y=68, label="Local alpha = 0.1", size = 3.2)+
    annotate(geom="text", x=45000, y=64, label="Global alpha = 0.1", size = 3.2)
)
dev.off()




