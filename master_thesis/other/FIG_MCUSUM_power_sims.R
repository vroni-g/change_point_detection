#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# Input Data is not available in GitHub Repository but can be produced with
# MCUSUM_power_sims.R

# Plot simulation results
#****************************
library(tidyverse)

#**************** MAGNITUDE *********************
#********************************************************

mag_mean <- readRDS("mag_power_mcusum_meanchange_a10p_50each.rds")
mag_var <- readRDS("mag_power_mcusum_varchange_a10p_50each.rds")
mag_inter <- readRDS("mag_power_mcusum_interchange_a10p_50each.rds")
mag_slope <- readRDS("mag_power_mcusum_slopechange_a10p_50each.rds")


mag <- tibble::tibble(Magnitude = mag_mean$change_magnitude, Mean = mag_mean$V1,
                             Variance = mag_var$V1, Intercept = mag_inter$V1,
                             Slope = mag_slope$V1)
mag_long <- mag %>%
  pivot_longer(Mean:Slope)
mag_long$value <- as.numeric(mag_long$value)

pdf("magnitude_power_mcusum.pdf", height=6, width=10)
print(
ggplot(data = mag_long, aes(x=Magnitude, y=value, col=name)) +
        geom_line() +
        theme( panel.background = element_blank(),
               panel.border = element_rect(colour = "black", fill = NA),
               panel.grid = element_blank(),
               axis.line = element_line(colour = "black"),
               legend.position = c(0.37, 0.9),
               legend.key = element_blank(),
               legend.key.size = unit(1,"line"),
               legend.text = element_text(size=7),
               legend.direction = "horizontal",
               legend.title = element_text(size=9)) +
        xlab('Magnitude of Change at Timepoint 19 of 38') +
        ylab('CUSUM Power') +
        scale_color_manual(values=c("#DE3B77","#1E88E5", "#FFC107","#383C18")) +
        labs(colour = "Breakpoint Type") +
        guides(colour = guide_legend(title.position="top", title.hjust = 0.5))
)
dev.off()


#**************** TIMEPOINT *********************
#********************************************************

# MCUSUM
#**********************************

time_mean <- readRDS("time_power_cusum_meanchange_a10p_50each.rds")
time_var <- readRDS("time_power_cusum_varchange_a10p_50each.rds")
time_inter <- readRDS("time_power_cusum_interchange_a10p_50each.rds")
time_slope <- readRDS("time_power_cusum_slopechange_a10p_50each.rds")


time<- tibble::tibble(Timepoint = c(2:36), Mean = time_mean$V1,
                             Variance = time_var$V1, Intercept = time_inter$V1,
                             Slope = time_slope$V1)
time_long <- time %>%
  pivot_longer(Mean:Slope)
time_long$value <- as.numeric(time_long$value)

pdf("time_power_mcusum_50.pdf", height=8, width=12)
print(
  ggplot(data = time_long, aes(x=Timepoint, y=value, col=name)) +
    geom_line() +
    theme( panel.background = element_blank(),
           panel.border = element_rect(colour = "black", fill = NA),
           panel.grid = element_blank(),
           axis.line = element_line(colour = "black"),
           legend.position = c(0.9, 0.8),
           legend.key = element_blank(),
           legend.key.size = unit(1,"line"),
           legend.text = element_text(size=9),
           legend.title = element_text(size=11)) +
    xlab('Timepoints of Change with Magnitude 0.5 (except slope: 2.5)') +
    ylab('Modified CUSUM Power') +
    scale_color_manual(values=c("#DE3B77","#1E88E5", "#FFC107","#383C18")) +
    labs(colour = "Types of Change") +
    guides(colour = guide_legend(title.position="top", title.hjust = 0.5))
)
dev.off()

