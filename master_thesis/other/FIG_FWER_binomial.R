#*******************************************************************************
# Master Thesis: "Correcting for Multiple Testing in Change Point Detection 
#  of Global Vegetation Trends"
#*******************************************************************************
# Veronika Grupp, M.Sc. Geoinformatics and Remote Sensing, FSU Jena, 2021
#*******************************************************************************

# generate figure of Family Wise Error Rate by number of independent simultaneous tests
#*******************************************************
library(tidyverse)
sizes <- 1:150
binoms <- lapply(sizes, function(x) 1 - pbinom(0, size=x, prob=0.05)) %>% unlist

pdf_fn <- "FWER_binomial.pdf"
pdf(pdf_fn, height=6, width=8)
plot(sizes,binoms, type = "l",
     xlab = 'Number of independent Tests',
     ylab = 'Familywise Error Rate (FWER)')
dev.off()
