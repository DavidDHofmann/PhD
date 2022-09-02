################################################################################
#### Table of Inter-Patch Connectivity
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
setwd(wd)

# Load required packages
library(tidyverse)

# Load results on interpatch connectivity
dat <- read_rds("03_Data/03_Results/BootstrappedInterpatchConnectivity.rds")

sub <- subset(dat, FloodLevel == "Min")
sub <- select(sub, From, To, Freq)

t1 <- pivot_wider(sub, names_from = From, values_from = Freq)
t2 <- pivot_wider(sub, names_from = From, values_from = Freq)
