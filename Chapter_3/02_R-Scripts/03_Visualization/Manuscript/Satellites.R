################################################################################
#### Covariate Table
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/02_Academia/02_PhD/Chapter_3"
wd <- "D:/SwitchDrive/02_Academia/02_PhD/Chapter_3"
setwd(wd)

# Load required packages
library(tidyverse)    # To wrangle data
library(kableExtra)   # To generate a nice table

# Load covariate table
dat <- read_csv("03_Data/01_RawData/Satellites.csv")

# Make a nice table
kbl(dat, booktabs = T, format = "latex", escape = F, align = "llcc") %>%
  writeLines("04_Manuscript/Figures/Satellites.tex")
