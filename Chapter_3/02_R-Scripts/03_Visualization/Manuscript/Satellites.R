################################################################################
#### Covariate Table
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(tidyverse)    # To wrangle data
library(kableExtra)   # To generate a nice table

# Load covariate table
dat <- read_csv("03_Data/01_RawData/Satellites.csv")

# Make a nice table
kbl(dat, booktabs = T, format = "latex", escape = F, align = "llcc") %>%
  writeLines("04_Manuscript/Figures/Satellites.tex")
