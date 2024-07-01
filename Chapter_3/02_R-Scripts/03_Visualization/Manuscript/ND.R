################################################################################
#### Normalized Differences Table
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
dat <- read_csv("03_Data/01_RawData/ND.csv")

# Make a nice table
kbl(dat[-1, ], booktabs = T, format = "latex", digits = 3, align = "lcc") %>%
  add_header_above(c("", "Bands" = 2)) %>%
  writeLines("04_Manuscript/Figures/ND.tex")
