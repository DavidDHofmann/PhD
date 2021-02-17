################################################################################
#### Preparation of Data for Dryad
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load packages
library(tidyverse)

# Load GPS data of dispersers
data <- read_csv("03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv")
data <- subset(data, State == "Disperser")

# Generate random numbers
data$identifier <- sample(1:nrow(data), replace = F)

# Sort data by identifier
data <- arrange(data, identifier)

# Prepare the same dataframe with fewer columns
data_dryad <- select(data, identifier, lon = x, lat = y)

# Store dataframes
write_csv(data_dryad, "03_Data/03_Results/GPS_Dryad.csv")
write_csv(data, "03_Data/03_Results/GPS_Full.csv")
