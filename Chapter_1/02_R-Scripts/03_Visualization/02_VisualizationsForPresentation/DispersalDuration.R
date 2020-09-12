################################################################################
#### Plot of Dispersal Durations
################################################################################
# Description: A simple plot of the Dispersal Durations

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)  # For data wrangling
library(lubridate)  # To handle dates nicely
library(ggpubr)     # For nice plots

################################################################################
#### Data Preparation
################################################################################
# Load required data
cut <- read_csv("03_Data/02_CleanData/00_General_Dispersers_Popecol_CutoffDates.csv")

# Add Two hours to the timestamps
cut$StartDate <- cut$StartDate + hours(2)
cut$EndDate <- cut$EndDate + hours(2)

# Calculate dispersal duration
durations <- cut %>%
  mutate(DispersalDuration = EndDate - StartDate) %>%
  group_by(DogName) %>%
  summarize(DispersalDuration = as.numeric(sum(DispersalDuration)))

# Visualize
pdf("test.pdf", width = 9, height = 5)
ggdensity(
    durations
  , x     = "DispersalDuration"
  , fill  = "orange"
  , color = "orange"
  , add   = "mean"
  , rug   = TRUE
  , xlab  = "Dispersal Duration (days)"
  , ylab  = "Density"
)
dev.off()
