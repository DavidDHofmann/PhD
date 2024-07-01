################################################################################
#### Example Dispersal Tracks
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)
library(tidyverse)

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData")

# Load dispersal tracks
dat <- read_csv("00_General_Dispersers_Popecol.csv")
dat <- subset(dat, DogName %in% c("Abel"))

# Visualize them
ggplot(dat, aes(x = x, y = y, col = State, group = DogName)) +
  geom_path(lwd = 1) +
  scale_color_manual(values = c("orange", "white")) +
  theme_void()
ggsave("DispersalTrack.png", device = "png", bg = "transparent")
