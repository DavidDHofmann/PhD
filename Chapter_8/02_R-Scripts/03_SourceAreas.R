################################################################################
#### Source Areas
################################################################################
# Description: Generation of source area polygons

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load required data
water <- rast("03_Data/02_CleanData/WaterCover.tif")
shrub <- rast("03_Data/02_CleanData/ShrubCover.tif")
prote <- vect("03_Data/02_CleanData/Protected.shp")

# We only want to keep the maximum extent scenario for now
water <- water[["max"]]
shrub <- shrub[["max"]]

# Subset to protected areas of interest (we will use them as source areas)
area1 <- prote[prote$Desig == "National Park" & prote$Name == "Moremi" | prote$Name %in% c("NG/33", "NG/34"), ]
area2 <- prote[prote$Desig == "National Park" & prote$Name == "Chobe", ]
area3 <- prote[prote$Desig == "National Park" & prote$Name == "Hwange", ]
area4 <- prote[prote$Name %in% c("NG/26", "NG/29"), ]

# Dissolve
area1 <- aggregate(area1, dissolve = T)
area4 <- aggregate(area4, dissolve = T)

# Split up Moremi into multiple polygons (based on the water)
mask <- water == 1
mask <- subst(mask, 0, NA)
mask <- as.polygons(mask)
area1 <- erase(area1, mask)
area1 <- disagg(area1)
area1$Area <- expanse(area1)
area1 <- area1[area1$Area %in% tail(sort(area1$Area), 2), ]
area1 <- buffer(area1, width = +2000)
area1 <- buffer(area1, width = -1500)
area1 <- disagg(area1)

# Same for NG/26 NG/29
area4 <- erase(area4, mask)
area4 <- disagg(area4)
area4$Area <- expanse(area4)
area4 <- area4[area4$Area == max(area4$Area), ]
area4 <- buffer(area4, width = +3000)
area4 <- buffer(area4, width = -1500)
area4 <- disagg(area4)

# Make sure area 1 does not overlap with area 2
area1 <- area1 - area2

# Put everything back together
areas <- rbind(area1, area2, area3, area4)
values(areas) <- data.frame(ID = 1:length(areas))

# Give them names
areas$Name <- c("Moremi East", "Moremi West", "Chobe", "Hwange", "West-Delta")

# Let's visualize
plot(shrub)
plot(areas, add = T)
text(areas, label = areas$Name, cex = 0.5)

# Storre the polygons to file
writeVector(areas, "03_Data/02_CleanData/SourceAreas.shp", overwrite = T)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/03_SourceAreas.rds")
