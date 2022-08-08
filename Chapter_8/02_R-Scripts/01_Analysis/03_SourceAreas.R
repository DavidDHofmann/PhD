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

# We'll use geographical, in particular hydrological features to separate source
# areas
river <- vect("03_Data/02_CleanData/MajorRivers.shp")
river <- river[-2, ]
water <- vect("03_Data/02_CleanData/MajorWaters.shp")
water <- water[c(1, 2, 4), ]

# Buffer polygons slightly
river_buff <- buffer(river, width = 2000)
water_buff <- buffer(water, width = 500)

# Put them together
river_buff <- aggregate(river_buff)
water_buff <- aggregate(water_buff)

# We may also want to remove the close vicinity of Maun
maun <- rast("03_Data/02_CleanData/HumanInfluence.tif")
maun <- crop(maun, ext(23.3, 23.7, -20.5, -19.6))
maun <- classify(maun, rbind(c(0, 7, 0), c(7, Inf, 1)))
maun <- subst(maun, 0, NA)
maun <- as.polygons(maun)
maun <- buffer(maun, width = 7500)

# Drawn an ellipse over the okavango delta
ell <- ellipse(22.8, -19.3, 1.4, 1.25, phi = 45, spatial = T)
crs(ell) <- crs(water)

# Visualize it
plot(ell, col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
plot(water, add = T, col = "cornflowerblue", border = NA)
plot(river, add = T, col = "cornflowerblue")
plot(maun, add = T, col = "red")

# This looks good. Let's buffer the rivers and generate some potential source
# areas
ell <- ell - water_buff - river_buff - maun
ell <- disagg(ell)

# We only want to keep the biggest 5 source areas!
ell$Area <- expanse(ell)
ell <- ell[ell$Area %in% sort(ell$Area, decreasing = T)[1:5], ]
ell$ID <- 1:length(ell)

# Let's visualize the remaining areas
plot(ell, col = adjustcolor(sample(rainbow(length(ell))), alpha.f = 0.2))
text(ell, "ID")

# Store them
writeVector(ell, "03_Data/02_CleanData/SourceAreas.shp", overwrite = T)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/03_SourceAreas.rds")

# Print to terminal
cat("Done :) \n")
