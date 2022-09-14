################################################################################
#### Source Areas
################################################################################
# Description: Generation of source area polygons

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")
# setwd("C:/Users/david/switchdrive/University/15. PhD/Chapter_8")

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
roads <- vect("03_Data/02_CleanData/Roads.shp")

# Buffer polygons slightly
river_buff <- buffer(river, width = 2000)
water_buff <- buffer(water, width = 500)

# Simplify geometries
river_buff <- aggregate(river_buff)
water_buff <- aggregate(water_buff)

# Make sure geometries are valid
river_buff <- makeValid(river_buff)
water_buff <- makeValid(water_buff)

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
plot(roads, add = T, col = "black")
plot(water_buff, add = T, border = "gray")
plot(river_buff, add = T, border = "gray")
plot(maun, add = T, col = adjustcolor("red", alpha.f = 0.2), border = NA)

# This looks good. Let's buffer the rivers and generate some potential source
# areas
src <- ell - water_buff - river_buff - maun
src <- disagg(src)
plot(src)

# We only want to keep the biggest 5 source areas!
src$Area <- expanse(src)
src      <- src[src$Area %in% sort(src$Area, decreasing = T)[1:5], ]
src$ID   <- 1:length(src)

# Crop roads
roads_crop <- crop(roads[roads$fclass == "trunk", ], ext(22, 24, -21, -19.5))

# Let's visualize the remaining areas
plot(src, col = adjustcolor(sample(rainbow(length(src))), alpha.f = 0.2))
plot(roads_crop, add = T)
text(src, "ID")

# Let's remove the part of the source areas in the south that intersect with the
# road
src      <- src - buffer(roads_crop, width = 1000)
src      <- disagg(src)
src$Area <- expanse(src)
src      <- src[src$Area %in% sort(src$Area, decreasing = T)[1:5], ]
src$ID   <- 1:length(src)

# Let's visualize the remaining areas
plot(src, col = adjustcolor(sample(rainbow(length(src))), alpha.f = 0.2))
plot(roads_crop, add = T)
text(src, "ID")

# We also want to generate source points representing immigrants into the
# ecosystem
out <- buffer(ell, width = -5000)
out <- src - out
src <- src - out

# We use this as one source area
out$ID <- max(src$ID) + (1:length(out))
out$Area <- expanse(out)
src <- rbind(src, out)

# Plot again
plot(src, col = adjustcolor(sample(rainbow(length(src))), alpha.f = 0.2))
text(src, "ID")

# Store them
writeVector(src, "03_Data/02_CleanData/SourceAreas.shp", overwrite = T)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/03_SourceAreas.rds")

# Print to terminal
cat("Done :) \n")
