############################################################
#### Preparation of the Reference Raster (250m)
############################################################
# Description: In order to facilitate the manipulation of other raster files I
# create a reference raster according to which I can crop and resample all other
# rasters. This will ensure that all rasters ultimately have the same
# resolution, origin and extent, which is necessary if we want to use the
# rasters for predictions.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(raster)
library(rgdal)
library(rworldmap)
library(davidoff)

# Load in the KAZA shapefile and use it as a reference
kaza <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")

# We want an extent that is slightly bigger than the KAZA extent. Let's
# manipulate the numbers to get the desired extent.
new_ext <- extent(kaza) + c(-0.5, +0.25, -0.63, +0.4)

# Now we create a raster that is defined by this new extent
r <- raster(new_ext)
crs(r) <- CRS("+init=epsg:4326")

# Fill the raster with some random values (nicer to plot)
values(r) <- runif(ncell(r))

# Plot the raster and the kaza together.
plot(r)
plot(kaza, border = "red", add = TRUE)

# Now create a 250m raster with the same projection and extent
r250  <- raster(r)

# Using the function we can now enter the desired resoluton
resolution250 <- metersToDegrees(250)

# And adjust the rasters resolution accordingly.
res(r250) <- resolution250

# We also fill it with random data so we can eventually plot the layer
values(r250) <- runif(ncell(r250))

# Visualize again
plot(r250)
plot(kaza, add = T, border = "red")

# Save the rasters to file
writeRaster(r250, "03_Data/02_CleanData/00_General_Raster250.tif", overwrite = TRUE)
