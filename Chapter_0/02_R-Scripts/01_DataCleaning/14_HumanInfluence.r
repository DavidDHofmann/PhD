############################################################
#### Creating the Human Influence Layers (Base, Density & Buffer)
############################################################
# Description: Here I combine the layers from Geofabrik, Croplands, Gabriele and
# Globelands to create two layers. First, a binary human influence layer that
# depicts in each cell whether there are humans or not. Second, I create a human
# density map, indicating whether human influence is strong or weak.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(raster)
library(rgeos)
library(rgdal)
library(tidyverse)
library(spatstat)
library(maptools)
library(gdalUtils)

# Make use of multicores
beginCluster()

############################################################
#### Human Base and Density Layers
############################################################
# Import data representing human influence
# buildings   <- shapefile("04_AnthropogenicFeatures_Buildings(Points)")
population  <- raster("03_Data/02_CleanData/04_AnthropogenicFeatures_HumanDensity_Facebook.tif")
farms_gabs  <- raster("03_Data/02_CleanData/04_AnthropogenicFeatures_Farms_Gabriele.tif")
farms_crops <- raster("03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_Croplands.tif")
farms_glob  <- raster("03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_Globelands.tif")
roads       <- raster("03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.tif")

# Combine all the layers (make sure that the farms are not double counted)
merged <- population + roads + max(farms_crops, farms_gabs, farms_glob)

# Let's check how the resulting values are distributed (we only care about
# values larger than 0)
hist(values(merged)[values(merged) > 0])

# Apparently there are some heavy outliers. Let's remove values beyond 50 (set
# them to 50)
values(merged)[values(merged) > 50] <- 50

# Look at the histogram again
hist(values(merged)[values(merged) > 0])

# Load the shapefiles for the areas in which we know there are only camps but no
# human influence otherwise
delete <- shapefile("03_Data/01_RawData/DAVID/AnthropogenicInfluence(Delete)")

# Use the shapefile to delete the buildings that are inside the polygon
merged <- mask(merged, delete, inverse = TRUE, updatevalue = 0)

# Save the raster to file
writeRaster(merged
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Density).tif"
  , overwrite = TRUE
)

# Reclassify the merged layer so it contains only 0's and 1's (Human influence
# no/yes)
merged <- raster::reclassify(merged, c(0, max(values(merged)), 1))
writeRaster(merged
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Base).tif"
  , overwrite = TRUE
)

############################################################
#### Buffered Human Density
############################################################
# Reload the human density layer
dens <- raster("03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Density).tif")

# Create focal matrices for the desired buffer (5km)
w5000 <- focalWeight(dens, d = 1 / 111 * 5.0, type = "circle")

# Apply the matrices to sum up focal cells
buff5000 <- focal(dens
  , w         = w5000 / max(w5000)
  , FUN       = sum
  , pad       = TRUE
  , padValue  = 0
)

# Take the log of the layer (add one to assure that we only get positive
# numbers)
log5000 <- log(buff5000 + 1)

# Plot the result
plot(log5000)

# Write the buffered layers to file
writeRaster(log5000
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Buffer5000).tif"
  , overwrite = TRUE
)

############################################################
#### Distance To Humans
############################################################
# We also need a raster that depicts the distance to the nearest human inhabited
# cell. Let's load the human base layer again and coerce inhabited cells to ppp
humansBase <- "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Base).tif" %>%
  raster() %>%
  rasterToPoints(., fun = function(x){x == 1}, spatial = TRUE) %>%
  spTransform(., CRS("+init=epsg:32734")) %>%
  as(., "ppp")

# Load the reference raster to prepare a distance to humans layer
DistanceToHumans <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Now replace the values of the raster with the distances to the nearest human
# inhabited cell
values(DistanceToHumans) <- DistanceToHumans %>%
  as(., "SpatialPoints") %>%
  spTransform(., CRS("+init=epsg:32734")) %>%
  as(., "ppp") %>%
  nncross(., humansBase) %>%
  .[["dist"]]

# Plot the result
plot(DistanceToHumans)

# Store the rasters to file
writeRaster(DistanceToHumans
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToHumans.tif"
  , overwrite = TRUE
)

# Terminate the cluster
endCluster()
