################################################################################
#### Preparation of Human Influence Layers
################################################################################
# Description: Preparation of the human influence layers. Details for the
# generation of those layers is given in Hofmann et al. 2021.

# Clear R's brain
rm(list = ls())

# Load required packages
library(terra)     # To handle spatial data

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Import the human influence layer and the reference raster
human <- rast("03_Data/01_RawData/DAVID/HumanInfluence.grd")
r <- rast("03_Data/02_CleanData/Raster.tif")

# Crop humans to reference raster
cat("Preparing human influence layer...\n")
human <- crop(human, r)

# We will only work with the 5km buffered data
human <- human$Buffer_5000
names(human) <- "HumanInfluence"

# Write the raster to file
writeRaster(human, "03_Data/02_CleanData/Humans.tif", overwrite = T)

# Let's also load some roads that we can plot later
roads <- vect("03_Data/01_RawData/GEOFABRIK/Roads.shp")

# Check out the description of the classes as derived from the OSM webpage
# (https://wiki.openstreetmap.org/wiki/Key:highway).
legend <- read.csv2(
    file             = "03_Data/01_RawData/GEOFABRIK/RoadsDescription.csv"
  , stringsAsFactors = F
)

# Keep only the largest roads (1-4) and their links (9-12)
roads <- subset(roads, roads$fclass %in% legend$Value[c(1:4, 9:12)])

# Crop roads to our study extent
s <- vect("03_Data/02_CleanData/Shapefile.gpkg")
roads <- crop(roads, ext(s))

# Save the cropped shapefile
writeVector(roads, "03_Data/02_CleanData/Roads.gpkg", overwrite = T)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/03_HumanInfluence_SessionInfo.rds")
cat("Done :)\n")
