############################################################
#### Preparation of the Protected Areas Layer
############################################################
# Description: I explore the land use types as downloaded from the Peace Parks
# website (http://new-ppfmaps.opendata.arcgis.com/datasets/ppf-protected-areas-
# detailed?geometry=-13.87%2C-25.558%2C69.846%2C-11.001). In addition, the
# classifications are simplified and reduced. Finally, protected areas are
# rasterized.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load packages
library(raster)
library(RColorBrewer)
library(rgdal)
library(rworldmap)
library(gdalUtils)
library(tidyverse)
library(cleangeo)
library(terra)

############################################################
#### Reduction to 3 Designations
############################################################
# Import file
prot <- readOGR("03_Data/01_RawData/PEACEPARKS/PPF_Protected_Areas_Detailed.shp")

# Crop the data according to the reference shapefile
r <- readOGR("03_Data/02_CleanData/00_General_Shapefile.shp")
prot <- crop(prot, r)

# Keep only the attributes of interest and rename them neatly
prot@data <- dplyr::select(prot@data
  , Name    = Name
  , Desig   = Designatio
  , IUCN    = IUCN
  , Country = Country
)

# Load the reclassification table
Desigs <- "03_Data/01_RawData/PEACEPARKS/Reclassification.csv" %>%
  read_csv() %>%
  set_names(., c("Nr", "Old", "New", "Comment"))

# Look at the table
Desigs

# Join the dataframes
prot@data <- left_join(prot@data, Desigs, by = c("Desig" = "Old"))

# Game reserves in Botswana serve the same purpose as national parks. Let's thus
# reclassify them accordingly
prot$New[prot$Desig == "Game Reserve" & prot$Country == "Botswana"] <-
  "National Park"

# Remove columns that we dont need anymore
prot@data <- prot@data %>% dplyr::select(-c("Desig", "Nr", "Comment"))

# Rename the remaining columns
prot@data <- prot@data %>% rename(Desig = New)

# Delete the objects and attributes that are not needed
prot <- subset(prot, prot$Desig != "Delete")

# Check validity of shapefile
gIsValid(prot, reason = TRUE)

# Plot the stuff with each designation in a different colour. Note that the
# country borders will look pretty nasty due to the low resolution.
map <- getMap("coarse")
u <- unique(prot$Desig)
m <- match(prot$Desig, u)
n <- length(unique(prot$Desig))
plot(prot, col = brewer.pal(n, "Greens")[m])
plot(map, add = TRUE)
text(subset(prot, prot$Desig == "National Park"), 'Name', cex = 0.5, halo = TRUE)

# Store the layer
writeOGR(prot
  , "03_Data/02_CleanData"
  , "02_LandUseTypes_Protected_PeaceParks(3Classes)"
  , driver = "ESRI Shapefile"
  , overwrite = TRUE
)

############################################################
#### Reduction to one Designation
############################################################
# We decided to only care about whether a region is protected or not. I will
# thus reclassify all designations to just protected
prot$Desig <- "Protected"

# Save the file
writeOGR(prot
  , "03_Data/02_CleanData"
  , "02_LandUseTypes_Protected_PeaceParks(1Class)"
  , driver = "ESRI Shapefile"
  , overwrite = TRUE
)

############################################################
#### Rasterize Protected Areas Layer
############################################################
# Now we also rasterize the layer to the 250m resolution reference layer
r250 <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Assign an arbitrary value of 1 to protected areas
prot$Value <- 1

# Rasterize protected areas
prot_r <- terra::rasterize(
    x           = vect(prot)
  , y           = rast(r250)
  , field       = "Value"
  , background  = 0
)

# Write raster to file
terra::writeRaster(
    prot_r
  , filename = "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(1Class).tif"
  , overwrite = T
)
