################################################################################
#### Define Main Study Area
################################################################################
# Description: Define main stuy area of the project. Create shapefile for it and
# visualize all.

# Clean environment
rm(list = ls())

# Load required packages
library(raster)            # For handling spatial data
library(rgdal)             # For loading and storing spatial data
library(rgeos)             # For manipulating spatial data
library(tidyverse)         # For data wrangling
library(sf)                # For plotting spatial data

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

################################################################################
#### Prepare Data
################################################################################
# Prepare polygon for study area
ext <- extent(22, 27, -21, -18)
ext <- as(ext, "SpatialPolygons")
crs(ext) <- CRS("+init=epsg:4326")

# Let's load one of the floodmaps and extend it to the desired extent
r <- raster("03_Data/01_RawData/DAVID/ReferenceRaster.tif")
r <- crop(r, ext, snap = "out")

# Check the resolution of the final map (should be 250 meters)
res(r) * 111000

# Again create a shapefile that perfectly aligns with the borders of the raster
s <- extent(r)
s <- as(s, "SpatialPolygons")
crs(s) <- crs(r)

# Add some artificial data
s$Name <- "StudyArea"

# Let's load and clean the shapefile of africa
africa <- readOGR("03_Data/01_RawData/ESRI/Africa.shp")

# Remove the small islands (only keep africa + madagascar)
keep <- aggregate(africa, dissolve = T)
keep <- gBuffer(keep, width = 0.1)
keep <- disaggregate(keep)
keep$Area <- gArea(keep, byid = T)
keep <- keep[keep$Area %in% sort(keep$Area, decreasing = T)[1:2], ]
africa <- africa[keep, ]

# Let's also load some water areas (these are only for plotting)
water <- readOGR("03_Data/01_RawData/GEOFABRIK/Water.shp")

# Specify the areas that we want to keep for plotting later
object1 <- water[grepl(water@data$name, pattern = "Okavango.*Delta"), ][2, ]
object2 <- water[grepl(water@data$name, pattern = "Linyanti.*Delta"), ]
object3 <- water[grepl(water@data$name, pattern = "Garangwe.*Pan"), ][2, ]
object4 <- water[grepl(water@data$name, pattern = "Lake.*Ngami"), ]

# Put all objects together
water <- rbind(object1, object2, object3, object4)

# Crop data to the study area
water <- crop(water, s)

# Load data on protected Areas
prot <- readOGR("03_Data/01_RawData/PEACEPARKS/PPF_Protected_Areas_Detailed.shp")

# Use the reference shapefile to crop the areas
prot <- crop(prot, s)

# Keep only the attributes of interest and rename them nicely
prot@data <- dplyr::select(prot@data
  , Name    = Name
  , Desig   = Designatio
  , IUCN    = IUCN
  , Country = Country
)

# We now want to simplify the protection categories. We created a
# reclassification table for this, so let's use it
desigs <- read_csv("03_Data/01_RawData/PEACEPARKS/Reclassification.csv")
names(desigs) <- c("Nr", "Old", "New", "Comment")

# Join the dataframes
prot@data <- left_join(prot@data, desigs, by = c("Desig" = "Old"))

# Check out the distribution of the new categories
table(prot$New)

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

# Let's assign a value to each class
values <- data.frame(
    Desig  = c("Forest Reserve", "Protected", "National Park")
  , Values = 1:3
)

# Join the values to the shapefile
prot@data <- left_join(prot@data, values, by = "Desig")
prot$Desig <- factor(prot$Desig, levels = c("National Park", "Protected", "Forest Reserve"))

# Finally, we can load some data on the distribution of wild dogs
dogs <- readOGR("03_Data/01_RawData/IUCN/data_0.shp")

# Visualize
ggplot() +
  geom_sf(data = st_as_sf(africa), fill = "gray", lwd = 0) +
  geom_sf(data = st_as_sf(dogs), fill = "gray30", col = NA) +
  geom_sf(data = st_as_sf(s), fill = "red", col = "red", lwd = 0.5, alpha = 0.3, lty = 2) +
  theme_minimal()
ggplot() +
  geom_sf(data = st_as_sf(s)) +
  geom_sf(data = st_as_sf(prot), aes(fill = Desig), col = NA, alpha = 0.7) +
  geom_sf(data = st_as_sf(water), fill = "cornflowerblue", col = NA) +
  scale_fill_brewer(palette = "Greens", direction = -1, name = "Protection Status") +
  coord_sf(xlim = c(22, 27), ylim = c(-21, -18)) +
  theme_minimal()

################################################################################
#### Store Cleaned Data
################################################################################
# Let's make sure the crs of all data is the same
all <- list(r, s, water, prot, dogs)
lapply(all, crs)

# Store the data to file
writeRaster(r
  , filename  = "03_Data/02_CleanData/00_General_Raster.tif"
  , overwrite = T
)
writeOGR(s
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_Shapefile"
  , driver    = "ESRI Shapefile"
  , overwrite = T
)
writeOGR(water
  , dsn       = "03_Data/02_CleanData"
  , layer     = "03_LandscapeFeatures_MajorWaters_GEOFABRIK.shp"
  , driver    = "ESRI Shapefile"
  , overwrite = T
)
writeOGR(prot
  , dsn       = "03_Data/02_CleanData"
  , layer     = "02_LandUse_ProtectedAreas_PEACEPARKS.shp"
  , driver    = "ESRI Shapefile"
  , overwrite = T
)
writeOGR(dogs
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_WildDogs_IUCN"
  , driver    = "ESRI Shapefile"
  , overwrite = T
)
