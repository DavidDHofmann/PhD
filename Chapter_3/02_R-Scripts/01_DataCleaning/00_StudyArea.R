################################################################################
#### Define Main Study Area
################################################################################
# Description: Define main stuy area of the project and create a shapefile for
# it.

# Clean environment
rm(list = ls())

# Load required packages
library(terra)            # For handling spatial data
library(tidyterra)        # For tidy terra handling
library(tidyverse)         # For data wrangling
library(sf)                # For plotting spatial data
library(ggpubr)            # To arrange multiple plots

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

################################################################################
#### Prepare Data
################################################################################
# Prepare a polygon for the study area
ext <- c(21, 27, -21, -17) %>%
  ext() %>%
  as.polygons() %>%
  set.crs("+init=epsg:4326")

# Let's load the original reference raster and crop it to the new extent
r <- "03_Data/01_RawData/DAVID/ReferenceRaster.tif" %>%
  rast() %>%
  crop(ext, snap = "out")

# Check the resolution of the final map (should be 250 meters)
res(r) * 111000

# Again create a shapefile that perfectly aligns with the borders of the raster
s <- r %>%
  ext() %>%
  as.polygons() %>%
  set.crs("+init=epsg:4326") %>%
  setValues(., values = data.frame(Name = "StudyArea"))

# Let's load and clean the shapefile of africa
cat("Preparing shapefiles for africa...\n")
africa <- vect("03_Data/01_RawData/ESRI/Africa.shp")
keep <- africa %>%
  buffer(width = 0) %>%
  aggregate() %>%
  buffer(width = 10) %>%
  disagg() %>%
  set.crs("+init=epsg:4326") %>%
  setValues(., values = data.frame(Area = expanse(.)))
keep <- keep[keep$Area %in% sort(keep$Area, decreasing = T)[1:2], ]
africa <- africa[keep, ]

# Let's also load some water areas (these are only for plotting)
cat("Preparing shapefiles of major water areas...\n")
water <- vect("03_Data/01_RawData/GEOFABRIK/Water.shp")

# Specify the areas that we want to keep for plotting later
object1 <- water[grepl(water$name, pattern = "Okavango.*Delta"), ][2, ]
object2 <- water[grepl(water$name, pattern = "Linyanti.*Delta"), ]
object3 <- water[grepl(water$name, pattern = "Garangwe.*Pan"), ][2, ]
object4 <- water[grepl(water$name, pattern = "Lake.*Ngami"), ]

# Put all objects together
water <- rbind(object1, object2, object3, object4)
plot(water, col = "cornflowerblue", border = NA)
plot(africa, add = T)

# Crop data to the study area
water <- crop(water, s)

# Load data on protected Areas
cat("Preparing shapefile of protected areas...\n")
prot <- "03_Data/01_RawData/PEACEPARKS/PPF_Protected_Areas_Detailed.shp" %>%
  vect() %>%
  crop(s) %>%
  select(
      Name    = Name
    , Desig   = Designatio
    , IUCN    = IUCN
    , Country = Country
  )

# We now want to simplify the protection categories. We created a
# reclassification table for this, so let's use it
desigs <- "03_Data/01_RawData/PEACEPARKS/Reclassification.csv" %>%
  read_csv() %>%
  setNames(c("Nr", "Old", "New", "Comment"))

# Join the dataframes
values(prot) <- left_join(values(prot), desigs, by = c("Desig" = "Old"))

# Game reserves in Botswana serve the same purpose as national parks. Let's thus
# reclassify them accordingly
prot$New[prot$Desig == "Game Reserve" & prot$Country == "Botswana"] <-
  "National Park"

# Remove columns that we dont need anymore
prot <- prot %>%
  select(-c("Desig", "Nr", "Comment")) %>%
  rename(Desig = New) %>%
  filter(Desig != "Delete")

# Let's assign a value to each class
values <- data.frame(
    Desig  = c("Forest Reserve", "Protected", "National Park")
  , Values = 1:3
)

# Join the values to the shapefile
values(prot) <- left_join(values(prot), values, by = "Desig")
prot <- mutate(prot, Desig = factor(Desig, levels = c("National Park", "Protected", "Forest Reserve")))

# Also create a rasterized version
protr <- rasterize(prot, r, field = "Values", background = 0)

# Load fillages
vill <- vect("03_Data/01_RawData/DAVID/Villages.shp")

# Finally, we can load some data on the distribution of wild dogs for
# visualization
dogs <- vect("03_Data/01_RawData/IUCN/data_0.shp")

# Prepare some plots
p1 <- ggplot() +
  geom_sf(data = st_as_sf(africa), fill = "gray90", col = "white", lwd = 0.1) +
  geom_sf(data = st_as_sf(dogs), fill = "gray80", col = NA) +
  geom_sf(data = st_as_sf(s), fill = "cornflowerblue", col = "cornflowerblue", lwd = 0.5, alpha = 0.5, lty = 2) +
  theme_void()
p2 <- ggplot() +
  geom_sf(data = st_as_sf(s), fill = "transparent") +
  geom_sf(data = st_as_sf(prot), aes(fill = Desig), col = NA, alpha = 0.7) +
  geom_sf(data = st_as_sf(water), fill = "cornflowerblue", col = NA) +
  scale_fill_brewer(palette = "Greens", direction = -1, name = "Protection Status") +
  coord_sf(xlim = c(22, 27), ylim = c(-21, -18)) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# Arrange plots
ggarrange(p1, p2, nrow = 2)

################################################################################
#### Store Cleaned Data
################################################################################
# Store the data to file
cat("Storing all data to file...\n")
writeRaster(r, "03_Data/02_CleanData/Raster.tif", overwrite = T)
writeVector(s, "03_Data/02_CleanData/Shapefile.gpkg", overwrite = T)
writeVector(africa, "03_Data/02_CleanData/Africa.gpkg", overwrite = T)
writeVector(water, "03_Data/02_CleanData/MajorWaters.gpkg", overwrite = T)
writeVector(prot, "03_Data/02_CleanData/Protected.gpkg", overwrite = T)
writeRaster(protr, "03_Data/02_CleanData/Protected.tif", overwrite = T)
writeVector(dogs, "03_Data/02_CleanData/WildDogs.gpkg", overwrite = T)
writeVector(vill, "03_Data/02_CleanData/Villages.gpkg", overwrite = T)

################################################################################
#### Session Information
################################################################################
# Create directory for session info
dir.create("02_R-Scripts/99_SessionInformation", showWarnings = F)
dir.create("02_R-Scripts/99_SessionInformation/01_DataCleaning", showWarnings = F)

# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/00_StudyArea.rds")
cat("Done :)\n")
