############################################################
#### Preparation of all Plots for the Thesis
############################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
input   <- "/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/00_LargeExtent"
output  <- "/home/david/ownCloud/University/14. FS 19/Masterarbeit/05_Presentation"
setwd(input)

# Set a seed
set.seed(12345)

# Load required packages
library(NLMR)
library(lubridate)
library(gdistance)
library(scales)
library(raster)
library(tmap)
library(tidyverse)
library(glmmTMB)
library(tikzDevice)
library(Cairo)
library(xtable)
library(viridis)
library(rosm)
library(rgeos)
library(grid)
library(tmaptools)
library(latex2exp)
library(measurements)
library(gridExtra)
library(RColorBrewer)
library(lemon)
library(rasterVis)
library(prettymapr)
library(move)
library(moveVis)


################################################################################
#### NEEWWW
################################################################################
############################################################
#### Plotting Trajectories
############################################################
# Clear R's brain
rm(list = ls())

# Load required packages
options("rgdal_show_exportToProj4_warnings"="none")
library(raster)
library(tmap)
library(tidyverse)
library(Cairo)
library(viridis)
library(rosm)
library(rgeos)
library(tmaptools)
library(prettymapr)
library(rgeos)
library(rgdal)

# Set the working directory
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_0"
setwd(wd)

# Load tracks
tracks <- readOGR("03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).shp")

# Prepare the extent of our core area
bbox <- as(extent(c(23.25, 24, -19.8, -19.05)), "SpatialPolygons")
crs(bbox) <- CRS("+init=epsg:4326")
bbox <- as(bbox, "SpatialLines")

# Plot an empty map
setwd("/home/david/Schreibtisch")
png("99_Trajectories1.png", width = 1920, height = 1080)

# Get the bing map
prettymap(
  bmaps.plot(bbox(tracks)
    , type = "Aerial"
    , key = "Ai0A3oHx1c3rivWQfR2nAmSle5AVyo6RxjtzUFLkxTT2_qbxyNDSge0A5jxxnDgW"
    , stoponlargerequest = FALSE
  )
  , scale.style     = "ticks"
  , scale.label.col = "white"
  , scale.linecol   = "white"
  , scale.label.cex = 3
  , scale.lwd       = 5
)

# Store the plot
dev.off()

# Plot an empty map with the extent
png("99_Trajectories2.png", width = 1920, height = 1080)

# Get the bing map
prettymap(
  bmaps.plot(bbox(tracks)
    , type = "Aerial"
    , key = "Ai0A3oHx1c3rivWQfR2nAmSle5AVyo6RxjtzUFLkxTT2_qbxyNDSge0A5jxxnDgW"
    , stoponlargerequest = FALSE
  )
  , scale.style     = "ticks"
  , scale.label.col = "white"
  , scale.linecol   = "white"
  , scale.label.cex = 3
  , scale.lwd       = 5
)

# Add the extent
bbox %>%
  spTransform(., CRS("+init=epsg:3857")) %>%
  plot(
      col = "white"
    , add = TRUE
    , lwd = 5
  )

# Store the plot
dev.off()

# Plot all trajectories onto a bing map.
png("99_Trajectories3.png", width = 1920, height = 1080)

# Get the bing map
prettymap(
  bmaps.plot(bbox(tracks)
    , type = "Aerial"
    , key = "Ai0A3oHx1c3rivWQfR2nAmSle5AVyo6RxjtzUFLkxTT2_qbxyNDSge0A5jxxnDgW"
    , stoponlargerequest = FALSE
  )
  , scale.style     = "ticks"
  , scale.label.col = "white"
  , scale.linecol   = "white"
  , scale.label.cex = 3
  , scale.lwd       = 5
)

# Add the trajectories onto the satellite map
tracks %>%

  # Aggregate them by individual
  gLineMerge(., byid = T, id = tracks$id) %>%

  # Transform to the crs of the bmap
  spTransform(., CRS("+init=epsg:3857")) %>%

  # Plot the tracks
  plot(.
    , col = rev(viridis(length(.)))
    , add = TRUE
    , lwd = 3
  )

# Add the extent
bbox %>%
  spTransform(., CRS("+init=epsg:3857")) %>%
  plot(
      col = "white"
    , add = TRUE
    , lwd = 5
  )

# Store the plot
dev.off()

# We also want to plot all trajectories but with different colours for
# dispersers and residents
png("99_Trajectories4.png", width = 1920, height = 1080)

# Get the bing map
prettymap(
  bmaps.plot(bbox(tracks)
    , type = "Aerial"
    , key = "Ai0A3oHx1c3rivWQfR2nAmSle5AVyo6RxjtzUFLkxTT2_qbxyNDSge0A5jxxnDgW"
    , stoponlargerequest = FALSE
  )
  , scale.style     = "ticks"
  , scale.label.col = "white"
  , scale.linecol   = "white"
  , scale.label.cex = 3
  , scale.lwd       = 5
)

# Add the extent
bbox %>%
  spTransform(., CRS("+init=epsg:3857")) %>%
  plot(
      col = "white"
    , add = TRUE
    , lwd = 5
  )

# Add the resident trajectories onto the satellite map
tracks %>%

  # Subset to residents
  subset(., state == "Resident") %>%

  # Transform to the crs of the bmap
  spTransform(., CRS("+init=epsg:3857")) %>%

  # Plot the tracks
  plot(.
    , col = "white"
    , add = TRUE
    , lwd = 3
  )

# Add the resident trajectories onto the satellite map
tracks %>%

  # Subset to residents
  subset(., state == "Disperser") %>%

  # Transform to the crs of the bmap
  spTransform(., CRS("+init=epsg:3857")) %>%

  # Plot the tracks
  plot(.
    , col = "orange"
    , add = TRUE
    , lwd = 3
  )

# Store the plot
dev.off()











################################################################################
#### NEEWWW
################################################################################
################################################################################
#### Visualization of Net Squared Displacement
################################################################################
# Description: Visualization of NSD.

# Clear R's brain
rm(list = ls())

# Load packages
library(tidyverse)
library(raster)
library(lubridate)
library(ggdark)

# Set the working directory
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_0"
setwd(wd)

# Load gps data of Abel
dat <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv() %>%
  subset(DogName == "Abel")

# Expand dispersal period a bit
first <- max(dat$Timestamp[dat$State == "Disperser"])
last <- max(dat$Timestamp[dat$State == "Disperser"]) + days(2)
dat$State[dat$Timestamp >= first & dat$Timestamp <= last] <- "Disperser"

# Identify first day of dispersal and subtract some time
first <- dat %>%
  subset(State == "Disperser") %>%
  dplyr::select(Timestamp) %>%
  head(1) %>%
  mutate(Delayed = Timestamp - days(15))

# Identify last day of dispersal and add some time
last <- dat %>%
  subset(State == "Disperser") %>%
  dplyr::select(Timestamp) %>%
  tail(1) %>%
  mutate(Delayed = Timestamp + days(15))

# Subset to these dates
dat <- subset(dat, Timestamp >= first$Delayed & Timestamp <= last$Delayed)

# Make it spatial
dat <- SpatialPointsDataFrame(
    coords      = dat[, c("x", "y")]
  , proj4string = CRS("+init=epsg:4326")
  , data        = dat
)

# Calculate NSD
dat <- spTransform(dat, CRS("+init=epsg:32734"))
dat$NSD <- (dat$x - dat$x[1])**2 + (dat$y - dat$y[1])**2 / 1000
dat <- as.data.frame(dat)

# Visualize
p <- ggplot(dat, aes(x = Timestamp, y = NSD, color = factor(State)))+
  geom_path(aes(group = 1)) +
  geom_point() +
  scale_color_manual(values = c("orange", "white"), name = "State") +
  dark_theme_classic()

# Make background transparent
p <- p +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

# Store the plot
ggsave("NSD.png"
  , plot    = p
  , width   = 14
  , height  = 8
  , scale   = 0.4
  , bg      = "transparent"
)


























# Load custom functions
source("/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/Functions.r")

############################################################
#### Plot of Africa
############################################################
# Load required shapefiles
africa    <- shapefile("00_General_Africa")
kaza      <- shapefile("00_General_KAZA_KAZA")
dogs      <- shapefile("00_General_WildDogs_IUCN")
water     <- shapefile("03_LandscapeFeatures_MajorWaters_GEOFABRIK")
prot      <- shapefile("02_LandUseTypes_Protected_PeaceParks(3Classes)")

# Prepare the extent of our core area
bbox <- as(extent(c(23.25, 24, -19.8, -19.05)), "SpatialPolygons")
crs(bbox) <- CRS("+init=epsg:4326")
bbox <- as(bbox, "SpatialLines")

# Remove small islands for plotting
africa <- subset(africa, !(ID %in% c(27:41, 689, 690)))

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Add pseudo data
dogs$Pseudo <- "Remaining Wild Dog Populations"

# Plot Africa with distribution of dogs
p1 <- tm_shape(africa) +
    tm_polygons(
        col         = viridis(20)[5]
      , border.col  = "white"
      , lwd         = 0.4
    ) +
  tm_shape(dogs) +
    tm_polygons("Pseudo"
      , palette       = viridis(20)[20]
      , border.col    = "white"
      , title         = ""
    ) +
  tm_layout(
      bg                = "transparent"
    , legend.text.color = "white"
    , frame             = FALSE
)

# Plot africa with extent of kaza
p2 <- tm_shape(africa) +
    tm_polygons(
        col         = viridis(20)[5]
      , border.col  = "white"
      , lwd         = 0.4
    ) +
  tm_shape(kaza) +
    tm_polygons(
        col           = "white"
      , border.col    = "white"
      , title         = ""
      , lwd           = 3
      , alpha         = 0.7
    ) +
  tm_shape(kaza_ext) +
    tm_borders(
        col = "white"
      , lwd = 3
    ) +
  tm_layout(
      bg                = "transparent"
    , legend.text.color = "white"
    , frame             = FALSE
)

# Plot kaza with country borders
p3 <- tm_shape(water) +
    tm_polygons(
        col         = viridis(20)[5]
      , border.col  = viridis(20)[5]
    ) +
  tm_shape(kaza, is.master = TRUE) +
    tm_polygons(
        col         = "white"
      , border.col  = "white"
      , alpha       = 0.2
      , lwd         = 5
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "white"
      , lty = 2
      , lwd = 4
    ) +
  tm_layout(
      bg    = "transparent"
    , frame = FALSE
)

# Plot kaza with protected areas
p4 <- tm_shape(prot) +
    tm_polygons(
        col         = viridis(20)[11]
      , border.col  = viridis(20)[11]
      , lwd         = 0.2
    ) +
  tm_shape(water) +
    tm_polygons(
        col         = viridis(20)[5]
      , border.col  = viridis(20)[5]
    ) +
  tm_shape(kaza, is.master = TRUE) +
    tm_polygons(
        col         = "white"
      , border.col  = "white"
      , alpha       = 0.2
      , lwd         = 5
    ) +
  tm_layout(
      bg    = "transparent"
    , frame = FALSE
)

# Plot kaza with protected areas
p5 <- tm_shape(water) +
    tm_polygons(
        col         = viridis(20)[5]
      , border.col  = viridis(20)[5]
    ) +
  tm_shape(kaza, is.master = TRUE) +
    tm_polygons(
        col         = "white"
      , border.col  = "white"
      , alpha       = 0.2
      , lwd         = 5
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "white"
      , lty = 2
      , lwd = 4
    ) +
  tm_shape(bbox) +
    tm_lines(
        col = viridis(20)[20]
      , lwd = 2
    ) +
  tm_layout(
      bg    = "transparent"
    , frame = FALSE
)

# Store the plots
setwd(output)
png("99_Africa.png"
  , width   = 1080
  , height  = 1080
  , bg      = "transparent"
  , pointsize = 20
)
p1 + tm_layout(scale = 2)
dev.off()
setwd(input)

setwd(output)
png("99_Africa(WithKAZA).png"
  , width   = 1080
  , height  = 1080
  , bg      = "transparent"
  , pointsize = 20
)
p2 + tm_layout(scale = 2)
dev.off()
setwd(input)

setwd(output)
png("99_KAZA.png"
  , width   = 1080
  , height  = 1080
  , bg      = "transparent"
  , pointsize = 20
)
p3 + tm_layout(scale = 2)
dev.off()
setwd(input)

setwd(output)
png("99_KAZA(WithProtected).png"
  , width   = 1080
  , height  = 1080
  , bg      = "transparent"
  , pointsize = 20
)
p4 + tm_layout(scale = 2)
dev.off()
setwd(input)

setwd(output)
png("99_KAZA(WithStudyArea).png"
  , width   = 1080
  , height  = 1080
  , bg      = "transparent"
  , pointsize = 20
)
p5 + tm_layout(scale = 2)
dev.off()
setwd(input)

############################################################
#### Plotting Trajectories
############################################################
# Load tracks
tracks <- shapefile("00_General_Dispersers_Popecol(Regular).shp")

# Prepare the extent of our core area
bbox <- as(extent(c(23.25, 24, -19.8, -19.05)), "SpatialPolygons")
crs(bbox) <- CRS("+init=epsg:4326")
bbox <- as(bbox, "SpatialLines")

# Plot an empty map
setwd(output)
png("99_Trajectories1.png", width = 1920, height = 1080)

# Get the bing map
prettymap(
  bmaps.plot(bbox(tracks)
    , type = "Aerial"
    , key = "Ai0A3oHx1c3rivWQfR2nAmSle5AVyo6RxjtzUFLkxTT2_qbxyNDSge0A5jxxnDgW"
    , stoponlargerequest = FALSE
  )
  , scale.style     = "ticks"
  , scale.label.col = "white"
  , scale.linecol   = "white"
  , scale.label.cex = 3
  , scale.lwd       = 5
)

# Store the plot
dev.off()
setwd(input)

# Plot an empty map with the extent
setwd(output)
png("99_Trajectories2.png", width = 1920, height = 1080)

# Get the bing map
prettymap(
  bmaps.plot(bbox(tracks)
    , type = "Aerial"
    , key = "Ai0A3oHx1c3rivWQfR2nAmSle5AVyo6RxjtzUFLkxTT2_qbxyNDSge0A5jxxnDgW"
    , stoponlargerequest = FALSE
  )
  , scale.style     = "ticks"
  , scale.label.col = "white"
  , scale.linecol   = "white"
  , scale.label.cex = 3
  , scale.lwd       = 5
)

# Add the extent
bbox %>%
  spTransform(., CRS("+init=epsg:3857")) %>%
  plot(
      col = "white"
    , add = TRUE
    , lwd = 5
  )

# Store the plot
dev.off()
setwd(input)

# Plot all trajectories onto a bing map.
setwd(output)
png("99_Trajectories3.png", width = 1920, height = 1080)

# Get the bing map
prettymap(
  bmaps.plot(bbox(tracks)
    , type = "Aerial"
    , key = "Ai0A3oHx1c3rivWQfR2nAmSle5AVyo6RxjtzUFLkxTT2_qbxyNDSge0A5jxxnDgW"
    , stoponlargerequest = FALSE
  )
  , scale.style     = "ticks"
  , scale.label.col = "white"
  , scale.linecol   = "white"
  , scale.label.cex = 3
  , scale.lwd       = 5
)

# Add the trajectories onto the satellite map
tracks %>%

  # Aggregate them by individual
  aggregate(., by = "id") %>%

  # Transform to the crs of the bmap
  spTransform(., CRS("+init=epsg:3857")) %>%

  # Plot the tracks
  plot(.
    , col = rev(viridis(nrow(.)))
    , add = TRUE
    , lwd = 3
  )

# Add the extent
bbox %>%
  spTransform(., CRS("+init=epsg:3857")) %>%
  plot(
      col = "white"
    , add = TRUE
    , lwd = 5
  )

# Store the plot
dev.off()
setwd(input)

# We also want to plot all trajectories but with different colours for
# dispersers and residents
setwd(output)
png("99_Trajectories4.png", width = 1920, height = 1080)

# Get the bing map
prettymap(
  bmaps.plot(bbox(tracks)
    , type = "Aerial"
    , key = "Ai0A3oHx1c3rivWQfR2nAmSle5AVyo6RxjtzUFLkxTT2_qbxyNDSge0A5jxxnDgW"
    , stoponlargerequest = FALSE
  )
  , scale.style     = "ticks"
  , scale.label.col = "white"
  , scale.linecol   = "white"
  , scale.label.cex = 3
  , scale.lwd       = 5
)

# Add the extent
bbox %>%
  spTransform(., CRS("+init=epsg:3857")) %>%
  plot(
      col = "white"
    , add = TRUE
    , lwd = 5
  )

# Add the resident trajectories onto the satellite map
tracks %>%

  # Subset to residents
  subset(., state == "Resident") %>%

  # Transform to the crs of the bmap
  spTransform(., CRS("+init=epsg:3857")) %>%

  # Plot the tracks
  plot(.
    , col = viridis(20)[6]
    , add = TRUE
    , lwd = 3
  )

# Add the resident trajectories onto the satellite map
tracks %>%

  # Subset to residents
  subset(., state == "Disperser") %>%

  # Transform to the crs of the bmap
  spTransform(., CRS("+init=epsg:3857")) %>%

  # Plot the tracks
  plot(.
    , col = viridis(20)[20]
    , add = TRUE
    , lwd = 3
  )

# Store the plot
dev.off()
setwd(input)

############################################################
#### Covariates
############################################################
# Prepare a third bounding box for the extent of Maun
bbox <- extent(c(21.74909, 24.30089, -20.64901, -18.14901))

# Prepare the extent of our core area
core <- as(extent(c(23.25, 24, -19.8, -19.05)), "SpatialPolygons")
crs(core) <- CRS("+init=epsg:4326")
core <- as(core, "SpatialLines")

# Load required data
water   <- stack("01_LandCoverClasses30_Globeland(Merged).grd")
trees   <- raster("01_LandCover_TreeCover_MODIS.tif")
shrubs  <- raster("01_LandCover_NonTreeVegetation_MODIS.tif")
prot    <- raster("02_LandUseTypes_Protected_PeaceParks(1Class).tif")
humans  <- raster("04_AnthropogenicFeatures_HumanInfluence(Buffer5000).tif")
roads1  <- raster("04_AnthropogenicFeatures_DistanceToRoads.tif")
roads2  <- shapefile("04_AnthropogenicFeatures_Roads_GEOFABRIK")

# Crop all data to an extent that is slightly larger than the bbox
water   <- crop(water, extent(bbox) * 1.2) %>% calc(., fun = mean)
trees   <- crop(trees, extent(bbox) * 1.2) %>% calc(., fun = function(x){100*x})
shrubs  <- crop(shrubs, extent(bbox) * 1.2) %>% calc(., fun = function(x){100*x})
prot    <- crop(prot, extent(bbox) * 1.2)
humans  <- crop(humans, extent(bbox) * 1.2)
roads1  <- crop(roads1, extent(bbox) * 1.2) %>% calc(., fun = function(x){x/1000})
roads2  <- crop(roads2, extent(bbox) * 1.2)

# Remove outliers in the tree layer
upper <- quantile(values(trees), 0.99, na.rm = TRUE)
lower <- quantile(values(trees), 0.01, na.rm = TRUE)
values(trees)[values(trees) > upper] <- upper
values(trees)[values(trees) < lower] <- lower

# Remove outliers in the shrub layer
upper <- quantile(values(shrubs), 0.99, na.rm = TRUE)
lower <- quantile(values(shrubs), 0.01, na.rm = TRUE)
values(shrubs)[values(shrubs) > upper] <- upper
values(shrubs)[values(shrubs) < lower] <- lower

# Add some pseudo data to the core study area lines
core$Pseudo <- ""

# Plot the layers
p1 <- tm_shape(water, bbox = bbox) +
  tm_raster(
      palette = "viridis"
    , style   = "cont"
    , title   = "Inundation\nFrequency"
    , labels = c("Low", "", "High")
  ) +
  tm_layout(
      legend.position   = c("left", "bottom")
    , legend.bg.color   = "black"
    , legend.text.col   = "white"
    , legend.title.col  = "white"
    , legend.title.size = 0.8
    , legend.width      = 0.5
    , frame             = "white"
    , frame.lwd         = 3
    , legend.show       = FALSE
  ) +
  tm_compass(
      color.dark = "white"
    , text.color = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
)

p2 <- tm_shape(trees, bbox = bbox) +
  tm_raster(
      palette = "viridis"
    , style   = "cont"
    , title   = "Tree Cover (%)"
    , breaks = c(0, 7, 14)
  ) +
  tm_layout(
      legend.position   = c("left", "bottom")
    , legend.bg.color   = "black"
    , legend.text.col   = "white"
    , legend.title.col  = "white"
    , legend.title.size = 0.8
    , legend.width      = 0.5
    , frame             = "white"
    , frame.lwd         = 3
    , legend.show       = FALSE
  ) +
  tm_compass(
      color.dark = "white"
    , text.color = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
)

p3 <- tm_shape(shrubs, bbox = bbox) +
  tm_raster(
      palette = "viridis"
    , style   = "cont"
    , title   = "Shrub/Grassland\nCover (%)"
    , breaks  = c(0, 50, 100)
  ) +
  tm_layout(
      legend.position   = c("left", "bottom")
    , legend.bg.color   = "black"
    , legend.text.col   = "white"
    , legend.title.col  = "white"
    , legend.title.size = 0.8
    , legend.width      = 0.5
    , frame             = "white"
    , frame.lwd         = 3
    , legend.show       = FALSE
  ) +
  tm_compass(
      color.dark = "white"
    , text.color = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
)

p4 <- tm_shape(prot, bbox = bbox) +
  tm_raster(
    , palette = viridis(10)[c(1, 8)]
    , style   = "cat"
    , title   = "Protection\nStatus"
    , labels  = c("Unprotected", "Proctected")
  ) +
  tm_layout(
      legend.position   = c("left", "bottom")
    , legend.bg.color   = "black"
    , legend.text.col   = "white"
    , legend.title.col  = "white"
    , legend.title.size = 0.8
    , legend.width      = 0.5
    , frame             = "white"
    , frame.lwd         = 3
    , legend.show       = FALSE
  ) +
  tm_compass(
      color.dark = "white"
    , text.color = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
)

p5 <- tm_shape(humans, bbox = bbox) +
  tm_raster(
      palette = "viridis"
    , style   = "cont"
    , title   = "Human\nInfluence"
    , breaks  = c(0, 6, 12)
  ) +
  tm_layout(
      legend.position   = c("left", "bottom")
    , legend.bg.color   = "black"
    , legend.text.col   = "white"
    , legend.title.col  = "white"
    , legend.title.size = 0.8
    , legend.width      = 0.5
    , frame             = "white"
    , frame.lwd         = 3
    , legend.show       = FALSE
  ) +
  tm_compass(
      color.dark = "white"
    , text.color = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
)

p6 <- tm_shape(roads1, bbox = bbox) +
  tm_raster(
      palette = "viridis"
    , style   = "cont"
    , title   = "Distance to\nRoads (km)"
    , breaks  = c(0, 60, 120)
  ) +
  tm_shape(roads2, bbox = bbox) +
    tm_lines(col = "white") +
  tm_layout(
      legend.position   = c("left", "bottom")
    , legend.bg.color   = "black"
    , legend.text.col   = "white"
    , legend.title.col  = "white"
    , legend.title.size = 0.8
    , legend.width      = 0.5
    , frame             = "white"
    , frame.lwd         = 3
    , legend.show       = FALSE
  ) +
  tm_compass(
      color.dark = "white"
    , text.color = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
)


# Put all plots together
p1 <- p1 + tm_layout(scale = 2)
p2 <- p2 + tm_layout(scale = 2)
p3 <- p3 + tm_layout(scale = 2)
p4 <- p4 + tm_layout(scale = 2)
p5 <- p5 + tm_layout(scale = 2)
p6 <- p6 + tm_layout(scale = 2)

p <- tmap_arrange(p1, p2, p3, p4, p5, p6, ncol = 3, asp = 1.1)

setwd(output)
png("99_Covariates.png"
  , width   = 1080
  , height  = 1080 * 3.5 / 6
  , bg      = "transparent"
)
p
dev.off()
setwd(input)

############################################################
#### How Least Cost Paths and Corridors Work
############################################################
# Specify size of raster
n <- 100

# Simulate Gaussian Field
map1 <- nlm_gaussianfield(ncol = n, nrow = n)

# Simulate an area that is avoided in the center
map2 <- nlm_distancegradient(ncol = n, nrow = n, origin = c(n/2, n/2, n/2, n/2))

# Put the two maps together
map <- map1 + map2

# Write function to normalize a map
norm <- function(x){
  (x - minValue(x))/(maxValue(x) - minValue(x))
}

# Normalize map
map <- norm(map)

# Let's use this as our "premeability" map, where larger values indicate easier
# permeability. We now convert the raster to a transition layer
trans <- transition(map, transitionFunction = mean, directions = 16)

# Apply correction
trans <- geoCorrection(trans, type = "c", multpl = FALSE)

# Create points to connect
A <- data.frame(x = 5, y = 50) %>% SpatialPoints()
B <- data.frame(x = 95, y = 50) %>% SpatialPoints()
points <- rbind(A, B)
points$ID <- c("A", "B")

# Find least cost path
path <- shortestPath(trans
  , origin  = points[1, ]
  , goal    = points[2, ]
  , output  = "SpatialLines"
)

# We can also calculate the cost map from the two locations
A_cost <- accCost(trans, points[1, ])
B_cost <- accCost(trans, points[2, ])

# Put the maps together
cost <- A_cost + B_cost

# Identify a desired percentile
quant <- quantile(cost, probs = 0.1, na.rm = TRUE)

# Define new map
corr <- cost

# Reclassify values above the quantile to get
values(corr) <- ifelse(values(corr) > quant, NA, values(corr))

# Add pseudo data
path$Pseudo <- "Least-Cost Path"

# Prepare Plot of source points
p1 <- tm_shape(norm(map)) +
    tm_raster(
        palette = "viridis"
      , style   = "cont"
      , title   = "Permeability"
      , labels  = c("Low", "", "High")
      , legend.show = FALSE
    ) +
  tm_shape(points) +
    tm_dots(
        col = "white"
      , size = 4
    ) +
    tm_text("ID", size = 1.5) +
  tm_layout(
      legend.title.color  = "white"
    , legend.text.color   = "white"
    , legend.position     = c("left", "bottom")
)

# Plot of LCPs
p2 <- tm_shape(norm(map)) +
    tm_raster(
        palette = "viridis"
      , style   = "cont"
      , title   = "Permeability"
      , labels  = c("Low", "", "High")
      , legend.show = FALSE
    ) +
  tm_shape(path, legend.show = FALSE) +
    tm_lines("Pseudo"
      , palette         = c("white")
      , lwd             = 4
      , legend.col.show = FALSE
    ) +
  tm_shape(points) +
    tm_dots(
        col   = "white"
      , size  = 4
    ) +
    tm_text("ID", size = 1.5) +
  tm_layout(
      legend.title.color  = "white"
    , legend.text.color   = "white"
    , legend.position     = c("left", "bottom")
)

# Plot the cost maps
p3 <- tm_shape(norm(A_cost)) +
    tm_raster(
        palette = "-viridis"
      , style   = "cont"
      , title   = "Cumulative Cost"
      , labels  = c("Low", "", "High")
      , legend.show = FALSE
    ) +
  tm_shape(points[1, ]) +
    tm_dots(
        col = "white"
      , size = 4
    ) +
    tm_text("ID", size = 1.5) +
  tm_layout(
      legend.title.color  = "white"
    , legend.text.color   = "white"
    , legend.position     = c("left", "bottom")
)

p4 <- tm_shape(norm(B_cost)) +
    tm_raster(
        palette = "-viridis"
      , style   = "cont"
      , title   = "Cumulative Cost"
      , labels  = c("Low", "", "High")
      , legend.show = FALSE
    ) +
  tm_shape(points[2, ]) +
    tm_dots(
        col = "white"
      , size = 4
    ) +
    tm_text("ID", size = 1.5) +
  tm_layout(
      legend.title.color  = "white"
    , legend.text.color   = "white"
    , legend.position     = c("left", "bottom")
)

p5 <- tm_shape(norm(cost)) +
    tm_raster(
        palette = "-viridis"
      , style   = "cont"
      , title   = "Cumulative Cost"
      , labels  = c("Low", "", "High")
      , legend.show = FALSE
    ) +
  tm_shape(points) +
    tm_dots(
        col = "white"
      , size = 4
    ) +
    tm_text("ID", size = 1.5) +
  tm_layout(
      legend.title.color  = "white"
    , legend.text.color   = "white"
    , legend.position     = c("left", "bottom")
)

p6 <- tm_shape(norm(corr)) +
    tm_raster(
        palette = "-viridis"
      , style   = "cont"
      , title   = "Corridor"
      , labels  = c("Low-Cost", "", "High-Cost")
      , legend.show = FALSE
    ) +
  tm_shape(points) +
    tm_dots(
        col = "white"
      , size = 4
    ) +
    tm_text("ID", size = 1.5) +
  tm_layout(
      legend.title.color  = "white"
    , legend.text.color   = "white"
    , legend.position     = c("left", "bottom")
    , bg.color            = "black"
)

# Store the plots
setwd(output)
png("99_LeastCostPaths(Example1).png"
  , width   = 720
  , height  = 720
  , bg      = "transparent"
  , pointsize = 30
)
p1
dev.off()
setwd(input)

setwd(output)
png("99_LeastCostPaths(Example2).png"
  , width   = 720
  , height  = 720
  , bg      = "transparent"
  , pointsize = 30
)
p2
dev.off()
setwd(input)

setwd(output)
png("99_LeastCostPaths(Example3).png"
  , width   = 720
  , height  = 720
  , bg      = "transparent"
  , pointsize = 30
)
p3
dev.off()
setwd(input)

setwd(output)
png("99_LeastCostPaths(Example4).png"
  , width   = 720
  , height  = 720
  , bg      = "transparent"
  , pointsize = 30
)
p4
dev.off()
setwd(input)

setwd(output)
png("99_LeastCostPaths(Example5).png"
  , width   = 720
  , height  = 720
  , bg      = "transparent"
  , pointsize = 30
)
p5
dev.off()
setwd(input)

setwd(output)
png("99_LeastCostPaths(Example6).png"
  , width   = 720
  , height  = 720
  , bg      = "transparent"
  , pointsize = 30
)
p6
dev.off()
setwd(input)

############################################################
#### How We Selected Source Points
############################################################
# Load required data
perm        <- raster("99_PermeabilityMap.tif")
kaza        <- shapefile("00_General_KAZA_KAZA") %>% as(., "SpatialLines")
africa      <- shapefile("00_General_Africa") %>% as(., "SpatialLines")
africa_crop <- shapefile("00_General_Africa") %>% crop(., kaza)
prot1       <- shapefile("02_LandUseTypes_Protected_PeaceParks(1Class)")
prot2       <- shapefile("99_SourceAreas")
points3     <- shapefile("99_SourcePoints")

# We only keep the countries of interest in the cropped africa file
africa_crop <- subset(africa_crop, COUNTRY %in% c(
    "Angola"
  , "Namibia"
  , "Botswana"
  , "Zimbabwe"
  , "Zambia")
)

# Create a regular grid of points (not subsetted to protected areas yet)
points1 <- perm %>%
  aggregate(., fact = (100 * 1000) / 250, fun = mean) %>%
  rasterToPoints(., spatial = TRUE)

# Subset to points within protected areas
index <- gContains(prot2, points1, byid = TRUE) %>% rowSums(.)
points2 <- points1[index > 0, ]

# Prepare a plot of the protected areas
p1 <- tm_shape(perm) +
    tm_raster(
      palette     = "black"
    , legend.show = FALSE
  ) +
  tm_shape(prot1) +
    tm_polygons(
      col           = viridis(20)[5]
    , border.col    = viridis(20)[5]
  ) +
  tm_shape(kaza) +
    tm_lines(
        col = "white"
      , lwd = 2
    ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_layout(
      bg.color            = "black"
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = "black"
    , legend.position     = c("left", "top")
)

# Prepare a plot of the protected areas
p2 <- tm_shape(perm) +
    tm_raster(
      palette     = "black"
    , legend.show = FALSE
  ) +
  tm_shape(prot2) +
    tm_polygons(
      col           = viridis(20)[5]
    , border.col    = viridis(20)[5]
  ) +
  tm_shape(kaza) +
    tm_lines(
        col = "white"
      , lwd = 2
    ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_layout(
      bg.color            = "black"
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = "black"
    , legend.position     = c("left", "top")
)

# Prepare a plot of the protected areas with all points
p3 <- tm_shape(perm) +
    tm_raster(
      palette     = "black"
    , legend.show = FALSE
  ) +
  tm_shape(prot2) +
    tm_polygons(
      col           = viridis(20)[5]
    , border.col    = viridis(20)[5]
  ) +
  tm_shape(kaza) +
    tm_lines(
        col = "white"
      , lwd = 2
    ) +
  tm_shape(points1) +
    tm_dots(
        col   = viridis(20)[20]
      , size  = 0.1
    ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_layout(
      bg.color            = "black"
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = "black"
    , legend.position     = c("left", "top")
)

# Prepare a plot of the protected areas with points inside protected
p4 <- tm_shape(perm) +
    tm_raster(
      palette     = "black"
    , legend.show = FALSE
  ) +
  tm_shape(prot2) +
    tm_polygons(
      col           = viridis(20)[5]
    , border.col    = viridis(20)[5]
  ) +
  tm_shape(kaza) +
    tm_lines(
        col = "white"
      , lwd = 2
    ) +
  tm_shape(points2) +
    tm_dots(
        col   = viridis(20)[20]
      , size  = 0.1
    ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_layout(
      bg.color            = "black"
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = "black"
    , legend.position     = c("left", "top")
)

# Prepare a plot of the protected areas with points inside protected + centroids
p5 <- tm_shape(perm) +
    tm_raster(
      palette     = "black"
    , legend.show = FALSE
  ) +
  tm_shape(prot2) +
    tm_polygons(
      col           = viridis(20)[5]
    , border.col    = viridis(20)[5]
  ) +
  tm_shape(kaza) +
    tm_lines(
        col = "white"
      , lwd = 2
    ) +
  tm_shape(points3) +
    tm_dots(
        col   = viridis(20)[20]
      , size  = 0.1
    ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_layout(
      bg.color            = "black"
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = "black"
    , legend.position     = c("left", "top")
)

# Store the plots
setwd(output)
png("99_SourcePoints1.png"
  , width   = 1080 * 4 / 3.5
  , height  = 1080
  , bg = "transparent"
)
p1 + tm_layout(scale = 3)
dev.off()
setwd(input)

setwd(output)
png("99_SourcePoints2.png"
  , width   = 1080 * 4 / 3.5
  , height  = 1080
  , bg = "transparent"
)
p2 + tm_layout(scale = 3)
dev.off()
setwd(input)

setwd(output)
png("99_SourcePoints3.png"
  , width   = 1080 * 4 / 3.5
  , height  = 1080
  , bg = "transparent"
)
p3 + tm_layout(scale = 3)
dev.off()
setwd(input)

setwd(output)
png("99_SourcePoints4.png"
  , width   = 1080 * 4 / 3.5
  , height  = 1080
  , bg = "transparent"
)
p4 + tm_layout(scale = 3)
dev.off()
setwd(input)

setwd(output)
png("99_SourcePoints5.png"
  , width   = 1080 * 4 / 3.5
  , height  = 1080
  , bg = "transparent"
)
p5 + tm_layout(scale = 3)
dev.off()
setwd(input)

############################################################
#### The ORI Algorithm
############################################################
# Load the required data
setwd("/home/david/Downloads/MODIS/MCD43A4/Masks")
wet_mask  <- shapefile("MaskWater.shp")
dry_mask  <- shapefile("MaskDryland.shp")
setwd(input)
kaza      <- shapefile("00_General_KAZA_KAZA")
water     <- shapefile("03_LandscapeFeatures_MajorWaters_GEOFABRIK")
africa    <- shapefile("00_General_Africa")

# Also load one of the classified MODIS layers
setwd("/home/david/Downloads/MODIS/MCD43A4/Output")
flood     <- raster(dir(pattern = ".tif$")[500], bands = 1)
values(flood)[values(flood) == 255] <- 2

# Crop the africa layer to the kaza polygon
africa <- crop(africa, kaza)

# Put the masks together into one dataframe
wet_mask$LandCover <- "Water"
dry_mask$LandCover <- "Dryland"
masks <- rbind(wet_mask, dry_mask)

# Define a bounding box for which we calculated dynamic floodmaps
bbox <- extent(c(21.74909, 24.30089, -20.64901, -18.14901))
bbox <- as(bbox, "SpatialPolygons")
crs(bbox) <- CRS("+init=epsg:4326")
bbox_ext <- bb(bbox)

# Download Satellite Image as background map
map <- read_osm(bbox_ext
  , type                = "bing"
  , current.projection  = CRS("+init=epsg:4326")
)

# Plot the masks on top of the satellite image
p1 <- tm_shape(map) +
  tm_rgb() +
  tm_shape(masks) +
    tm_polygons("LandCover"
      , palette       = c(viridis(20)[c(16, 5)])
      , border.alpha  = 0
      , title         = "Polygons"
    ) +
  tm_compass(
      color.dark = "white"
    , text.color = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
  ) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , projection          = "+init=epsg:4326"
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
      legend.position = c("right", "top")
    , legend.bg.color = "white"
  ) +
  tm_credits("(a)"
  , col       = "white"
  , position  = c("left", "bottom")
  , size      = 1.5
)

# Plot the classified map
p2 <- tm_shape(map) +
  tm_rgb() +
  tm_shape(flood) +
    tm_raster(
      palette         = viridis(20, alpha = 0.4)[c(5, 16)]
    , title           = "Classification"
    , labels          = c("Water", "Dryland")
    , style           = "cat"
    , legend.reverse  = TRUE
    ) +
  tm_compass(
      color.dark = "white"
    , text.color = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
  ) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , projection          = "+init=epsg:4326"
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
      legend.position = c("right", "top")
    , legend.bg.color = "white"
  ) +
  tm_credits("(b)"
  , col       = "white"
  , position  = c("left", "bottom")
  , size      = 1.5
)

# Put the two plots together
p3 <- tmap_arrange(p1, p2)

# Prepare a plot showing the floodmapping extent
p4 <- tm_shape(kaza) +
    tm_polygons(col = "black", border.alpha = 0) +
  tm_shape(subset(water)[1, ]) +
    tm_polygons(col = "white", border.alpha  = 0) +
  tm_shape(bbox) +
    tm_borders(col = "white") +
  tm_layout(bg.color = "transparent", frame = "transparent")

# Store the plot
setwd(output)
CairoPDF("99_Floodmapping.pdf", width = 9, height = 4.24)
p3
# print(p4, vp = viewport(0.085, 0.13, width = 0.2, height = 0.2))
dev.off()
setwd(input)

############################################################
#### Permeability Model Results
############################################################
# Load required data
model <- read_rds("99_ModelSelection.rds")

# Select best model
best <- model$Model[[1]]

# Summary of the best model
summary(best)

# Prepare and store a plot of the coefficients
coeffs <- getCoeffs(best)[-1, ]
coeffs$Covariate[coeffs$Covariate == "Shrubs"] <- "Shrubs/Grassland"

# Prepare plot
p <- ggplot(data = coeffs, aes(y = Covariate, x = Coefficient)) +

  # Add points for the estimated coefficients
  geom_point(shape = 1, size = 1, col = "white") +

  # Add errorbars
  geom_errorbarh(aes(
      xmin = Coefficient - 1.96 * SE
    , xmax = Coefficient + 1.96 * SE
    , height = 0.2)
    , col = "white"
  ) +

  # Add a vertical line that indicates no effect
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +

  # Reverse the scale of the y-axis
  scale_y_discrete(limits = rev(c(
      "cos(ta_)"
    , "log(sl_)"
    , "Water"
    , "DistanceToWater"
    , "Trees"
    , "Shrubs/Grassland"
    , "HumansBuff5000"
  ))) +

  # Use a simple black and white theme
  theme_minimal() +

  # Scale axis to -1, 1
  xlim(c(-1, 1)) +

  # Change label of x-axis
  labs(x = expression(beta*"-Coefficient")) +

  # Cap the axes
  coord_capped_cart(left = "both", bottom = "both") +

  theme(
    plot.background = element_rect(fill = "transparent")
  , axis.text       = element_text(color = "white")
  , axis.title      = element_text(color = "white")
  , axis.ticks      = element_line(color = "white")
  , axis.line       = element_line(color = "white")
  , panel.grid      = element_line(color = "transparent")
)

# Write the plot to file
setwd(output)
ggsave("99_PermeabilityModel.png"
  , plot    = p
  , width   = 14
  , height  = 12
  , scale   = 0.25
  , bg      = "transparent"
)
dev.off()
setwd(input)

############################################################
#### Permeability Model Validation
############################################################
# Load required data
validation <- read_rds("99_ModelValidation.rds")
dat_pref <- read_rds("99_ModelValidation(Data).rds")[[1]]
dat_rand <- read_rds("99_ModelValidation(Data).rds")[[2]]

# Put the k-fold cross validation data from observed and random preferences
# together
dat_pref$Group <- "Realized"
dat_rand$Group <- "Random"
dat <- rbind(dat_pref, dat_rand)

# We want to plot this data and add the information from the validation table
# too. Let's prepare a column that indicates the text that we want to plot on
# top of the data. Let's first round the values from the validation table
validation[, c(2:4)] <- round(validation[, c(2:4)], 2)

# Match the Grouping names of the validation table to the groups above
validation$Group <- as.factor(c("Realized", "Random"))

# Create a dataframe which we use to annotate the two facets
text <- data.frame(
    Group = c("Realized", "Random")
  , Text = c(
      "$\\bar{r}_s = -0.45$, $95%-CI = -0.48$, $-0.42$"
    , "$\\bar{r}_s = 0.04$, $95%-CI = 0.01$, $0.08$"
  )
)

# Reorder the factors in the Group variable
dat$Group <- factor(dat$Group, levels = c("Realized", "Random"))

# Plot the data
p1 <- ggplot(data = dat, aes(x = Rank, y = Frequency)) +
  geom_jitter(alpha = 0.2, size = 0.5) +
  geom_smooth(method = "loess") +
  theme_classic() +
  facet_wrap("Group") +
  geom_text(data = text
    , mapping = aes(
        x = -Inf
      , y = -Inf
      , label = TeX(Text, output = "character")
    )
    , hjust   = -0.05
    , vjust   = -0.5
    , parse   = TRUE
    , size    = 3
)

# Write the plot to file
setwd(output)
ggsave("99_PermeabilityValidation.png"
  , plot    = p1
  , width   = 16
  , height  = 9
  , scale   = 0.5
)
p1
dev.off()
setwd(input)

############################################################
#### Permeability Surface
############################################################
# Load required data
permeability  <- raster("99_PermeabilityMap.tif")
kaza          <- shapefile("00_General_KAZA_KAZA") %>% as(., "SpatialLines")
africa        <- shapefile("00_General_Africa") %>% as(., "SpatialLines")
africa_crop   <- shapefile("00_General_Africa") %>% crop(., kaza)
tracks        <- shapefile("00_General_Dispersers_Popecol(Regular)")

# Subset to dispersers only
tracks <- subset(tracks, state == "Disperser")

# Select only countries of interest
africa_crop <- subset(africa_crop, COUNTRY %in% c(
    "Angola"
  , "Namibia"
  , "Botswana"
  , "Zimbabwe"
  , "Zambia")
)

# Prepare the extent of our core area
core <- as(extent(c(23.25, 24, -19.8, -19.05)), "SpatialLines")
crs(core) <- CRS("+init=epsg:4326")

# Rescale between 0 and 1 (necessary to get nice names)
permeability <- calc(permeability, fun = function(x){
  (x - min(x)) / (max(x) - min(x))
})

# Add some pseudodata to the kaza shapefile
kaza$Pseudo <- ""

# Create a boundary box for a smaller extent that we want to plot. We want to
# have the same width/height ratio as in the permeability map
extent <- extent(permeability)
ratio <- (extent@xmax - extent@xmin) / (extent@ymax - extent@ymin)
xmin <- 22.0
xmax <- 24.5
ymin <- -20.2
ymax <- ymin + (xmax - xmin) / ratio
bbox <- as(extent(xmin, xmax, ymin, ymax), "SpatialLines")
crs(bbox) <- CRS("+init=epsg:4326")

# Prepare a dataframe with some locations that we want to plot
villages <- data.frame(
    x = c(23.419332, 24.178607, 24.271072, 23.653653)
  , y = c(-19.999149, -19.430306, -19.192169, -19.002809)
  , Name = c("Maun", "Sankuyo", "Mababe", "Khwai")
)
coordinates(villages) <- c("x", "y")
crs(villages) <- CRS("+init=epsg:4326")

waters <- data.frame(
    x = c(22.668000, 23.611706)
  , y = c(-19.106023, -18.270367)
  , Name = c("Okavango\nDelta", "Linyanti\nSwamp")
)
coordinates(waters) <- c("x", "y")
crs(waters) <- CRS("+init=epsg:4326")

# Prepare plot with permeability only
p1 <- tm_shape(permeability) +
    tm_raster(
        palette     = "viridis"
      , style       = "cont"
      , title       = "Permeability"
      , labels      = c("Low", "", "High")
    ) +
  tm_shape(kaza) +
    tm_lines(
        col = "white"
      , lwd = 2
    ) +
  tm_shape(africa) +
    tm_lines(
        col = "white"
      , lwd = 1
      , lty = 2
    ) +
  tm_shape(africa_crop) +
    tm_text("COUNTRY"
      , col       = "white"
      , just      = "bottom"
      , fontface  = 2
    ) +
  # tm_shape(bbox) +
  #   tm_lines(
  #       col = "black"
  #     , lwd = 2
  #     ) +
  tm_layout(
        legend.position     = c("left", "top")
      , legend.text.color   = "white"
      , legend.title.color  = "white"
      , legend.bg.color     = viridis(1, begin = 0)
  ) +
  tm_scale_bar(
        position    = "left"
      , text.size   = 0.5
      , text.color  = "white"
      , width = 0.125
  ) +
  tm_compass(
        text.color = "white"
      , color.dark = "white"
      , text.size  = 1.5
      , size       = 2
)

# Prepare plot with extent of central region
p2 <- tm_shape(permeability) +
    tm_raster(
        palette     = "viridis"
      , style       = "cont"
      , title       = "Permeability"
      , labels      = c("Low", "", "High")
    ) +
  tm_shape(kaza) +
    tm_lines(
        col = "white"
      , lwd = 2
    ) +
  tm_shape(africa) +
    tm_lines(
        col = "white"
      , lwd = 1
      , lty = 2
    ) +
  tm_shape(africa_crop) +
    tm_text("COUNTRY"
      , col       = "white"
      , just      = "bottom"
      , fontface  = 2
    ) +
  tm_shape(bbox) +
    tm_lines(
        col = "black"
      , lwd = 2
      ) +
  tm_layout(
        legend.position     = c("left", "top")
      , legend.text.color   = "white"
      , legend.title.color  = "white"
      , legend.bg.color     = viridis(1, begin = 0)
  ) +
  tm_scale_bar(
        position    = "left"
      , text.size   = 0.5
      , text.color  = "white"
      , width = 0.125
  ) +
  tm_compass(
        text.color = "white"
      , color.dark = "white"
      , text.size  = 1.5
      , size       = 2
)

# Prepare a plot of the Okavango Delta only
p3 <- tm_shape(permeability, bbox = bbox) +
    tm_raster(
        palette     = "viridis"
      , style       = "cont"
      , title       = "Permeability"
      , labels      = c("Low", "", "High")
    ) +
  tm_shape(villages) +
    tm_text("Name"
    , col   = "white"
    , size  = 0.9
  ) +
  tm_shape(waters) +
    tm_text("Name"
    , col       = "white"
    , fontface  = 3
  ) +
  tm_layout(
        legend.position     = c("left", "top")
      , legend.text.color   = "white"
      , legend.title.color  = "white"
      , legend.bg.color     = viridis(1, begin = 0)
  ) +
  tm_scale_bar(
        position    = "left"
      , text.size   = 0.5
      , text.color  = "white"
      , width = 0.125
  ) +
  tm_compass(
        text.color = "white"
      , color.dark = "white"
      , text.size  = 1.5
      , size       = 2
)

# Prepare a plot of the Okavango Delta only
p4 <- tm_shape(permeability, bbox = bbox) +
    tm_raster(
        palette     = "viridis"
      , style       = "cont"
      , title       = "Permeability"
      , labels      = c("Low", "", "High")
    ) +
  tm_shape(tracks) +
    tm_lines(
        col = "white"
      , lwd = 0.7
    ) +
  tm_shape(villages) +
    tm_text("Name"
    , col   = "white"
    , size  = 0.9
  ) +
  tm_shape(waters) +
    tm_text("Name"
    , col       = "white"
    , fontface  = 3
  ) +
  tm_layout(
        legend.position     = c("left", "top")
      , legend.text.color   = "white"
      , legend.title.color  = "white"
      , legend.bg.color     = viridis(1, begin = 0)
  ) +
  tm_scale_bar(
        position    = "left"
      , text.size   = 0.5
      , text.color  = "white"
      , width = 0.125
  ) +
  tm_compass(
        text.color = "white"
      , color.dark = "white"
      , text.size  = 1.5
      , size       = 2
)

# Store the plot
setwd(output)
png("99_PermeabilityMap1.png"
  , width   = 1080 * 4 / 3.5
  , height  = 1080
  , bg      = "transparent"
)
p1 + tm_layout(scale = 3)
dev.off()
setwd(input)

setwd(output)
png("99_PermeabilityMap2.png"
  , width   = 1080 * 4 / 3.5
  , height  = 1080
  , bg      = "transparent"
)
p2 + tm_layout(scale = 3)
dev.off()
setwd(input)

setwd(output)
png("99_PermeabilityMap3.png"
  , width   = 1080 * 4 / 3.5
  , height  = 1080
  , bg      = "transparent"
)
p3 + tm_layout(scale = 3)
dev.off()
setwd(input)

setwd(output)
png("99_PermeabilityMap4.png"
  , width   = 1080 * 4 / 3.5
  , height  = 1080
  , bg      = "transparent"
)
p4 + tm_layout(scale = 3)
dev.off()
setwd(input)

############################################################
#### Least Cost Paths & Corridors
############################################################
# Load required data
perm        <- raster("99_PermeabilityMap.tif")
corrs       <- raster("99_LeastCostCorridors.tif")
paths       <- shapefile("99_LeastCostPaths")
paths_rast  <- raster("99_LeastCostPaths.tif")
kaza        <- shapefile("00_General_KAZA_KAZA") %>% as(., "SpatialLines")
africa      <- shapefile("00_General_Africa") %>% as(., "SpatialLines")
africa_crop <- shapefile("00_General_Africa") %>% crop(., kaza)
prot        <- shapefile("99_SourceAreas.shp")
points      <- shapefile("99_SourcePoints.shp")
nati        <- shapefile("02_LandUseTypes_Protected_PeaceParks(3Classes)")

# Subset to national parks
nati <- subset(nati, Desig == "National Park")

# Subset to national parks that we want to plot
nati <- subset(nati, Name %in% c("Mavinga", "Luengue-Luiana", "Kafue"
  , "Hwange", "Central Kalahari", "Chobe", "Moremi", "Matusadona", "Khaudum"))

# There is a double entry for Kafue, get rid of the erronous one
nati$Area <- gArea(nati, byid = TRUE)
nati <- subset(nati, Area != min(Area))

# Create a separate shapefile for the text. We have to change some of the
# coordinates to make sure that they don't overlap
nati_text <- nati
nati_text$x <- coordinates(nati_text)[, 1]
nati_text$y <- coordinates(nati_text)[, 2]
nati_text <- nati_text@data
nati_text$y[nati_text$Name == "Kafue"] <-
  nati_text$y[nati_text$Name == "Kafue"] + 0.5
nati_text$y[nati_text$Name == "Chobe"] <-
  nati_text$y[nati_text$Name == "Chobe"] - 0.1
nati_text$y[nati_text$Name == "Matusadona"] <-
  nati_text$y[nati_text$Name == "Matusadona"] - 0.1
coordinates(nati_text) <- c("x", "y")
crs(nati_text) <- CRS("+init=epsg:4326")

# Check how they align
plot(nati)
points(nati_text)

# Add "NP" to the text (on a new line)
head(nati_text)
nati_text$Name <- paste0(nati_text$Name, "\nNP")

# We only keep the countries of interest in the cropped africa file
africa_crop <- subset(africa_crop, COUNTRY %in% c(
    "Angola"
  , "Namibia"
  , "Botswana"
  , "Zimbabwe"
  , "Zambia")
)

# Add some pseudodata to the paths
paths$Pseudo <- ""

# Take sqrt of LCPs to make lower frequented LPCs more visible
paths_rast <- sqrt(paths_rast)

# Rescale between 0 and 1
corrs <- calc(corrs, fun = function(x){
  (x - min(x)) / (max(x) - min(x))
})
paths_rast <- calc(paths_rast, fun = function(x){
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
})

# Prepare a plot of the LCPs only
p1 <- tm_shape(corrs) +
    tm_raster(
      palette     = "black"
    , legend.show = FALSE
  ) +
  tm_shape(prot) +
    tm_polygons(
      col           = viridis(20)[5]
    , border.col    = viridis(20)[5]
  ) +
  tm_shape(paths_rast) +
    tm_raster(
        palette = viridis(20)[c(7:20)]
      , style   = "cont"
      , title   = "Least-Cost Paths"
      , labels  = c("Low-Frequency", "", "High-Frequency")
    ) +
  tm_shape(points) +
    tm_dots(
        col   = viridis(20)[20]
      , size  = 0.1
    ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_layout(
      bg.color            = "black"
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = "black"
    , legend.position     = c("left", "top")
)

# Prepare a plot of the LCPs with country and kaza border
p2 <- tm_shape(corrs) +
    tm_raster(
      palette     = "black"
    , legend.show = FALSE
  ) +
  tm_shape(prot) +
    tm_polygons(
      col           = viridis(20)[5]
    , border.col    = viridis(20)[5]
  ) +
  tm_shape(paths_rast) +
    tm_raster(
        palette = viridis(20)[c(7:20)]
      , style   = "cont"
      , title   = "Least-Cost Paths"
      , labels  = c("Low-Frequency", "", "High-Frequency")
    ) +
  tm_shape(points) +
    tm_dots(
        col   = viridis(20)[20]
      , size  = 0.1
    ) +
  tm_shape(kaza) +
    tm_lines(
        col = "white"
      , lwd = 2
    ) +
  tm_shape(africa) +
    tm_lines(
        col = "white"
      , lwd = 1
      , lty = 2
    ) +
  tm_shape(africa_crop) +
    tm_text("COUNTRY"
      , col   = "white"
      , just  = "bottom"
    ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_layout(
      bg.color            = "black"
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = "black"
    , legend.position     = c("left", "top")
)

# Prepare the plot for LCCs only
p3 <- tm_shape(corrs) +
    tm_raster(
        palette = "viridis"
      , style   = "cont"
      , title   = "Least-Cost Corridors"
      , labels = c("Low-Frequency", "", "High-Frequency")
    ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_layout(
      legend.position     = c("left", "top")
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = viridis(1, begin = 0)
)

# Prepare plot of LCCs with country borders and kaza
p4 <- tm_shape(corrs) +
    tm_raster(
        palette = "viridis"
      , style   = "cont"
      , title   = "Least-Cost Corridors"
      , labels = c("Low-Frequency", "", "High-Frequency")
    ) +
  tm_shape(nati) +
    tm_borders(
        col   = "white"
      , alpha = 0.2
    ) +
  tm_shape(kaza) +
    tm_lines(
        col = "white"
      , lwd = 2
    ) +
  tm_shape(nati_text) +
    tm_text("Name"
      , col       = "white"
      , alpha     = 0.5
      , fontface  = 3
      , size      = 0.5
    ) +
  tm_shape(africa) +
    tm_lines(
        col     = "white"
      , lwd     = 1
      , lty     = 2
    ) +
  tm_shape(africa_crop) +
    tm_text("COUNTRY"
      , col     = "white"
      , just    = "bottom"
      , shadow  = TRUE
    ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_layout(
      legend.position     = c("left", "top")
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = viridis(1, begin = 0)
)

# Prepare plot of LCCs with country borders
p5 <- tm_shape(corrs) +
    tm_raster(
        palette = "viridis"
      , style   = "cont"
      , title   = "Least-Cost Corridors"
      , labels = c("Low-Frequency", "", "High-Frequency")
    ) +
  tm_shape(kaza) +
    tm_lines(
        col = "white"
      , lwd = 2
    ) +
  tm_layout(
      legend.position     = c("left", "top")
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = viridis(1, begin = 0)
    , legend.show         = FALSE
)

# Store the plots
setwd(output)
png("99_LeastCostPaths1.png"
  , width = 1080 * 4 / 3.5
  , height = 1080
  , bg = "transparent"
)
p1 + tm_layout(scale = 3)
dev.off()
setwd(input)

setwd(output)
png("99_LeastCostPaths2.png"
  , width = 1080 * 4 / 3.5
  , height = 1080
  , bg = "transparent"
)
p2 + tm_layout(scale = 3)
dev.off()
setwd(input)

# Store the plot
setwd(output)
png("99_LeastCostCorrs1.png"
  , width = 1080 * 4 / 3.5
  , height = 1080
  , bg = "transparent"
)
p3 + tm_layout(scale = 3)
dev.off()
setwd(input)

setwd(output)
png("99_LeastCostCorrs2.png"
  , width = 1080 * 4 / 3.5
  , height = 1080
  , bg = "transparent"
)
p4 + tm_layout(scale = 3)
dev.off()
setwd(input)

setwd(output)
png("99_LeastCostCorrs3.png"
  , width = 1080 * 4 / 3.5
  , height = 1080
  , bg = "transparent"
)
p5 + tm_layout(scale = 3)
dev.off()
setwd(input)

############################################################
#### Simulations
############################################################
# Load required data
paths       <- raster("99_DispersalSimulation.tif")
perm        <- raster("99_PermeabilityMap.tif") %>% crop(., paths)
kaza        <- shapefile("00_General_KAZA_KAZA") %>% as(., "SpatialLines")
africa      <- shapefile("00_General_Africa") %>% as(., "SpatialLines")
africa_crop <- shapefile("00_General_Africa") %>% crop(., paths)
tracks      <- shapefile("00_General_Dispersers_Popecol(Regular).shp")

# Subset tracks
tracks <- subset(tracks, state == "Disperser")

# Prepare the extent of our core area
bbox <- as(extent(c(23.25, 24, -19.8, -19.05)), "SpatialPolygons")
crs(bbox) <- CRS("+init=epsg:4326")
bbox <- as(bbox, "SpatialLines")

# Take sqrt for better visibility
paths <- sqrt(paths)

# Prepare a background
bg <- paths
values(bg) <- 1

# Plot simulated paths
p1 <- tm_shape(bg) +
    tm_raster(
        palette     = viridis(20)[1]
      , legend.show = FALSE
    ) +
  tm_shape(africa) +
    tm_lines(
        col = "white"
      , lwd = 1
      , lty = 2
    ) +
  tm_shape(bbox) +
    tm_lines(
        col = "white"
      , lwd = 1
      , lty = 1
    ) +
  tm_shape(africa_crop) +
    tm_text("COUNTRY"
      , col   = "white"
      , just  = "left"
    ) +
  tm_layout(
      , legend.text.color       = "white"
      , legend.title.color      = "white"
      , legend.position = c("left", "bottom")
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_scale_bar(
      position    = "right"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
)

# Plot simulated paths
p2 <- tm_shape(bg) +
    tm_raster(
        palette     = "black"
      , legend.show = FALSE
    ) +
  tm_shape(paths) +
    tm_raster(
        palette = "viridis"
      , style   = "cont"
      , title   = "Simulated\nTrajectories"
      , labels  = c("Low-Frequency", "", "High-Frequency")
    ) +
  tm_shape(africa) +
    tm_lines(
        col = "white"
      , lwd = 1
      , lty = 2
    ) +
  tm_shape(bbox) +
    tm_lines(
        col = "white"
      , lwd = 1
      , lty = 1
    ) +
  tm_shape(africa_crop) +
    tm_text("COUNTRY"
      , col   = "white"
      , just  = "left"
    ) +
  tm_layout(
      , legend.text.color       = "white"
      , legend.title.color      = "white"
      , legend.position = c("left", "bottom")
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_scale_bar(
      position    = "right"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
)

# Plot simulated paths
p3 <- tm_shape(bg) +
    tm_raster(
        palette     = "black"
      , legend.show = FALSE
    ) +
  tm_shape(paths) +
    tm_raster(
        palette = "viridis"
      , style   = "cont"
      , title   = "Simulated\nTrajectories"
      , labels  = c("Low-Frequency", "", "High-Frequency")
    ) +
  tm_shape(tracks) +
    tm_lines(
        col = "white"
      , lwd = 0.5
    ) +
  tm_shape(africa) +
    tm_lines(
        col = "white"
      , lwd = 1
      , lty = 2
    ) +
  tm_shape(bbox) +
    tm_lines(
        col = "white"
      , lwd = 1
      , lty = 1
    ) +
  tm_shape(africa_crop) +
    tm_text("COUNTRY"
      , col   = "white"
      , just  = "left"
    ) +
  tm_layout(
      , legend.text.color       = "white"
      , legend.title.color      = "white"
      , legend.position = c("left", "bottom")
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_scale_bar(
      position    = "right"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
)

# Plot permeability surface
p4 <- tm_shape(perm) +
    tm_raster(
        palette = "viridis"
      , style   = "cont"
      , title   = "Permeability"
      , labels  = c("Low", "", "High")
    ) +
  tm_shape(africa) +
    tm_lines(
        col = "white"
      , lwd = 1
      , lty = 2
    ) +
  tm_shape(bbox) +
    tm_lines(
        col = "white"
      , lwd = 1
      , lty = 1
    ) +
  tm_shape(africa_crop) +
    tm_text("COUNTRY"
      , col   = "white"
      , just  = "left"
    ) +
  tm_layout(
      , legend.text.color       = "white"
      , legend.title.color      = "white"
      , legend.position = c("left", "bottom")
  ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
    , text.size   = 1.5
    , size        = 2
  ) +
  tm_scale_bar(
      position    = "right"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
)

# Store the plots
setwd(output)
png("99_Simulations1.png", width = 1280, height = 720, bg = "transparent")
p1 + tm_layout(scale = 3)
dev.off()
setwd(input)

setwd(output)
png("99_Simulations2.png", width = 1280, height = 720, bg = "transparent")
p2 + tm_layout(scale = 3)
dev.off()
setwd(input)

# Store the plots
setwd(output)
png("99_Simulations3.png", width = 1280, height = 720, bg = "transparent")
p3 + tm_layout(scale = 3)
dev.off()
setwd(input)

# Store the plots
setwd(output)
png("99_Simulations4.png", width = 1280, height = 720, bg = "transparent")
p4 + tm_layout(scale = 3)
dev.off()
setwd(input)

############################################################
#### Simulation Animation
############################################################
# Load tracks
sim <- shapefile("99_DispersalSimulation")

# Copy the object
tracks <- sim

# Convert to dataframe
tracks <- tracks@data

# Add counter to dataframe
tracks$No1 <- seq(from = 1, length = nrow(tracks), by = 2)
tracks$No2 <- seq(from = 2, length = nrow(tracks), by = 2)

# Put start and endpoints below each other
tracks1 <- tracks[, c("id", "x1_", "y1_", "No1")]
tracks2 <- tracks[, c("id", "x2_", "y2_", "No2")]

# Rename colums
names(tracks1) <- c("id", "x", "y", "No")
names(tracks2) <- c("id", "x", "y", "No")

# Bind the two datasets
tracks <- arrange(rbind(tracks1, tracks2), No)

# Remove duplicates
tracks$remove <- FALSE
for (i in 2:nrow(tracks)){
  if (
    tracks$x[i] == tracks$x[i-1] &
    tracks$y[i] == tracks$y[i-1] &
    tracks$id[i] == tracks$id[i-1]
  ){
    tracks$remove[i] <- TRUE
  }
}
tracks <- subset(tracks, !remove)

# We need to add timestamps. Lets create a vector
tracks <- split.data.frame(tracks, tracks$id)
Timestamp <- rep(ymd_hms("2000-01-01 00:00:00"), nrow(tracks[[1]]))
for (i in 1:length(Timestamp)){
  Timestamp[i] <- Timestamp[i] + hours(4) * i
}

# Add this vector to all datasets
tracks <- lapply(tracks, function(x){
  x$Timestamp <- Timestamp
  return(x)
})

# Bind them
tracks <- do.call(rbind, tracks)
head(tracks)

# Create a move object
move <- move(
    x       = tracks$x
  , y       = tracks$y
  , time    = tracks$Timestamp
  , data    = tracks
  , animal  = tracks$id
  , proj    = CRS("+init=epsg:32734")
)
move$colour <- "white"
move <- spTransform(move, CRS("+init=epsg:4326"))

# Align tracks
m <- align_move(move, res = 2, digit = 0, unit = "hours")

# Prepare frames
frames <- frames_spatial(
    m
  , tail_colour   = "white"
  , trace_colour  = "white"
  , trace_show    = TRUE
  , path_legend   = FALSE
  , path_size     = 0.50
  , tail_size     = 0.25
  , ext           = extent(sim)
  , map_service   = "mapbox"
  , map_type      = "satellite"
  , map_token     = "pk.eyJ1IjoiZG9keDkiLCJhIjoiY2p3dnltejJjMGR4YjN5bXp0ZjA2ZXBzMCJ9.4hirgQ-1SfJ2KHI7SR54cQ") %>%
  add_progress()
  # add_northarrow(colour = "white", size = 3) %>%
  # add_scalebar(colour = "white", distance = 50) %>%
  # add_labels(x = "Longitude", y = "Latitude")

# Look at the frames
frames[[200]]

# Store the animation
setwd(output)
animate_frames(
    frames
  , out_file = "Simulation2.mov"
  , overwrite = TRUE
  , width = 1980
  , height = 1080
  , end_pause = 2
)
setwd(input)

############################################################
#### Movement Model
############################################################
# Load required data
move <- readRDS("99_MovementModel.rds")

# Select best model
best <- move$Model[[1]]

# Get coefficients from best model
coeffs <- getCoeffs(best)[-1, ]
coeffs$Covariate[coeffs$Covariate == "Shrubs"] <- "Shrubs/Grassland"

# Prepare plot
p <- ggplot(data = coeffs, aes(y = Covariate, x = Coefficient)) +

  # Add points for the estimated coefficients
  geom_point(shape = 1, size = 1, col = "white") +

  # Add errorbars
  geom_errorbarh(aes(
      xmin = Coefficient - 1.96 * SE
    , xmax = Coefficient + 1.96 * SE
    , height = 0.2)
    , col = "white"
  ) +

  # Add a vertical line that indicates no effect
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +

  # Reverse the scale of the y-axis
  scale_y_discrete(limits = rev(c(
      "cos(ta_)"
    , "log(sl_)"
    , "Water"
    , "DistanceToWater"
    , "Trees"
    , "Shrubs/Grassland"
    , "HumansBuff5000"
    , "log(sl_):Water"
    , "cos(ta_):DistanceToWater"
    , "cos(ta_):HumansBuff5000"
  ))) +

  # Use a simple black and white theme
  theme_minimal() +

  # Scale axis to -1, 1
  xlim(c(-1, 1)) +

  # Change label of x-axis
  labs(x = expression(beta*"-Coefficient")) +

  # Cap the axes
  coord_capped_cart(left = "both", bottom = "both") +

  theme(
    plot.background = element_rect(fill = "transparent")
  , axis.text       = element_text(color = "white")
  , axis.title      = element_text(color = "white")
  , axis.ticks      = element_line(color = "white")
  , axis.line       = element_line(color = "white")
  , panel.grid      = element_line(color = "transparent")
)

# Write the plot to file
setwd(output)
ggsave("99_MovementModel.png"
  , plot    = p
  , width   = 14
  , height  = 12
  , scale   = 0.25
  , bg      = "transparent"
)
dev.off()
setwd(input)

# Prepare raster for contourplot of interactions
r1 <- visInt2(best, "log(sl_)", "Water")
r2 <- visInt2(best, "cos(ta_)", "DistanceToWater")
r3 <- visInt2(best, "cos(ta_)", "HumansBuff5000")

# Prepare the plots
p1 <- levelplot(r1
  , contour = TRUE
  , xlab    = "log(sl_)"
  , ylab    = "Water"
  , margin  = FALSE
  , main    = "(b1)"
  , par.settings = viridisTheme
)

p2 <- levelplot(r2
  , contour = TRUE
  , xlab    = "cos(ta_)"
  , ylab    = "DistanceToWater"
  , margin  = FALSE
  , main    = "(b2)"
  , par.settings = viridisTheme
)

p3 <- levelplot(r3
  , contour = TRUE
  , xlab    = "cos(ta_)"
  , ylab    = "HumansBuff5000"
  , margin  = FALSE
  , main    = "(b3)"
  , par.settings = viridisTheme
)

# Put the plots together and store them
setwd(output)
CairoPDF("99_MovementModel(Interactions).pdf"
  , width     = 14/1.5
  , height    = 5/1.5
  , pointsize = 40
)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()
setwd(input)

############################################################
#### Social Model
############################################################
# Load required data
social <- readRDS("99_SocialModel.rds")

# Select best model
best <- social

# Get coefficients from best model
coeffs <- getCoeffs(best)[-1, ]
coeffs$Covariate[coeffs$Covariate == "Shrubs"] <- "Shrubs/Grassland"

# Prepare plot
p <- ggplot(data = coeffs, aes(y = Covariate, x = Coefficient)) +

  # Add points for the estimated coefficients
  geom_point(shape = 1, size = 1, col = "white") +

  # Add errorbars
  geom_errorbarh(aes(
      xmin = Coefficient - 1.96 * SE
    , xmax = Coefficient + 1.96 * SE
    , height = 0.2)
    , col = "white"
  ) +

  # Add a vertical line that indicates no effect
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +

  # Reverse the scale of the y-axis
  scale_y_discrete(limits = rev(c(
      "cos(ta_)"
    , "log(sl_)"
    , "Water"
    , "DistanceToWater"
    , "Trees"
    , "Shrubs/Grassland"
    , "HumansBuff5000"
    , "ProbEncounter"
    , "OwnHomerange"
    , "ForeignHomeranges"
    , "ProbEncounter:OwnHomerange"
    , "ProbEncounter:ForeignHomeranges"
  ))) +

  # Use a simple black and white theme
  theme_minimal() +

  # Scale axis to -1, 1
  xlim(c(-1.7, 1.7)) +

  # Change label of x-axis
  labs(x = expression(beta*"-Coefficient")) +

  # Cap the axes
  coord_capped_cart(left = "both", bottom = "both") +

  theme(
    plot.background = element_rect(fill = "transparent")
  , axis.text       = element_text(color = "white")
  , axis.title      = element_text(color = "white")
  , axis.ticks      = element_line(color = "white")
  , axis.line       = element_line(color = "white")
  , panel.grid      = element_line(color = "transparent")
)

# Write the plot to file
setwd(output)
ggsave("99_SocialModel.png"
  , plot    = p
  , width   = 16
  , height  = 12
  , scale   = 0.25
  , bg      = "transparent"
)
dev.off()
setwd(input)

############################################################
#### Social Model: Example Landscape
############################################################
# load required data
setwd("/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/01_SmallExtent")
uds <- stack("05_SocialFeatures_SocialLandscape(UtilisationDistributions).grd")
hrs <- shapefile("05_SocialFeatures_SocialLandscape(HomeRanges)")

# Select a date to plot
names(uds)
hrs_sub <- subset(hrs, Date == "2017-07-01")
uds_sub <- crop(uds[["X2017.07.01"]], hrs_sub)

# Prepare plot of the uds
p1 <- tm_shape(uds_sub) +
  tm_raster(
      palette = "-Spectral"
    , style = "cont"
    , title = "Probability of Encounter"
    , labels = c("Low", "", "High")
  ) +
  tm_layout(
      , legend.text.color       = "white"
      , legend.title.color      = "white"
      , legend.position = c("left", "bottom")
  ) +
  tm_compass(
      text.color    = "white"
    , color.dark    = "white"
    , color.light   = "white"
  ) +
  tm_scale_bar(
      position    = "right"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
)

# Prepare a plot of the HRs
p2 <- tm_shape(uds_sub) +
    tm_raster(
      palette     = "white"
    , legend.show = FALSE
  ) +
  tm_shape(hrs_sub) +
    tm_polygons(
      "Pack"
    , palette = "viridis"
    , alpha = 0.5
  ) +
  tm_layout(
      , legend.text.color       = "black"
      , legend.title.color      = "black"
      , legend.position = c("left", "bottom")
  ) +
  tm_compass(
      text.color    = "black"
    , color.dark    = "black"
    , color.light   = "black"
  ) +
  tm_scale_bar(
      position    = "right"
    , text.size   = 0.5
    , text.color  = "black"
    , width       = 0.125
)

# Store the plot
setwd(output)
png("99_SocialLandscape1.png"
  , width   = 1080
  , height  = 1080
  , bg      = "transparent"
  , pointsize = 20
)
p1 + tm_layout(scale = 2)
dev.off()
setwd(input)

setwd(output)
png("99_SocialLandscape2.png"
  , width   = 1080
  , height  = 1080
  , bg      = "transparent"
  , pointsize = 20
)
p2 + tm_layout(scale = 2)
dev.off()
setwd(input)
