############################################################
#### Preparation of all Plots for the Thesis
############################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Set a seed
set.seed(12345)

# Load required packages
library(rgdal)
library(adehabitatLT)
library(NLMR)
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
library(ggmap)
library(mapview)
library(knitr)
library(kableExtra)
library(davidoff)

############################################################
#### General Statistics
############################################################
# Load required data and do some cleaning
info <- "03_Data/03_Results/99_GPSInfo.csv" %>%
  read_csv() %>%
  select(
      CoalitionID           = DogName
    , Sex                   = Sex
    , PackAffiliation       = CurrentPack
    , NumberFixes           = NoFixesTotal
    , NumberFixesDispersal  = NoFixesDispersal
    , DaysDispersing        = DaysDispersing
    , EculideanDistance     = StraightLineDistanceTravelled
    , CumulativeDistance    = DistanceTravelled
  ) %>%
  mutate(
      NumberFixes = formatC(NumberFixes, big.mark ="'")
    , CumulativeDistance = formatC(round(CumulativeDistance), big.mark ="'")
  )

# Convert to xtable
table <- xtable(info, auto = T, digits = 0)

# Write the table to a .tex table
print(table
  , floating            = FALSE
  , latex.environments  = NULL
  , booktabs            = TRUE
  , include.rownames    = FALSE
  , type                = "latex"
  , file                = "04_Manuscript/99_GPSInfo.tex"
)

############################################################
#### Plot of Africa
############################################################
# Load required data
africa <- "03_Data/02_CleanData/00_General_Africa.shp" %>% readOGR()
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>% readOGR()
dogs <- "03_Data/02_CleanData/00_General_WildDogs_IUCN.shp" %>% readOGR()

# Remove small islands for plotting
africa <- subset(africa, !(ID %in% c(27:41, 689, 690)))

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Prepare a map of Africa
p1 <- tm_shape(dogs) +
    tm_polygons(
        col         = "gray40"
      , border.col  = "black"
      , lwd         = 0.5
      , legend.show = F
    ) +
  tm_shape(africa, is.master = T) +
    tm_borders(
        col = "black"
      , lwd = 0.7
    ) +
  tm_shape(kaza_ext) +
    tm_borders(
        col = "red"
      , lty = 3
      , lwd = 1.5
    ) +
  tm_layout(
      asp         = 0.8
    , frame       = "black"
    , frame.lwd   = 3
    , legend.show = FALSE
  ) +
  tm_credits("(a)"
    , position  = c("left", "bottom")
    , size      = 1.5
    , col       = "black"
)

# Look at the plot
p1

############################################################
#### Plot of the KAZA
############################################################
# Load required data
africa <- "03_Data/02_CleanData/00_General_Africa.shp" %>% readOGR(.)
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>% readOGR(.)
prot <- "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(3Classes).shp" %>% readOGR(.)
water <- "03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters_GEOFABRIK.shp" %>% readOGR(.)

# Reclassify forest reserves as protected
prot$Desig[prot$Desig == "Forest Reserve"] <- "Protected"

# Specify the extent of the core study area
core <- extent(23.25, 24, -19.8, -19.05) %>% as(., "SpatialPolygons")
crs(core) <- CRS("+init=epsg:4326")

# Create labels for countries
labels_countries <- data.frame(
    x = c(20.39, 23.94, 20.07, 25.69, 28.22)
  , y = c(-15.28, -19.94, -19.39, -15.22, -18.9)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.8, 25.6)
  , y     = c(-19, -18.2, -17, -20.7)
  , Label = c(
    "Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans"
  )
)
coordinates(labels_waters) <- c("x", "y")
crs(labels_waters) <- CRS("+init=epsg:4326")

# Create labels for some national parks
labels_nationalparks <- data.frame(
    x = c(26.56, 28.61, 21.15, 25.87, 20.38, 23.58, 23.71, 24.51, 20.78)
  , y = c(-19.08, -17.05, -17.26, -14.66, -16.08, -21.4, -19.29, -18.65, -18.81)
  , Label = paste0(c(
      "Hwange", "Matusadona", "Luengue-Luiana", "Kafue", "Mavinga"
    , "Central Kalahari", "Moremi", "Chobe", "Khaudum"
  ), "\nNP")
)
coordinates(labels_nationalparks) <- c("x", "y")
crs(labels_nationalparks) <- CRS("+init=epsg:4326")

# Prepare a map of Kaza
p2 <- tm_shape(subset(prot, Desig == "Protected")) +
    tm_polygons(
        col           = "Desig"
      , palette       = viridis(20)[15]
      , border.col    = darken(viridis(20)[15], 1.2)
      , border.alpha  = 0.8
      , lwd           = 0.5
      , legend.show   = F
    ) +
  tm_shape(subset(prot, Desig == "National Park")) +
    tm_polygons(
        col           = "Desig"
      , palette       = viridis(20)[12]
      , border.col    = darken(viridis(20)[12], 1.2)
      , border.alpha  = 0.8
      , lwd           = 0.5
      , legend.show   = F
    ) +
  tm_shape(water) +
    tm_polygons(
        col           = "cornflowerblue"
      , border.col    = "darkblue"
      , border.alpha  = 0.6
      , lwd           = 0.2
    ) +
  tm_shape(kaza, is.master = T) +
    tm_borders(
        col = "blue"
      , lty = 1
      , lwd = 2
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "black"
    ) +
  tm_shape(core) +
    tm_borders(
        col = "white"
      , lty = 3
      , lwd = 1.5
    ) +
  tm_shape(labels_countries) +
    tm_text("Label"
      , fontface  = 2
      , size      = 1.5
    ) +
  tm_shape(labels_waters) +
    tm_text("Label"
      , fontface  = 3
      , size      = 0.5
      , col       = "darkblue"
    ) +
  tm_shape(labels_nationalparks) +
    tm_text("Label"
      , col       = "darkgreen"
      , fontface  = 3
      , size      = 0.5
    ) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
        legend.outside          = TRUE
      , legend.outside.position = "left"
      , legend.text.color       = "black"
      , legend.title.color      = "black"
      , legend.stack            = "vertical"
      , legend.text.size        = 0.8
      , legend.bg.color         = "white"
      , legend.frame            = "black"
      , frame                   = "black"
      , frame.lwd               = 3
      , asp                     = 1.2
  ) +
  tm_scale_bar(
        position  = "left"
      , text.size = 0.5
      , width     = 0.125
  ) +
  tm_compass(color.light = "black") +
  tm_credits("(b)"
    , position  = c("left", "top")
    , size      = 1.5
    , col       = "black"
)

# Look at the plot
p2

############################################################
#### Prepare a common Legend
############################################################
# We now prepare a legend that the two plots share. However, we have to assign
# to one of the two plots. Let's use p2 for this
p2 <- p2 + tm_add_legend(
    type = "fill"
  , labels = c(
      "Wild Dog Populations"
    , "Major Water Areas"
    , "National Parks"
    , "Other Protected Areas"
  )
  , col = c(
      "gray40"
    , "cornflowerblue"
    , viridis(20)[12]
    , viridis(20)[15]
  )) + tm_add_legend(
    type    = "line"
  , labels  = c("KAZA-TFCA Borders", "Country Borders")
  , col     = c("blue", "black")
  , lty     = c(1, 1)
  , lwd     = c(2, 2)
) + tm_layout(legend.frame.lwd = 2, legend.text.size = 1.05)

# Combine plots and store them
CairoPDF("04_Manuscript/99_StudyArea.pdf", width = 9, height = 5.25)
p2
print(p1, vp = viewport(0.158, 0.341, width = 0.56, height = 0.56))
dev.off()

############################################################
#### Plot of the Dispersal Trajectories
############################################################
# Define the extent that we want to plot the trajectories on
extent <- extent(22, 27, -20.5, -17.5) %>% as(., "SpatialPolygons")
crs(extent) <- CRS("+init=epsg:4326")

# Specify the extent of the core study area
core <- extent(23.25, 24, -19.8, -19.05) %>% as(., "SpatialPolygons")
crs(core) <- CRS("+init=epsg:4326")

# Load required data
africa <- "03_Data/02_CleanData/00_General_Africa.shp" %>%
  readOGR(.) %>%
  crop(., extent(extent) + c(-1, 1, -1, 1))
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>%
  readOGR(.) %>%
  crop(., extent(extent) + c(-1, 1, -1, 1))
movement <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv(.) %>%
  subset(., State == "Disperser")

# Prepare spatial dots from the gps fixes
movement_p <- movement
coordinates(movement_p) <- c("x", "y")
crs(movement_p) <- CRS("+init=epsg:4326")

# Prepare spatial lines from the gps fixes
movement_l <- movement %>%
  group_by(DogName) %>%
  nest(.) %>%
  mutate(Lines = map(data, function(x){
    points <- SpatialPoints(x[, c("x", "y")])
    lines <- spLines(points)
  }))
movement_l <- do.call(rbind, movement_l$Lines)
movement_l$DogName <- unique(movement$DogName)
crs(movement_l) <- CRS("+init=epsg:4326")

# Get a satellite image for the extent of which we will plot the trajectories
# For whatever reason they switched from "raster" to "stars" objects
map <- extent %>%
  read_osm(., type = "bing", zoom = 8) %>%
  as(., "Raster") %>%
  projectRaster(., crs = CRS("+init=epsg:4326"), method = "ngb") %>%
  trim(.)

# Check the map
plotRGB(map)

# Create country labels on specified position
labels_countries <- data.frame(
    x     = c(22.5, 24, 25, 26.5, 26)
  , y     = c(-17.7, -17.8, -18.3, -19, -17.6)
  , Label = c("Angola", "Namibia", "Botswana", "Zimbabwe", "Zambia")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.8, 25.4)
  , y     = c(-19, -18.3, -17, -20.7)
  , Label = c(
    "Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Nwetwe\nPan"
  )
)
coordinates(labels_waters) <- c("x", "y")
crs(labels_waters) <- CRS("+init=epsg:4326")

# I will use a specific order of the viridis colors
order <- c(20, 17, 2, 3, 19, 4, 1, 5, 17, 16, 6, 15, 14, 7, 8, 9, 10, 11, 12, 13)

# Prepare a map of the trajectories
p3 <- tm_shape(map) +
    tm_rgb() +
  tm_shape(africa) +
    tm_borders(
        col = "white"
      , lwd = 1
    ) +
  tm_shape(movement_p) +
    tm_dots(
        col     = "DogName"
      , palette = viridis(20)[order]
      , size    = 0.05
      , shape   = 1
    ) +
  tm_shape(movement_l) +
    tm_lines(
        col     = "DogName"
      , palette = viridis(20)[order]
      , lwd     = 1.5
    ) +
  tm_shape(core) +
    tm_borders(
        col = "white"
      , lty = 2
    ) +
  tm_shape(labels_countries) +
    tm_text(
        "Label"
      , col       = "white"
      , fontface  = 2
    ) +
  tm_shape(labels_waters) +
    tm_text(
        "Label"
      , col       = "white"
      , fontface  = 3
      , size      = 0.6
    ) +
  tm_compass(
      color.dark  = "white"
    , text.color  = "white"
    , position    = c("right", "bottom")
  ) +
  tm_scale_bar(
      width       = 0.125
    , text.color  = "white"
    , position    = c("left", "bottom")
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
    legend.show = FALSE
)

# Store the plot
CairoPDF("04_Manuscript/99_Trajectories.pdf", width = 10, height = 6.5)
p3
dev.off()

############################################################
#### Plot of NSD
############################################################
# Load data
dat <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv() %>%
  subset(DogName == "Odzala" & Timestamp > "2018-07-01 00:00:00")

# Convert to utm
coordinates(dat) <- c("x", "y")
crs(dat) <- CRS("+init=epsg:4326")
dat <- spTransform(dat, CRS("+init=epsg:32734"))
dat <- as.data.frame(dat, xy = T)

# Specify reference points
x1 <- dat$x[1]
y1 <- dat$y[1]

# Calculate NSD
dat$NSD <- (dat$x - x1) ** 2 + (dat$y - y1) ** 2

# Prepare plot
p <- ggplot(dat, aes(x = Timestamp, y = NSD, color = State)) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.3) +
  scale_color_manual(values = c("blue", "black")) +
  theme_linedraw() +
  aes(group = NA) +
  labs(
      color = "Behavioral Mode"
    , x = "Date"
    , y = "Net Squared Displacement (NSD)"
  ) +
  theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    , legend.position = c(.95, .95)
    , legend.justification = c("right", "top")
    , legend.box.just = "right"
    , legend.margin = margin(6, 6, 6, 6)
  ) +
  theme(panel.grid = element_blank()) +
  scale_x_datetime(labels = date_format("%Y-%m"))

# Store the plot
CairoPDF("04_Manuscript/99_NSD.pdf", width = 8, height = 4.5)
p
dev.off()

############################################################
#### Plot of NSD of All Dispersers
############################################################
# Load data
dat <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv()

# We only want individuals that dispersed at some point
disps <- subset(dat, State == "Disperser") %>%
  select(DogName) %>%
  distinct()
dat <- subset(dat, DogName %in% disps$DogName)

# Convert to utm
coordinates(dat) <- c("x", "y")
crs(dat) <- CRS("+init=epsg:4326")
dat <- spTransform(dat, CRS("+init=epsg:32734"))
dat <- as.data.frame(dat, xy = T)

# Nest
dat <- dat %>%
  group_by(DogName) %>%
  nest()

# Specify reference points
dat <- mutate(dat, data = map(data, function(z){
  x1 <- z$x[1]
  y1 <- z$y[1]
  z$NSD <- (z$x - x1) ** 2 + (z$y - y1) ** 2
  return(z)
}))

# Unnest
dat <- unnest(dat)

# Prepare plot
p <- ggplot(dat, aes(x = Timestamp, y = NSD, color = State)) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.3) +
  scale_color_manual(values = c("blue", "black")) +
  theme_linedraw() +
  aes(group = NA) +
  labs(
      color = "Behavioral Mode"
    , x = "Date"
    , y = "Net Squared Displacement (NSD)"
  ) +
  theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    , legend.position = c(.05, .95)
    , legend.justification = c("right", "top")
    , legend.box.just = "right"
    , legend.margin = margin(6, 6, 6, 6)
  ) +
  theme(panel.grid = element_blank()) +
  scale_x_datetime(labels = date_format("%Y-%m")) +
  facet_wrap(. ~ DogName, scales = "free")

# Store the plot
CairoPDF("04_Manuscript/99_NSD(All).pdf", width = 32, height = 18)
p
dev.off()

############################################################
#### Covariates
############################################################
# Prepare a bounding box for the desired extent
bbox <- extent(c(21.74909, 24.30089, -20.64901, -18.14901))

# Prepare the extent of our core area
core <- as(extent(c(23.25, 24, -19.8, -19.05)), "SpatialPolygons")
crs(core) <- CRS("+init=epsg:4326")
core <- as(core, "SpatialLines")

# Load required data
water   <- "03_Data/02_CleanData/01_LandCover_Water(Merged).tif" %>% stack()
trees   <- "03_Data/02_CleanData/01_LandCover_TreeCover_MODIS.tif" %>% raster()
shrubs  <- "03_Data/02_CleanData/01_LandCover_NonTreeVegetation_MODIS.tif" %>% raster()
prot    <- "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(1Class).tif" %>% raster()
humans  <- "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Buffer5000).tif" %>% raster()
roads1  <- "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToRoads.tif" %>% raster()
roads2  <- "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.shp" %>% readOGR()

# Crop all data to an extent that is slightly larger than the bbox
water   <- crop(water, extent(bbox) * 1.2) %>% calc(., fun = mean)
trees   <- crop(trees, extent(bbox) * 1.2) %>% calc(., fun = function(x){100*x})
shrubs  <- crop(shrubs, extent(bbox) * 1.2) %>% calc(., fun = function(x){100*x})
prot    <- crop(prot, extent(bbox) * 1.2)
humans  <- crop(humans, extent(bbox) * 1.2)
roads1  <- crop(roads1, extent(bbox) * 1.2) %>% calc(., fun = function(x){x/1000})
roads2  <- crop(roads2, extent(bbox) * 1.2)

# Remove outliers in the tree layer (for better visuals)
upper <- quantile(values(trees), 0.99, na.rm = TRUE)
lower <- quantile(values(trees), 0.01, na.rm = TRUE)
values(trees)[values(trees) > upper] <- upper
values(trees)[values(trees) < lower] <- lower

# Remove outliers in the shrub layer (for better visuals)
upper <- quantile(values(shrubs), 0.99, na.rm = TRUE)
lower <- quantile(values(shrubs), 0.01, na.rm = TRUE)
values(shrubs)[values(shrubs) > upper] <- upper
values(shrubs)[values(shrubs) < lower] <- lower

# Rescale all continuous maps to a range between 0 and 1
water <- calc(water, fun = function(x){(x - min(x)) / (max(x) - min(x))})
trees <- calc(trees, fun = function(x){(x - min(x)) / (max(x) - min(x))})
shrubs <- calc(shrubs, fun = function(x){(x - min(x)) / (max(x) - min(x))})
roads1 <- calc(roads1, fun = function(x){(x - min(x)) / (max(x) - min(x))})
humans <- calc(humans, fun = function(x){(x - min(x)) / (max(x) - min(x))})

# Add some pseudo data to the core study area lines
core$Pseudo <- ""

# Plot the layers
p1 <- tm_shape(water, bbox = bbox) +
  tm_raster(
      palette         = "viridis"
    , style           = "cont"
    , title           = "Inundation\nFrequency"
    , labels          = c("Low", "", "High")
    , legend.reverse  = T
  ) +
  tm_shape(core) +
    tm_lines("Pseudo"
      , palette         = c("white")
      , lwd             = 1.5
      , lty             = 2
      , title.col       = "Core Study Area"
      , legend.col.show = FALSE
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
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_credits("(a)"
    , position  = c("center", "top")
    , size      = 1.5
    , col       = "white"
)

p2 <- tm_shape(trees, bbox = bbox) +
  tm_raster(
      palette         = "viridis"
    , style           = "cont"
    , title           = "Tree Cover (%)"
    , labels          = c("Low", "", "High")
    , legend.reverse  = T
  ) +
  tm_shape(core) +
    tm_lines("Pseudo"
      , palette         = c("white")
      , lwd             = 1.5
      , lty             = 2
      , title.col       = "Core Study Area"
      , legend.col.show = FALSE
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
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_credits("(b)"
    , position  = c("center", "top")
    , size      = 1.5
    , col       = "white"
)

p3 <- tm_shape(shrubs, bbox = bbox) +
  tm_raster(
      palette         = "viridis"
    , style           = "cont"
    , title           = "Shrub/Grassland\nCover (%)"
    , labels          = c("Low", "", "High")
    , legend.reverse  = T
  ) +
  tm_shape(core) +
    tm_lines("Pseudo"
      , palette         = c("white")
      , lwd             = 1.5
      , lty             = 2
      , title.col       = "Core Study Area"
      , legend.col.show = FALSE
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
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_credits("(c)"
    , position  = c("center", "top")
    , size      = 1.5
    , col       = "white"
)

p4 <- tm_shape(prot, bbox = bbox) +
  tm_raster(
    , palette         = viridis(10)[c(1, 8)]
    , style           = "cat"
    , title           = "Protection\nStatus"
    , labels          = c("Pastoral", "Proctected")
    , legend.reverse  = T
  ) +
  tm_shape(core) +
    tm_lines("Pseudo"
      , palette         = c("white")
      , lwd             = 1.5
      , lty             = 2
      , title.col       = "Core Study Area"
      , legend.col.show = FALSE
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
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_credits("(d)"
    , position  = c("center", "top")
    , size      = 1.5
    , col       = "white"
)

p5 <- tm_shape(humans, bbox = bbox) +
  tm_raster(
      palette         = "viridis"
    , style           = "cont"
    , title           = "Human\nInfluence"
    , labels          = c("Low", "", "High")
    , legend.reverse  = T
  ) +
  tm_shape(core) +
    tm_lines("Pseudo"
      , palette         = c("white")
      , lwd             = 1.5
      , lty             = 2
      , title.col       = "Core Study Area"
      , legend.col.show = FALSE
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
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_credits("(e)"
    , position  = c("center", "top")
    , size      = 1.5
    , col       = "white"
)

p6 <- tm_shape(roads1, bbox = bbox) +
  tm_raster(
      palette         = "viridis"
    , style           = "cont"
    , title           = "Distance to\nRoads (km)"
    , labels          = c("Short", "", "Long")
    , legend.reverse  = T
  ) +
  tm_shape(roads2, bbox = bbox) +
    tm_lines(col = "white") +
  tm_shape(core) +
    tm_lines("Pseudo"
      , palette         = c("white")
      , lwd             = 1.5
      , lty             = 2
      , title.col       = "Core Study Area"
      , legend.col.show = FALSE
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
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_credits("(f)"
    , position  = c("center", "top")
    , size      = 1.5
    , col       = "white"
)

# Put all plots together
p <- tmap_arrange(p1, p2, p3, p4, p5, p6, ncol = 3, asp = 1.1)

# Store them
CairoPDF("04_Manuscript/99_Covariates.pdf", width = 10.8, height = 6.3)
p
dev.off()

############################################################
#### The ORI Algorithm
############################################################
# Load the required data
wet_mask  <- "03_Data/01_RawData/MODIS/MCD43A4/02_Masks/MaskWater.shp" %>% readOGR()
dry_mask  <- "03_Data/01_RawData/MODIS/MCD43A4/02_Masks/MaskDryland.shp" %>% readOGR()
kaza      <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>% readOGR()
water     <- "03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters_GEOFABRIK.shp" %>% readOGR()
africa    <- "03_Data/02_CleanData/00_General_Africa.shp" %>% readOGR()

# Also load one of the classified MODIS layers
flood <- dir(
      path        = "03_Data/02_CleanData/00_Floodmaps/01_Original"
    , pattern     = ".tif$"
    , full.names  = T
  )[500] %>%
  raster()

# Make it binary
values(flood)[values(flood) == 255] <- 1

# Crop the africa layer to the kaza polygon
africa <- crop(africa, kaza)

# Put the masks together into one dataframe
wet_mask$LandCover <- "Water"
dry_mask$LandCover <- "Dryland"
masks <- rbind(dry_mask, wet_mask)

# Define a bounding box for which we calculated dynamic floodmaps
bbox <- extent(c(21.74909, 24.30089, -20.64901, -18.14901))
bbox <- as(bbox, "SpatialPolygons")
crs(bbox) <- CRS("+init=epsg:4326")
bbox_ext <- bb(bbox)

# Download Satellite Image as background map
map <- read_osm(bbox_ext, type = "bing") %>%
  as(., "Raster") %>%
  projectRaster(., crs = CRS("+init=epsg:4326"), method = "ngb") %>%
  trim(.)

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

# Store the plot
CairoPDF("04_Manuscript/99_Floodmapping.pdf", width = 9, height = 4.45)
p3
dev.off()

############################################################
#### The ORI Algorithm: Validation
############################################################
# Load required data
own <- "03_Data/02_CleanData/00_Floodmaps/03_Validation/2018.07.20.tif" %>%
  raster()
ori <- "03_Data/02_CleanData/00_Floodmaps/01_Original/2018.07.20.tif" %>%
  raster()

# We only care about water and dryland, but not about clouds. Let's reclassify
# the images so that water becomes 1, dryland 0 and clouds 0 as well
rcl <- data.frame(old = c(0, 127, 255), new = c(1, 0, 0))
own <- reclassify(own, rcl)
ori <- reclassify(ori, rcl)

# Function to create difference map
falseMap <- function(x, y, z){
  x[[z]] - y[[z]]
}

# Apply it to our classification
false <- falseMap(own, ori, 1)
values(false)[values(false) == 0] <- NA

# Prepare the plots
p1 <- tm_shape(ori) +
    tm_raster(
        style   = "cat"
      , palette = c("white", "blue")
      , labels  = c("Dryland", "Water")
      , title   = ""
    ) +
  tm_compass(
      color.dark  = "black"
    , color.light = "black"
    , text.color  = "black"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "black"
  ) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
      title.position  = c("center", "top")
    , legend.width    = 0.5
    , legend.bg.color = "white"
    , legend.position = c("left", "bottom")
  ) +
  tm_credits("(a)"
  , position  = c("center", "top")
  , size      = 1.5
)

p2 <- tm_shape(own) +
    tm_raster(
        style   = "cat"
      , palette = c("white", "blue")
      , labels  = c("Dryland", "Water")
      , title   = ""
    ) +
  tm_compass(
      color.dark  = "black"
    , color.light = "black"
    , text.color  = "black"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "black"
  ) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
      title.position  = c("center", "top")
    , legend.width    = 0.5
    , legend.bg.color = "white"
    , legend.position = c("left", "bottom")
  ) +
  tm_credits("(b)"
  , position  = c("center", "top")
  , size      = 1.5
)

p3 <- tm_shape(false) +
    tm_raster(
        style   = "cat"
      , palette = c("red", "darkblue")
      , labels  = c("False Negatives", "False Positives")
      , title   = ""
    ) +
  tm_compass(
      color.dark  = "black"
    , color.light = "black"
    , text.color  = "black"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "black"
  ) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
      title.position  = c("center", "top")
    , legend.width    = 0.5
    , legend.bg.color = "white"
    , legend.position = c("left", "bottom")
  ) +
  tm_credits("(c)"
  , position  = c("center", "top")
  , size      = 1.5
)

# Put the plots together
p <- tmap_arrange(p1, p2, p3, ncol = 3)

# Store the final plot
CairoPDF("04_Manuscript/99_FloodmappingValidation.pdf", width = 9, height = 2.8)
p
dev.off()

############################################################
#### Road Classes
############################################################
# Load required data
legend <- "03_Data/01_RawData/GEOFABRIK/RoadsDescription.csv" %>%
  read_csv2(locale = locale(encoding = "latin1"))

# We don't need to keep the key as it doesn't change
legend <- select(legend, -c(Key))

# Rename columns
names(legend) <- c("Group", "Subgroup", "Description")

# Write the table to a .tex table
print(xtable(legend)
  , floating            = FALSE
  , latex.environments  = NULL
  , booktabs            = TRUE
  , include.rownames    = FALSE
  , type                = "latex"
  , file                = "04_Manuscript/99_RoadsDescription.tex"
)

############################################################
#### Proxy for Human Influence
############################################################
# Load required data
population  <- "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanDensity_Facebook.tif" %>% raster()
farms_gabs  <- "03_Data/02_CleanData/04_AnthropogenicFeatures_Farms_Gabriele.tif" %>% raster()
farms_crops <- "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_Croplands.tif" %>% raster()
farms_glob  <- "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_Globelands.tif" %>% raster()
roads       <- "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.tif" %>% raster()
influence   <- "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Buffer5000).tif" %>% raster()
delete      <- "03_Data/01_RawData/DAVID/AnthropogenicInfluence(Delete).shp" %>% readOGR()

# Merge farm layers
farms <- max(farms_gabs, farms_crops, farms_glob)

# Prepare a bounding box for which we will plot the data
box <- extent(c(21.74909, 24.30089, -20.64901, -18.14901))

# Crop rasters to bounding box
farms       <- crop(farms, box)
population  <- crop(population, box)
roads       <- crop(roads, box)
influence   <- crop(influence, box)
delete      <- crop(delete, box)

# Remove human influence below the "delete" shapefile
farms[delete] <- 0
population[delete] <- 0

# Prepare a dataframe with some locations that we want to plot
villages <- data.frame(
    x = c(23.319332, 24.178607, 24.035000, 23.653653)
  , y = c(-19.799149, -19.430306, -19.000306, -19.002809)
  , Name = c("Maun", "Sankuyo", "Mababe", "Khwai")
)
coordinates(villages) <- c("x", "y")
crs(villages) <- CRS("+init=epsg:4326")

# Also prepare some pointers for the villages
l1 <- cbind(c(23.339332, 23.419332), c(-19.839149, -19.999149)) %>% spLines()
l2 <- cbind(c(24.050000, 23.890000), c(-19.430306, -19.427306)) %>% spLines()
l3 <- cbind(c(24.032000, 24.002000), c(-19.042809, -19.170306)) %>% spLines()
l4 <- cbind(c(23.653653, 23.753653), c(-19.042809, -19.152809)) %>% spLines()
villages_pointers <- rbind(l1, l2, l3, l4)
crs(villages_pointers) <- crs(villages)

# Plot of each covariate
p1 <- tm_shape((population + 1), bbox = box) +
    tm_raster(
        style   = "log10"
      , palette = viridis(50)
      # , palette = "Greys"
      , title   = "Human Density"
      , labels  = c("Low", "", "High")
    ) +
  tm_shape(villages) +
    tm_text("Name"
    , col   = "white"
    , size  = 0.4
  ) +
  tm_shape(villages_pointers) +
    tm_lines(
      col = "white"
    , lwd = 0.4
  ) +
  tm_compass(
      color.dark  = "white"
    , color.light = "white"
    , text.color  = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
  ) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
      title.position  = c("center", "top")
    , legend.width    = 0.5
    , legend.bg.color = "white"
    , legend.position = c("left", "bottom")
  ) +
  tm_credits("(a)"
  , position  = c("center", "top")
  , col       = "white"
  , size      = 1.5
)

p2 <- tm_shape(farms, bbox = box) +
    tm_raster(
        style   = "cat"
      , palette = viridis(20)[c(1, 20)]
      # , palette = c("white", "black")
      , title   = "Farming"
      , labels  = c("Absent", "Present")
    ) +
  tm_compass(
      color.dark  = "white"
    , color.light = "white"
    , text.color  = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
  ) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
      title.position  = c("center", "top")
    , legend.width    = 0.5
    , legend.bg.color = "white"
    , legend.position = c("left", "bottom")
  ) +
  tm_credits("(b)"
  , position  = c("center", "top")
  , col       = "white"
  , size      = 1.5
)

p3 <- tm_shape(roads, bbox = box) +
    tm_raster(
        style   = "cat"
      , palette = viridis(20)[c(1, 20)]
      # , palette = c("white", "black")
      , title   = "Roads"
      , labels  = c("Absent", "Present")
    ) +
  tm_compass(
      color.dark  = "white"
    , color.light = "white"
    , text.color  = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
  ) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
      title.position  = c("center", "top")
    , legend.width    = 0.5
    , legend.bg.color = "white"
    , legend.position = c("left", "bottom")
  ) +
  tm_credits("(c)"
  , position  = c("center", "top")
  , col       = "white"
  , size      = 1.5
)

p4 <- tm_shape(influence, bbox = box) +
    tm_raster(
        style   = "cont"
      , palette = viridis(20)
      # , palette = "Greys"
      , title   = "Human Influence"
      , labels  = c("Low", "", "High")
    ) +
  tm_shape(villages) +
    tm_text("Name"
    , col   = "white"
    , size  = 0.4
  ) +
  tm_shape(villages_pointers) +
    tm_lines(
      col = "white"
    , lwd = 0.4
  ) +
  tm_compass(
      color.dark  = "white"
    , color.light = "white"
    , text.color  = "white"
  ) +
  tm_scale_bar(
      breaks      = c(0, 25, 50)
    , text.color  = "white"
  ) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
      title.position  = c("center", "top")
    , legend.width    = 0.5
    , legend.bg.color = "white"
    , legend.position = c("left", "bottom")
  ) +
  tm_credits("(d)"
  , position  = c("center", "top")
  , col       = "white"
  , size      = 1.5
)

# Put all plots together
p <- tmap_arrange(p1, p2, p3, p4)

# Save it to file
CairoPDF("04_Manuscript/99_HumanInfluence.pdf", width = 8, height = 7.85)
p
dev.off()

############################################################
#### Permeability Model AICs
############################################################
# Load required data
aics <- read_rds("03_Data/03_Results/99_PermeabilityModelAICs.rds")

# Remove any entries with either "HumansAverage" or "HumansBase"
aics <- aics[!grepl(pattern = "HumansAverage|HumansBase", aics$Covariates), ]

# Also remove the ModelID column
aics$ModelID <- NULL

# Replace covariates with short code and add cos(ta) and log(sl) as covariatess
aics$Covariates <- aics$Covariates %>%
  gsub(pattern = "\\bWater\\b", replacement = "W") %>%
  gsub(pattern = "\\bShrubs\\b", replacement = "S") %>%
  gsub(pattern = "\\bTrees\\b", replacement = "T") %>%
  gsub(pattern = "\\bProtected\\b", replacement = "P") %>%
  gsub(pattern = "\\bHumansBuff5000\\b", replacement = "HI") %>%
  gsub(pattern = "\\bDistanceToWater\\b", replacement = "DTW") %>%
  gsub(pattern = "\\bDistanceToRoads\\b", replacement = "DTR") %>%
  gsub(pattern = "\\bRoadCrossing\\b", replacement = "RC") %>%
  gsub(pattern = ",", replacement = " +") %>%
  paste("cos(ta) + sl + log(sl) +", .)

# We will also add cos(ta) and log(sl) as covariates
aics

# Write the table to a .tex table
print(xtable(aics)
  , floating            = FALSE
  , latex.environments  = NULL
  , booktabs            = TRUE
  , include.rownames    = FALSE
  , type                = "latex"
  , file                = "04_Manuscript/99_PermeabilityModelAICs.tex"
)

############################################################
#### Permeability Model Results
############################################################
# Load required data
model <- read_rds("03_Data/03_Results/99_ModelSelection.rds")

# Select best model
best <- model$Model[[1]]

# Summary of the best model
summary(best)

# Prepare a table that shows the results
best_table <- getCoeffs(best, pvalue = TRUE) %>%

  # Calculate confidence intervals
  mutate(
      LCI = Coefficient - 1.96 * SE
    , UCI = Coefficient + 1.96 * SE
  )

# Rename Covariates
best_table$Covariate[best_table$Covariate == "Shrubs"] <- "Shrubs/Grassland"
best_table$Covariate[best_table$Covariate == "cos_ta_"] <- "cos(ta)"
best_table$Covariate[best_table$Covariate == "sl_"] <- "sl"
best_table$Covariate[best_table$Covariate == "log_sl_"] <- "log(sl)"
best_table$Covariate[best_table$Covariate == "HumansBuff5000"] <- "HumanInfluence"


# Print the result as a table to a tex file
options(scipen = 999)
print(xtable(best_table)
  , floating            = FALSE
  , latex.environments  = NULL
  , booktabs            = TRUE
  , include.rownames    = FALSE
  , type                = "latex"
  , file                = "04_Manuscript/99_PermeabilityModel.tex"
)

# Prepare and store a plot of the coefficients
coeffs <- getCoeffs(best, pvalue = T)[-1, ]

# Add stars indicating the significance
coeffs$Significance <- sapply(1:nrow(coeffs), function(x){
  if (coeffs$pvalue[x] <= 0.01){
    return("***")
  } else if (coeffs$pvalue[x] <= 0.05){
    return("**")
  } else if (coeffs$pvalue[x] <= 0.1){
    return("*")
  }
})

# Rename covariates nicely
coeffs$Covariate[coeffs$Covariate == "Shrubs"] <- "Shrubs/Grassland"
coeffs$Covariate[coeffs$Covariate == "cos_ta_"] <- "cos(ta)"
coeffs$Covariate[coeffs$Covariate == "sl_"] <- "sl"
coeffs$Covariate[coeffs$Covariate == "log_sl_"] <- "log(sl)"
coeffs$Covariate[coeffs$Covariate == "HumansBuff5000"] <- "HumanInfluence"
p1 <- showCoeffs(coeffs
  , shape     = 1
  , size      = 1
  , whiskers  = 1.96
  , order     = c(
      "cos(ta)"
    , "sl"
    , "log(sl)"
    , "Water"
    , "DistanceToWater"
    , "Trees"
    , "Shrubs/Grassland"
    , "HumanInfluence"
  )) +
  geom_text(aes(label = Significance, hjust = 0.5, vjust = -0.2), show.legend = F)

############################################################
#### Permeability Model Validation
############################################################
# Load required data
validation <- read_rds("03_Data/03_Results/99_ModelValidation.rds")
dat_pref <- read_rds("03_Data/03_Results/99_ModelValidation(Data).rds")[[1]]
dat_rand <- read_rds("03_Data/03_Results/99_ModelValidation(Data).rds")[[2]]

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
        paste0(
            "$\\bar{r}_s = "
          , validation$Mean[1]
          , "$, $95%-CI = "
          , validation$LCL[1]
          , "$, $"
          , validation$UCL[1], "$"
        )
      , paste0(
            "$\\bar{r}_s = "
          , validation$Mean[2]
          , "$, $95%-CI = "
          , validation$LCL[2]
          , "$, $"
          , validation$UCL[2], "$"
        )
    #   "$\\bar{r}_s = -0.45$, $95%-CI = -0.48$, $-0.42$"
    # , "$\\bar{r}_s = 0.04$, $95%-CI = 0.01$, $0.08$"
  )
)

# Reorder the factors in the Group variable
dat$Group <- factor(dat$Group, levels = c("Realized", "Random"))

# Plot the data
p2 <- ggplot(data = dat, aes(x = Rank, y = Frequency)) +
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

# Add titles to the plots
p1 <- p1 + ggtitle("(a)") + theme(plot.title = element_text(hjust = -0.8))
p2 <- p2 + ggtitle("(b)")

# Write the plot to file
CairoPDF(
    "04_Manuscript/99_PermeabilityResults.pdf"
  , width   = 8
  , height  = 2.8
)
grid.arrange(p1, p2, ncol = 2, widths = c(0.4, 0.6))
dev.off()

############################################################
#### Permeability Model Random Effects
############################################################
# Reload model
mod <- read_rds("03_Data/03_Results/99_ModelSelection.rds")
mod <- mod$Model[[1]]

# Check out coefficients per individual
coeffs <- coef(mod)
ranefs <- ranef(mod, condVar = T)
coeffs$cond$id
ranefs$cond$id

# Calculate the SD of the REs
apply(coeffs$cond$id, 2, sd)
apply(ranefs$cond$id, 2, sd)

# Note: ranef yields the difference between the individual specific effect and
# the mean level effect. coef, on the other hand, yields the individual specific
# effect. Thus, the followin two lines yield (approximately) the same
mean(coeffs$cond$id$cos_ta_) + ranefs$cond$id$cos_ta_[3]

# We now want to visualize the individual variation. There are two possibilities
# for this: lme4::dotplot() or a ggplot. The dotplot is easier, yet not
# customizable. Let's first do the dotplot, then recreate it in ggplot.
lme4:::dotplot.ranef.mer(ranef(mod)$cond)

# Maybe scalefree?
lme4:::dotplot.ranef.mer(ranef(mod)$cond, scales = list(x = list(relation = "free")))

# Prepare dataframe that we need to plot the same in ggplot
rfs <- ranefs$cond$id %>%
  rownames_to_column() %>%
  gather(key = Covariate, value = Mean, 2:9)
names(rfs)[1] <- "id"

# We need to add the conditional variance
condVar <- attributes(ranefs$cond$id)$condVar
names(condVar) <- attributes(ranefs$cond$id)$names
condVar <- as.data.frame(do.call(rbind, condVar))
names(condVar) <- attributes(ranefs$cond$id)$row.names
condVar <- rownames_to_column(condVar)
names(condVar)[1] <- "Covariate"
condVar <- gather(condVar, key = id, value = Variance, 2:17)

# Join data to rfs dataframe
rfs <- left_join(rfs, condVar)

# Rename stuff nicely
rfs$Covariate <- gsub(rfs$Covariate, pattern = "cos_ta_", replacement = "cos(ta)")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "log_sl_", replacement = "log(sl)")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "sl_", replacement = "sl")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "HumansBuff5000", replacement = "HumanInfluence")

# Make covariates a factor
rfs$Covariate <- factor(rfs$Covariate, levels = c(
  "cos(ta)", "sl", "log(sl)", "Water", "DistanceToWater", "Shrubs"
  , "Trees", "HumanInfluence"
))

# # # Create running number for individuals
# rfs$id <- as.numeric(as.factor(rfs$id))
# rfs$id <- sprintf("%02d", rfs$id)

# Visualize. Note that I am transforming the variance using mean - 2 *
# sqrt(Variance). This was taken from here: https://stackoverflow.com/questions
# /13847936/plot-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot-
# how-to-mak
p <- ggplot(rfs, aes(x = Mean, y = id)) +
  geom_point() +
  facet_wrap("Covariate", nrow = 2) +
  geom_errorbarh(aes(
      xmin = Mean - 2 * sqrt(Variance)
    , xmax = Mean + 2 * sqrt(Variance)
  ), colour = "black", height = 0) +
  xlim(-2.1, 2.1) +
  labs(y = "Coalition ID") +
  labs(x = expression(beta*"-Coefficient"))

# Store the plot
CairoPDF(
    "04_Manuscript/99_RandomEffects.pdf"
  , width   = 8
  , height  = 6
)
p
dev.off()

# Check out the variation in effects
as.data.frame(apply(coeffs$cond$id, 2, mean))
as.data.frame(apply(coeffs$cond$id, 2, sd))
as.data.frame(apply(coeffs$cond$id, 2, var))

# Put into a dataframe
dat <- data.frame(
    Covariate = names(coeffs$cond$id)
  , Mean_RE   = apply(coeffs$cond$id, 2, mean)
  , SD_RE     = apply(coeffs$cond$id, 2, sd)
)
dat <- dat[-1, ]
rownames(dat) <- NULL

# Convert to xtable
dat <- xtable(dat, auto = T, digits = 3)

# Write the table to a .tex table
print(dat
  , floating            = FALSE
  , latex.environments  = NULL
  , booktabs            = TRUE
  , include.rownames    = FALSE
)

# Note: Variance of some random terms is mereley 0. This shouldn't be an issue
# though -> Check this post: # https://stats.stackexchange.com/questions/115090/
# why-do-i-get-zero-variance-of-a-random-effect-in-my-mixed-model-despite-some-va

############################################################
#### Permeability Surface
############################################################
# Load required data
permeability  <- "03_Data/03_Results/99_PermeabilityMap.tif" %>%
  raster()
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>%
  readOGR() %>%
  as("SpatialLines")
africa <- "03_Data/02_CleanData/00_General_Africa.shp" %>%
  readOGR() %>%
  as("SpatialLines")
africa_crop <- "03_Data/02_CleanData/00_General_Africa.shp" %>%
  readOGR() %>%
  crop(kaza)

# Select only countries of interest
africa_crop <- subset(africa_crop, COUNTRY %in% c(
      "Angola"
    , "Namibia"
    , "Botswana"
    , "Zimbabwe"
    , "Zambia"
  )
)

# Prepare the extent of our core area
core <- as(extent(c(23.25, 24, -19.8, -19.05)), "SpatialLines")
crs(core) <- CRS("+init=epsg:4326")

# Rescale between 0 and 1 (necessary to get nice names)
permeability <- normalizeMap(permeability)

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

# We will indicate the Okavango delta and Linyanti riversystem
waters <- data.frame(
    x = c(22.368000, 23.611706)
  , y = c(-18.506023, -17.670367)
  , Name = c("Okavango\nDelta", "Linyanti\nSwamp")
)
coordinates(waters) <- c("x", "y")
crs(waters) <- CRS("+init=epsg:4326")

# Prepare lines that point from the text object to the water area
l1 <- cbind(c(22.468000, 22.668000), c(-18.806023, -19.106023)) %>% spLines()
l2 <- cbind(c(23.611706, 23.611706), c(-17.970367, -18.270367)) %>% spLines()
waters_pointers <- rbind(l1, l2)
crs(waters_pointers) <- crs(waters)

# In the second plot we will need different locations for the text
waters2 <- data.frame(
    x = c(22.668000, 23.611706)
  , y = c(-19.106023, -18.270367)
  , Name = c("Okavango\nDelta", "Linyanti\nSwamp")
)
coordinates(waters2) <- c("x", "y")
crs(waters2) <- CRS("+init=epsg:4326")

# Prepare a dataframe with some locations that we want to plot
villages <- data.frame(
    x = c(23.419332, 24.178607, 24.271072, 23.653653)
  , y = c(-19.999149, -19.430306, -19.192169, -19.002809)
  , Name = c("Maun", "Sankuyo", "Mababe", "Khwai")
)
coordinates(villages) <- c("x", "y")
crs(villages) <- CRS("+init=epsg:4326")

# Also prepare some pointers for the villages
l1 <- cbind(c(24.030000, 23.880000), c(-19.430306, -19.427306)) %>% spLines()
l2 <- cbind(c(24.135000, 24.020000), c(-19.195306, -19.188306)) %>% spLines()
l3 <- cbind(c(23.653653, 23.753653), c(-19.042809, -19.152809)) %>% spLines()
villages_pointers <- rbind(l1, l2, l3)
crs(villages_pointers) <- crs(villages)

# Prepare the plot of the entire region
p1 <- tm_shape(permeability) +
    tm_raster(
        palette         = "viridis"
      , style           = "cont"
      , title           = "Permeability"
      , labels          = c("Low", "", "High")
      , legend.reverse  = T
    ) +
  tm_grid(
      n.x = 5
    , n.y = 5
    , labels.inside.frame = FALSE
    , lines = FALSE
    , ticks = TRUE
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
  tm_shape(waters) +
    tm_text("Name"
    , col       = "white"
    , fontface  = 3
    , size      = 0.6
  ) +
  tm_shape(waters_pointers) +
    tm_lines(
      col = "white"
    , lwd = 0.7
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
      # , legend.title.size   = 0.7
      # , legend.text.size    = 0.5
      , legend.bg.color     = viridis(1, begin = 0)
  ) +
  # tm_add_legend(
  #     type    = "line"
  #   , labels  = "KAZA TFCA"
  #   , col     = "white"
  #   , lwd     = 2
  # ) +
  # tm_add_legend(
  #     type    = "line"
  #   , labels  = "Country Borders"
  #   , col     = "white"
  #   , lty     = 2
  #   , lwd     = 1
  # ) +
  tm_scale_bar(
        position    = "left"
      , text.size   = 0.5
      , text.color  = "white"
      , width = 0.125
  ) +
  tm_compass(
        text.color = "white"
      , color.dark = "white"
  ) #+
#   tm_credits("(a)"
#   , position  = c("right", "top")
#   , size      = 1.5
#   , col       = "white"
# )

# Prepare a plot of the Okavango Delta only
p2 <- tm_shape(permeability, bbox = bbox) +
    tm_raster(
        palette         = "viridis"
      , style           = "cont"
      , title           = "Permeability"
      , labels          = c("Low", "", "High")
      , legend.reverse  = T
    ) +
  tm_shape(core) +
    tm_lines(
        col = "white"
      , lty = 2
      , lwd = 2
    ) +
  tm_shape(villages) +
    tm_text("Name"
    , col   = "white"
    , size  = 0.9
  ) +
  tm_shape(villages_pointers) +
    tm_lines(
      col = "white"
    , lwd = 0.7
  ) +
  tm_shape(waters2) +
    tm_text("Name"
    , col       = "white"
    , fontface  = 3
  ) +
  tm_grid(
      n.x = 5
    , n.y = 5
    , labels.inside.frame = FALSE
    , lines = FALSE
    , ticks = TRUE
  ) +
  tm_layout(
        legend.position     = c("left", "top")
      , legend.text.color   = "white"
      , legend.title.color  = "white"
      # , legend.title.size   = 0.7
      # , legend.text.size    = 0.5
      , legend.bg.color     = viridis(1, begin = 0)
  ) +
  # tm_add_legend(
  #     type    = "line"
  #   , labels  = "Core Study Area"
  #   , col     = "white"
  #   , lty     = 1
  #   , lwd     = 1
  # ) +
  tm_scale_bar(
        position    = "left"
      , text.size   = 0.5
      , text.color  = "white"
      , width = 0.125
  ) +
  tm_compass(
        text.color = "white"
      , color.dark = "white"
  ) +
  tm_credits("(b)"
  , position  = c("right", "top")
  , size      = 1.5
  , col       = "white"
)

# Store the plots
CairoPDF("04_Manuscript/99_PermeabilityMap.pdf"
  , width   = 6
  , height  = 5.25
)
p1
dev.off()

CairoPDF(
    "04_Manuscript/99_PermeabilityMap2.pdf"
  , width   = 6
  , height  = 5.25
)
p2
dev.off()

############################################################
#### Permeability Comparison
############################################################
# Load required data
count_kaza <- read_rds("03_Data/03_Results/99_PermeabilityValues(KAZA).rds")
count_prot <- read_rds("03_Data/03_Results/99_PermeabilityValues(Prot).rds")

# Prepare tables
dat1 <- count_kaza %>%
  group_by(., Country, Group) %>%
  summarize(.
    , Median  = format(round(median(Permeability), 2), nsmall = 2)
    , IQR     = format(round(IQR(Permeability), 2), nsmall = 2)
  ) %>%
  ungroup() %>%
  mutate(Permeability = paste0(Median, " (", IQR, ")")) %>%
  select(-c(Median, IQR)) %>%
  spread(key = Group, value = Permeability) %>%
  setNames(c("Country", "Inside", "Outside"))

dat2 <- count_prot %>%
  group_by(., Country, Group) %>%
  summarize(.
    , Median  = format(round(median(Permeability), 2), nsmall = 2)
    , IQR     = format(round(IQR(Permeability), 2), nsmall = 2)
  ) %>%
  ungroup() %>%
  mutate(Permeability = paste0(Median, " (", IQR, ")")) %>%
  select(-c(Median, IQR)) %>%
  spread(key = Group, value = Permeability) %>%
  select(c(Country, Protected, Pastoral))

# Put them together
dat <- cbind(dat1, dat2[, 2:3])

# Add sums over countries
dat <- count_kaza %>%
  group_by(., Country) %>%
  summarize(.
    , Median  = format(round(median(Permeability), 2), nsmall = 2)
    , IQR     = format(round(IQR(Permeability), 2), nsmall = 2)
  ) %>%
  ungroup() %>%
  mutate(Permeability = paste0(Median, " (", IQR, ")")) %>%
  select(-c(Median, IQR)) %>%
  setNames(c("Country", "Overall")) %>%
  left_join(dat, ., by = "Country")

# Add sums over regions
dat1 <- count_kaza %>%
  group_by(Group) %>%
  summarize(.
    , Median  = format(round(median(Permeability), 2), nsmall = 2)
    , IQR     = format(round(IQR(Permeability), 2), nsmall = 2)
  ) %>%
  ungroup() %>%
  mutate(Permeability = paste0(Median, " (", IQR, ")")) %>%
  select(-c(Median, IQR)) %>%
  spread(key = Group, value = Permeability)

dat2 <- count_prot %>%
  group_by(Group) %>%
  summarize(.
    , Median  = format(round(median(Permeability), 2), nsmall = 2)
    , IQR     = format(round(IQR(Permeability), 2), nsmall = 2)
  ) %>%
  ungroup() %>%
  mutate(Permeability = paste0(Median, " (", IQR, ")")) %>%
  select(-c(Median, IQR)) %>%
  spread(key = Group, value = Permeability) %>%
  select(c(Protected, Pastoral))

dat3 <- data.frame(
    Median  = format(round(median(count_prot$Permeability), 2), nsmall = 2)
  , IQR     = format(round(IQR(count_prot$Permeability), 2), nsmall = 2)
  ) %>%
  mutate(Permeability = paste0(Median, " (", IQR, ")")) %>%
  select(-c(Median, IQR))

# Put the data into our dataframe
dat[nrow(dat) + 1, ] <- as.matrix(cbind("Overall", dat1, dat2, dat3))

# Create latex table with multiple groups
table <- kable(dat, "latex", booktabs = T, linesep = "") %>%
  kable_styling(latex_options = "hold_position") %>%
  add_header_above(c(" " = 1, "KAZA-TFCA" = 2, "Protection Status" = 2, " " = 1)) %>%
  row_spec(5, hline_after = T)

# Store the table to a latex file
write(table, "04_Manuscript/99_PermeabilityComparisons.tex")

############################################################
#### Plot Permeability Values Below Different Polygons
############################################################
# Reload the required data
dat_kaza <- read_rds("03_Data/03_Results/99_PermeabilityValues(KAZA).rds")
dat_nati <- read_rds("03_Data/03_Results/99_PermeabilityValues(Prot).rds")

# Prepare a new column for facetting
dat_kaza$Facet <- "Comparison 1"
dat_nati$Facet <- "Comparison 2"

# Put the two datasets together
dat <- rbind(dat_kaza, dat_nati)

# Prepare new labels that indicate the number of datapoints for each
# cateogry. We will use them as new axis labels
dat <- dat %>%
  group_by(., Group) %>%
  mutate(n = comma(n())) %>%
  mutate(Label = paste0(Group, "\n", "n=", n)) %>%
  transform(Label = as.factor(Label))

# Reorder factor levels
dat$Label <- factor(dat$Label, levels(dat$Label)[c(1, 3, 2, 5, 4)])

# Prepare a violin plot
p1 <- ggplot(data = dat, aes(x = Label, y = Permeability)) +

  # Add the violin
  geom_violin(width = 1) +

  # Add a boxplot inside the violin
  geom_boxplot(width = 0.1, color = "black", alpha = 0.2, outlier.alpha = 0) +

  # Use the classic theme
  theme_classic() +

  # Some beautifications
  theme(
      legend.position   = "none"
    , strip.background  = element_blank()
    , strip.text.x      = element_blank()
  ) +

  # Create facets by group
  facet_wrap(~ Facet, ncol = 1, scales = "free", strip.position = ) +

  # Change labs
  xlab("") +
  ylab("Permeability") +
  coord_flip() +

  # Add a title
  ggtitle("(b)") +

  # Change the size of the title
  theme(plot.title = element_text(size = 20))

# Store the plot
CairoPDF("04_Manuscript/99_PermeabilityComparisons.pdf", width = 8, height = 5)
p1
dev.off()

############################################################
#### Show How Least Cost Paths and Corridors Work
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

# Prepare Plot of source points
p1 <- tm_shape(norm(map)) +
    tm_raster(
        palette = "viridis"
      , style   = "cont"
      , title   = "Permeability"
      , labels  = c("Low", "", "High")
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
  ) +
  tm_credits("(a)"
  , position  = c("right", "top")
  , size      = 1.5
  , col       = "white"
)

# Plot the cost maps
p2 <- tm_shape(norm(A_cost)) +
    tm_raster(
        palette = "-magma"
      , style   = "cont"
      , title   = "Cumulative Cost"
      , labels  = c("Low", "", "High")
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
  ) +
  tm_credits("(b1)"
  , position  = c("right", "top")
  , size      = 1.5
  , col       = "white"
)

p3 <- tm_shape(norm(B_cost)) +
    tm_raster(
        palette = "-magma"
      , style   = "cont"
      , title   = "Cumulative Cost"
      , labels  = c("Low", "", "High")
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
  ) +
  tm_credits("(b2)"
  , position  = c("right", "top")
  , size      = 1.5
  , col       = "white"
)

p4 <- tm_shape(norm(cost)) +
    tm_raster(
        palette = "-magma"
      , style   = "cont"
      , title   = "Cumulative Cost"
      , labels  = c("Low", "", "High")
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
  ) +
  tm_credits("(c)"
  , position  = c("right", "top")
  , size      = 1.5
  , col       = "white"
)

p5 <- tm_shape(norm(corr)) +
    tm_raster(
        palette = "-magma"
      , style   = "cont"
      , title   = "Corridor"
      , labels  = c("Low-Cost", "", "High-Cost")
    ) +
  # tm_shape(path) +
  #   tm_lines(
  #     col = "black"
  #   ) +
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
  ) +
  tm_credits("(d)"
  , position  = c("right", "top")
  , size      = 1.5
  , col       = "white"
)

# Store the plots
CairoPDF("04_Manuscript/99_LeastCostPaths(Example1).pdf", width = 3, height = 3)
p1
dev.off()

CairoPDF("04_Manuscript/99_LeastCostPaths(Example2).pdf", width = 3, height = 6)
tmap_arrange(p2, p3, ncol = 1)
dev.off()

CairoPDF("04_Manuscript/99_LeastCostPaths(Example3).pdf", width = 3, height = 6)
tmap_arrange(p4, p5, ncol = 1)
dev.off()

############################################################
#### Least Cost Paths & Corridors
############################################################
# Load required data
perm        <- raster("03_Data/03_Results/99_PermeabilityMap.tif")
corrs       <- raster("03_Data/03_Results/99_LeastCostCorridors.tif")
paths       <- readOGR("03_Data/03_Results/99_LeastCostPaths.shp")
paths_rast  <- raster("03_Data/03_Results/99_LeastCostPaths.tif")
kaza        <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp") %>% as(., "SpatialLines")
africa      <- readOGR("03_Data/02_CleanData/00_General_Africa.shp") %>% as(., "SpatialLines")
africa_crop <- readOGR("03_Data/02_CleanData/00_General_Africa.shp") %>% crop(., kaza)
prot        <- readOGR("03_Data/03_Results/99_SourceAreas.shp")
points      <- readOGR("03_Data/03_Results/99_SourcePoints.shp")
nati        <- readOGR("03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(3Classes).shp")

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

# Replace zeros with NAs
values(paths_rast)[values(paths_rast) == 0] <- NA

# Rescale between 0 and 1
corrs <- calc(corrs, fun = function(x){
  (x - min(x)) / (max(x) - min(x))
})
paths_rast <- calc(paths_rast, fun = function(x){
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
})

# Prepare the plot of LCPs (we only need the corridors to make sure we get a
# nice plot border)
p1 <- tm_shape(corrs) +
    tm_raster(
      palette     = "white"
    , legend.show = FALSE
  ) +
  tm_shape(prot) +
    tm_polygons(
      col           = "grey65"
    , border.col    = "grey65"
  ) +
  tm_grid(
      n.x                 = 5
    , n.y                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_shape(paths_rast) +
    tm_raster(
        palette         = "-Spectral"
      , style           = "cont"
      , title           = "Least-Cost Paths"
      , labels          = c("Low-Frequency", "", "High-Frequency")
      , legend.reverse  = T
    ) +
  tm_shape(points) +
    tm_dots(
        col   = "black"
      , size  = 0.1
    ) +
  tm_shape(kaza) +
    tm_lines(
        col = "black"
      , lwd = 2
    ) +
  tm_shape(africa) +
    tm_lines(
        col = "black"
      , lwd = 1
      , lty = 2
    ) +
  tm_shape(africa_crop) +
    tm_text("COUNTRY"
      , col       = "black"
      , just      = "bottom"
      , fontface  = 2
    ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "black"
    , width       = 0.125
  ) +
  tm_compass(
      text.color    = "black"
    , color.dark    = "black"
    , color.light   = "black"
  ) +
  tm_credits("(a)"
  , position  = c("right", "top")
  , size      = 1.5
  , col       = "black"
)

# Store the plot
CairoPDF("04_Manuscript/99_LeastCostPaths.pdf", width = 6, height = 5.25)
p1
dev.off()

# Prepare the plot for LCCs
p2 <- tm_shape(corrs) +
    tm_raster(
        palette         = "-Spectral"
      , style           = "cont"
      , title           = "Least-Cost Corridors"
      , labels          = c("Low-Frequency", "", "High-Frequency")
      , legend.reverse  = T
    ) +
  tm_grid(
      n.x                 = 5
    , n.y                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_shape(nati) +
    tm_borders(
        col   = "black"
      , alpha = 0.6
    ) +
  tm_shape(kaza) +
    tm_lines(
        col = "black"
      , lwd = 2
    ) +
  tm_shape(nati_text) +
    tm_text("Name"
      , col       = "black"
      , alpha     = 0.6
      , fontface  = 3
      , size      = 0.5
      , shadow    = TRUE
    ) +
  tm_shape(africa) +
    tm_lines(
        col = "black"
      , lwd = 1
      , lty = 2
    ) +
  tm_shape(africa_crop) +
    tm_text("COUNTRY"
      , col       = "black"
      , just      = "bottom"
      , fontface  = 2
    ) +
  tm_layout(
      legend.text.color   = "white"
    , legend.title.color  = "white"
  ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color   = "white"
    , color.dark   = "white"
    , color.light  = "white"
  ) +
  tm_credits("(b)"
  , position  = c("right", "top")
  , size      = 1.5
  , col       = "white"
)

# Store the plot
CairoPDF("04_Manuscript/99_LeastCostCorrs.pdf", width = 6, height = 5.25)
p2
dev.off()

############################################################
#### Simulations
############################################################
# Load required data
paths       <- raster("03_Data/03_Results/99_DispersalSimulation.tif")
kaza        <- shapefile("03_Data/02_CleanData/00_General_KAZA_KAZA") %>% as(., "SpatialLines")
africa      <- shapefile("03_Data/02_CleanData/00_General_Africa") %>% as(., "SpatialLines")
africa_crop <- shapefile("03_Data/02_CleanData/00_General_Africa") %>% crop(., kaza)
prot        <- shapefile("03_Data/03_Results/99_SourceAreas.shp")
points      <- shapefile("03_Data/03_Results/99_SourcePoints.shp")
nati        <- shapefile("03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(3Classes)")

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

# Take sqrt of LCPs to make lower frequented LPCs more visible
paths <- sqrt(paths)

# Rescale between 0 and 1
paths <- calc(paths, fun = function(x){
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
})

# Plot simulated paths
p1 <- tm_shape(paths) +
    tm_raster(
        palette = "-Spectral"
      , style   = "cont"
      , title   = "Simulated\nTrajectories"
      , labels  = c("Low-Frequency", "", "High-Frequency")
    ) +
  tm_shape(points) +
    tm_dots(
        col   = viridis(20)[20]
      , size  = 0.1
    ) +
  tm_grid(
      n.x                 = 5
    , n.y                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_shape(nati) +
    tm_borders(
        col   = "white"
      , alpha = 0.6
    ) +
  tm_shape(kaza) +
    tm_lines(
        col = "white"
      , lwd = 2
    ) +
  tm_shape(nati_text) +
    tm_text("Name"
      , col       = "white"
      , alpha     = 0.6
      , fontface  = 3
      , size      = 0.5
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
  tm_layout(
      legend.text.color   = "white"
    , legend.title.color  = "white"
  ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color   = "white"
    , color.dark   = "white"
    , color.light  = "white"
)

# Store the plots
CairoPDF("04_Manuscript/99_Simulations.pdf", width = 6, height = 6.3)
p1
dev.off()

############################################################
#### Yearly Floodmap Cycle
############################################################
# Load required data
flood <- stack(
    "03_Data/02_CleanData/00_Floodmaps/01_Original/2018.03.14.tif"
  , "03_Data/02_CleanData/00_Floodmaps/01_Original/2018.04.15.tif"
  , "03_Data/02_CleanData/00_Floodmaps/01_Original/2018.05.17.tif"
  , "03_Data/02_CleanData/00_Floodmaps/01_Original/2018.06.18.tif"
  , "03_Data/02_CleanData/00_Floodmaps/01_Original/2018.07.20.tif"
  , "03_Data/02_CleanData/00_Floodmaps/01_Original/2018.08.21.tif"
  , "03_Data/02_CleanData/00_Floodmaps/01_Original/2018.09.22.tif"
  , "03_Data/02_CleanData/00_Floodmaps/01_Original/2018.10.24.tif"
  , "03_Data/02_CleanData/00_Floodmaps/01_Original/2018.11.25.tif"
)

# Get rid of clouds
rcl <- data.frame(old = c(0, 127, 255), new = c(1, 0, 0))
flood <- mcreclassify(flood, rcl)

# Get titles for images
names <- substr(names(flood), start = 2, stop = 11)
names <- gsub(names, pattern = "\\.", replacement = "-")

# Prepare plots
p <- list()
for (i in 1:nlayers(flood)){
  p[[i]] <- tm_shape(flood[[i]]) +
    tm_raster(
        palette = c("white", "blue")
      , style   = "cat"
      , legend.show = FALSE
    ) +
    tm_layout(
      , title = names[i]
      , title.position = c("left", "bottom")
    )
}

pc <- tmap_arrange(
    p[[1]]
  , p[[2]]
  , p[[3]]
  , p[[4]]
  , p[[5]]
  , p[[6]]
  , p[[7]]
  , p[[8]]
  , p[[9]]
)

CairoPDF("04_Manuscript/99_FloodPulse.pdf", width = 6, height = 6)
pc
dev.off()

# ############################################################
# #### Movement Model
# ############################################################
# # Load required data
# move <- readRDS("99_MovementModel.rds")
#
# # Select best model
# best <- move$Model[[1]]
#
# # Get coefficients from best model
# coeffs <- getCoeffs(best)[-1, ]
# coeffs$Covariate[coeffs$Covariate == "Shrubs"] <- "Shrubs/Grassland"
#
# # Prepare a table that shows the results
# best_table <- getCoeffs(best, pvalue = TRUE) %>%
#
#   # Calculate confidence intervals
#   mutate(
#       LCI = Coefficient - 1.96 * SE
#     , UCI = Coefficient + 1.96 * SE
#   )
#
# # Remove CIs
# best_table <- select(best_table, -c("LCI", "UCI"))
#
# # Write the table to a .tex table
# print(xtable(best_table, digits = 3)
#   , floating            = FALSE
#   , latex.environments  = NULL
#   , booktabs            = TRUE
#   , include.rownames    = FALSE
#   , type                = "latex"
# )
#
# # Prepare plot
# p <- showCoeffs(coeffs
#   , shape     = 1
#   , size      = 1
#   , whiskers  = 1.96
#   , order     = c(
#       "cos(ta_)"
#     , "log(sl_)"
#     , "Water"
#     , "DistanceToWater"
#     , "Trees"
#     , "Shrubs/Grassland"
#     , "HumansBuff5000"
#     , "log(sl_):Water"
#     , "cos(ta_):DistanceToWater"
#     , "cos(ta_):HumansBuff5000"
# ))
# p <- p + ggtitle("(a)") + theme(plot.title = element_text(size = 10, face = "bold"))
#
# # Store the plot
# setwd(output)
# CairoPDF("99_MovementModel.pdf", width = 5, height = 3)
# p
# dev.off()
# setwd(input)
#
# # Prepare raster for contourplot of interactions
# r1 <- visInt2(best, "log(sl_)", "Water")
# r2 <- visInt2(best, "cos(ta_)", "DistanceToWater")
# r3 <- visInt2(best, "cos(ta_)", "HumansBuff5000")
#
# # Prepare the plots
# p1 <- levelplot(r1
#   , contour = TRUE
#   , xlab    = "log(sl_)"
#   , ylab    = "Water"
#   , margin  = FALSE
#   , main    = "(b1)"
#   , par.settings = viridisTheme
# )
#
# p2 <- levelplot(r2
#   , contour = TRUE
#   , xlab    = "cos(ta_)"
#   , ylab    = "DistanceToWater"
#   , margin  = FALSE
#   , main    = "(b2)"
#   , par.settings = viridisTheme
# )
#
# p3 <- levelplot(r3
#   , contour = TRUE
#   , xlab    = "cos(ta_)"
#   , ylab    = "HumansBuff5000"
#   , margin  = FALSE
#   , main    = "(b3)"
#   , par.settings = viridisTheme
# )
#
# # Put the plots together and store them
# setwd(output)
# CairoPDF("99_MovementModel(Interactions).pdf"
#   , width     = 14/1.5
#   , height    = 5/1.5
#   , pointsize = 40
# )
# grid.arrange(p1, p2, p3, ncol = 3)
# dev.off()
# setwd(input)
#
# ############################################################
# #### Movement Model AICs
# ############################################################
# # Load required data
# models <- read_rds("99_MovementModel.rds")
#
# # Let's prepare a nice table to report in our results
# table <- cbind(
#     models$Covariates
#   , select(models, -c(Covariates, Formula, Model))
# )
# aics <- table
#
# # Remove any entries with either "HumansAverage" or "HumansBase"
# aics <- aics[!grepl(pattern = "HumansAverage|HumansBase", aics$Covariates), ]
#
# # Also remove the ModelID column
# aics$ModelID <- NULL
#
# # Rename the shrubs
# aics$Covariates <- str_replace(
#     string      = aics$Covariates
#   , pattern     = "Shrubs"
#   , replacement = "Shrubs/Grassland"
# )
#
# # Write the table to a .tex table
# print(xtable(aics)
#   , floating            = FALSE
#   , latex.environments  = NULL
#   , booktabs            = TRUE
#   , include.rownames    = FALSE
#   , type                = "latex"
# )
#
# ############################################################
# #### Social Model
# ############################################################
# # Load required data
# social <- readRDS("99_SocialModel.rds")
#
# # Select best model
# best <- social
#
# # Get coefficients from best model
# coeffs <- getCoeffs(best)[-1, ]
# coeffs$Covariate[coeffs$Covariate == "Shrubs"] <- "Shrubs/Grassland"
#
# # Prepare a table that shows the results
# best_table <- getCoeffs(best, pvalue = TRUE) %>%
#
#   # Calculate confidence intervals
#   mutate(
#       LCI = Coefficient - 1.96 * SE
#     , UCI = Coefficient + 1.96 * SE
#   )
#
# # Remove CIs
# best_table <- select(best_table, -c("LCI", "UCI"))
#
# # Write the table to a .tex table
# print(xtable(best_table, digits = 3)
#   , floating            = FALSE
#   , latex.environments  = NULL
#   , booktabs            = TRUE
#   , include.rownames    = FALSE
#   , type                = "latex"
# )
#
# # Prepare plot
# p <- showCoeffs(coeffs
#   , shape     = 1
#   , size      = 1
#   , whiskers  = 1.96
#   , xlim      = c(-1.75, 1.75)
#   , order     = c(
#       "cos(ta_)"
#     , "log(sl_)"
#     , "Water"
#     , "DistanceToWater"
#     , "Trees"
#     , "Shrubs/Grassland"
#     , "HumansBuff5000"
#     , "ProbEncounter"
#     , "OwnHomerange"
#     , "ForeignHomeranges"
#     , "ProbEncounter:OwnHomerange"
#     , "ProbEncounter:ForeignHomeranges"
# ))
# p <- p + ggtitle("(a)") + theme(plot.title = element_text(size = 10, face = "bold"))
#
# # Store the plot
# setwd(output)
# CairoPDF("99_SocialModel.pdf", width = 6, height = 3)
# p
# dev.off()
# setwd(input)
#
# # Prepare raster for contourplot of interactions
# r1 <- visInt2(best, "ProbEncounter", "OwnHomerange")
# r2 <- visInt2(best, "ProbEncounter", "ForeignHomeranges")
#
# # Prepare the plots
# p1 <- levelplot(r1
#   , contour = TRUE
#   , xlab    = "ProbEncounter"
#   , ylab    = "OwnHomerange"
#   , margin  = FALSE
#   , main    = "(b1)"
#   , par.settings = viridisTheme
# )
#
# p2 <- levelplot(r2
#   , contour = TRUE
#   , xlab    = "ProbEncounter"
#   , ylab    = "ForeignHomeranges"
#   , margin  = FALSE
#   , main    = "(b2)"
#   , par.settings = viridisTheme
# )
#
# # Put the plots together and store them
# setwd(output)
# CairoPDF("99_SocialModel(Interactions).pdf", width = 14/3*2, height = 5)
# grid.arrange(p1, p2, ncol = 2)
# dev.off()
# setwd(input)
#
# ############################################################
# #### Social Model: Example Landscape
# ############################################################
# # load required data
# setwd("/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/01_SmallExtent")
# uds <- stack("05_SocialFeatures_SocialLandscape(UtilisationDistributions).grd")
# hrs <- shapefile("05_SocialFeatures_SocialLandscape(HomeRanges)")
#
# # Select a date to plot
# names(uds)
# hrs_sub <- subset(hrs, Date == "2017-07-01")
# uds_sub <- crop(uds[["X2017.07.01"]], hrs_sub)
#
# # Prepare plot of the uds
# p1 <- tm_shape(uds_sub) +
#   tm_raster(
#       palette = "-Spectral"
#     , style = "cont"
#     , title = "Probability of Encounter"
#     , labels = c("Low", "", "High")
#   ) +
#   tm_layout(
#       , legend.text.color       = "white"
#       , legend.title.color      = "white"
#       , legend.position = c("left", "bottom")
#   ) +
#   tm_compass(
#       text.color    = "white"
#     , color.dark    = "white"
#     , color.light   = "white"
#   ) +
#   tm_scale_bar(
#       position    = "right"
#     , text.size   = 0.5
#     , text.color  = "white"
#     , width       = 0.125
#   ) +
#   tm_credits("(a)"
#   , position  = c("right", "top")
#   , size      = 1.5
#   , col       = "white"
# )
#
# # Prepare a plot of the HRs
# p2 <- tm_shape(uds_sub) +
#     tm_raster(
#       palette     = "white"
#     , legend.show = FALSE
#   ) +
#   tm_shape(hrs_sub) +
#     tm_polygons(
#       "Pack"
#     , palette = "viridis"
#     , alpha = 0.5
#   ) +
#   tm_layout(
#       , legend.text.color       = "black"
#       , legend.title.color      = "black"
#       , legend.position = c("left", "bottom")
#   ) +
#   tm_compass(
#       text.color    = "black"
#     , color.dark    = "black"
#     , color.light   = "black"
#   ) +
#   tm_scale_bar(
#       position    = "right"
#     , text.size   = 0.5
#     , text.color  = "black"
#     , width       = 0.125
#   ) +
#   tm_credits("(b)"
#   , position  = c("right", "top")
#   , size      = 1.5
#   , col       = "black"
# )
#
# setwd(output)
# CairoPDF("99_SocialLandscapeExample.pdf", width = 7.5, height = 4)
# tmap_arrange(p1, p2, ncol = 2)
# dev.off()
# setwd(output)
