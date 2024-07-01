################################################################################
#### Plot of Study Area
################################################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(rgdal)        # To load and store spatial data
library(raster)       # To manipulate raster data
library(rgeos)        # To manipulate spatial data
library(tmap)         # For beautiful maps
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(grid)         # To arrange multiple plots
library(Cairo)        # To store plots

################################################################################
#### Colors
################################################################################
# Create a color ramp
pal1 <- colorRampPalette(colors = c("black", "orange", "yellow"))

# Visualize colors
plot(1:20, 1:20, pty = 20, pch = 20, cex = 15, col = pal1(20))

################################################################################
#### Plot of the Study Area
################################################################################
# Load required data
africa <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
kaza <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
dogs <- readOGR("03_Data/02_CleanData/00_General_WildDogs_IUCN.shp")
prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
water <- readOGR("03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters_GEOFABRIK.shp")

# Remove 0 from water layer
water2[water2 == 0] <- NA

# Simplify Protection zones
prot$Desig[prot$Desig == "Forest Reserve"] <- "Protected"
prot$Desig <- as.factor(as.character(prot$Desig))

# Remove small islands for plotting
africa <- subset(africa, !(ID %in% c(27:41, 689, 690)))

# Create a buffered polygon of africa
africa2 <- gBuffer(africa, width = metersToDegrees(100))

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Create labels for countries
labels_countries <- data.frame(
    x = c(16, 24, 11, 38, 34)
  , y = c(-9, -30, -26, -13, -23)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Create lines pointing towards these countries
l1 <- rbind(c(17, -10), c(19, -12)) %>% spLines()
l2 <- rbind(c(24, -29), c(24, -23)) %>% spLines()
l3 <- rbind(c(13, -25), c(17, -22)) %>% spLines()
l4 <- rbind(c(35, -13), c(30, -14)) %>% spLines()
l5 <- rbind(c(33, -22), c(30, -19)) %>% spLines()
lines_countries <- rbind(l1, l2, l3, l4, l5)
crs(lines_countries) <- crs(labels_countries)

# Prepare a plot of the KAZA. Create labels for countries first.
labels_countries2 <- data.frame(
    x = c(20.39, 23.94, 20.07, 25.69, 28.22)
  , y = c(-15.28, -19.94, -19.39, -15.22, -18.9)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries2) <- c("x", "y")
crs(labels_countries2) <- CRS("+init=epsg:4326")

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


# Prepare a map of Africa
p1 <- tm_shape(africa2) +
    tm_polygons(col = "gray70", border.col = "gray70", lwd = 2) +
  # tm_grid(
  #     n.y                 = 5
  #   , n.x                 = 5
  #   , labels.inside.frame = false
  #   , lines               = true
  #   , ticks               = true
  #   , col = "gray12"
  # ) +
  tm_shape(africa) +
    tm_polygons(
        col = "gray40"
      , lwd = 0.7
      , border.col = "gray70"
    ) +
  tm_shape(dogs) +
    tm_polygons(
        col           = "orange"
      , alpha         = 0.8
      , border.alpha  = 0
    ) +
  tm_shape(kaza) +
    tm_borders(
        col = "white"
      , lty = 1
      , lwd = 1
    ) +
  tm_shape(kaza_ext) +
    tm_borders(
        col = pal1(20)[10]
      , lty = 3
      , lwd = 1.5
    ) +
  tm_shape(lines_countries) +
    tm_lines(
        col = "white"
    ) +
  tm_shape(labels_countries) +
    tm_text(
        "Label"
      , size = 0.6
      , col = "white"
    ) +
  tm_layout(
      asp         = 0.8
    , frame       = "black"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.color    = "black"
  ) +
  tm_credits("(a)"
    , position  = c("right", "top")
    , size      = 1.5
    , col       = "white"
)

# Prepare a map of Kaza
p2 <- tm_shape(prot) +
    tm_polygons(
      , col           = "gray20"
      , border.col    = "black"
      , border.alpha  = 0.8
      , lwd           = 0.5
      , legend.show   = F
    ) +
  tm_shape(water) +
    tm_polygons(
        col           = "cornflowerblue"
      , border.col    = "cornflowerblue"
      , border.alpha  = 0.6
      , lwd           = 0.2
    ) +
  tm_shape(kaza, is.master = T) +
    tm_borders(
        col = "white"
      , lty = 1
      , lwd = 2
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray70"
    ) +
  tm_shape(labels_countries2) +
    tm_text("Label"
      , col       = "gray70"
      , fontface  = 2
      , size      = 1.5
    ) +
  tm_shape(labels_waters) +
    tm_text("Label"
      , fontface  = 3
      , size      = 0.5
      , col       = "white"
    ) +
  tm_shape(labels_nationalparks) +
    tm_text("Label"
      , col       = "gray70"
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
    , bg.color = "black"
    , frame                   = "gray20"
    , frame.lwd               = 3
    , asp                     = 1.2
    , legend.outside          = TRUE
    , legend.outside.position = "left"
    , legend.text.color       = "white"
    , legend.title.color      = "black"
    , legend.stack            = "vertical"
    , legend.text.size        = 0.8
    , legend.bg.color         = "black"
  ) +
  tm_scale_bar(
        position  = c("right", "bottom")
      , text.size = 0.5
      , text.col  = "white"
      , width     = 0.125
  ) +
  tm_credits("(b)"
    , position  = c("left", "top")
    , size      = 1.5
    , col       = "white"
  ) +
  tm_compass(
      color.dark  = "white"
    , color.light = "white"
    , text.color  = "white"
    , position    = c("left", "bottom")
)

################################################################################
#### Prepare a common Legend
################################################################################
# We now prepare a legend that the two plots share. However, we have to assign
# to one of the two plots. Let's use p2 for this
p2 <- p2 + tm_add_legend(
    type    = "fill"
  , labels  = c(
      "Wild Dog Populations"
    , "Major Water Areas"
    , "Protected Areas"
  )
  , col = c(
      "orange"
    , "cornflowerblue"
    , "gray20"
  )) + tm_add_legend(
    type    = "line"
  , labels  = c("KAZA-TFCA Borders", "Country Borders")
  , col     = c("white", "gray70")
  , lty     = c(1, 1)
  , lwd     = c(2, 2)
) + tm_layout(legend.frame.lwd = 2, legend.text.size = 1.05)

# Combine plots and store them
CairoPDF("StudyArea.pdf", width = 9, height = 5.25, bg = "black")
p2
print(p1, vp = viewport(0.158, 0.341, width = 0.56, height = 0.56))
dev.off()
