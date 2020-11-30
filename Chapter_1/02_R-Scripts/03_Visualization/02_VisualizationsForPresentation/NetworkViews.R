################################################################################
#### Two Possible Network Views
################################################################################
# Description: Simple plot of a random network that I can put in the background
# of a presentation

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(igraph) # For network analysis
library(raster) # For spatial data

################################################################################
#### Two Network Views
################################################################################
# Old directory
wd <- "/media/david/My Passport/Backups/WildDogs/15. PhD/00_WildDogs"
setwd(wd)

# Load required data
areas  <- shapefile("03_Data/03_Results/99_SourceAreas2.shp")
points <- shapefile("03_Data/03_Results/99_SourcePoints2.shp")
kaza <- shapefile("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")

# Option 1: Source point view
net1 <- make_empty_graph(n = nrow(points))
lay1  <- coordinates(points)

# Option 2: Raster view
r <- raster("03_Data/02_CleanData/00_General_Raster250.tif")
r <- aggregate(r, fact = 50000 / 250, fun = "max")

# Coerce raster to network
net2 <- make_empty_graph(n = ncell(r))
lay2  <- as.matrix(as.data.frame(r, xy = T)[, c(1, 2)])

# Coerce raster to shapefile
r <- as(r, "SpatialPolygons")

# Visualize Option 1
png("test1.png", width = 1080, height = 1080, bg = "transparent")
plot(areas, col = "gray30", border = "gray50", bg = "transparent", lwd = 3)
plot(kaza, add = T, border = "gray80", lwd = 5)
plot(
    net1
  , layout       = lay1
  , vertex.size  = 18
  , vertex.label = NA
  , add          = T
  , rescale      = F
)
dev.off()

# Visualize Option 2
png("test2.png", width = 1080, height = 1080, bg = "transparent")
plot(r, col = "gray30", border = "gray50", bg = "transparent", lwd = 3)
plot(kaza, add = T, border = "gray80", lwd = 5)
plot(
    net2
  , layout       = lay2
  , vertex.size  = 15
  , vertex.label = NA
  , add          = T
  , rescale      = F
)
dev.off()
