################################################################################
#### Plot of Betweenness
################################################################################
# Description: Preliminary results plot of the betweenness in our study area

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/media/david/My Passport/Backups/WildDogs/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(parallel)       # To run on multiple cores
library(pbmcapply)      # To run on multiple cores with progress bar
library(igraph)         # For network analysis
library(viridis)        # For nice colors
library(RColorBrewer)   # For more nice colors
library(davidoff)       # Custom functions
library(rgdal)
library(tmap)
library(rgeos)
library(Cairo)

# Set a seed
set.seed(12345)

################################################################################
#### Load and Clean Data
################################################################################
# Load the simulated dispersal trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Subset to simulations of interest
sims <- subset(sims
  , StepNumber    <= 500
  & PointSampling == "Static"
)

# Load the reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Coarsen to 5km
r <- aggregate(r, fact = 10000 / 250)

# Prepare vertices
vertices  <- 1:ncell(r)

# Store layout
lay <- as.matrix(as.data.frame(r, xy = T)[, c(1, 2)])

# Fill the rasters with unique cell values (we will use the values as cell IDs)
values(r)  <- vertices

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# At each coordinate we now extract the cell IDs
visits <- data.frame(
    ID      = sims$ID
  , x       = coordinates(sims)[, 1]
  , y       = coordinates(sims)[, 2]
  , Raster  = raster::extract(r, sims)
)

# We now want to create a visitation history for each trajectory. This history
# depicts all transitions from one raster cell to another. To do this
# repeatedly, we write a function to retrieve the visitation history from a
# sequence of values
visitHist <- function(x, singlecount = F){
  transitions <- data.frame(from = lag(x), to = x) %>%
    group_by(from, to) %>%
    na.omit() %>%
    summarize(TotalConnections = n(), .groups = "drop")
  if (singlecount){
    transitions$TotalConnections = 1
  }
  return(transitions)
}

# Let's check what the function does exactly
visitHist(c(1, 2, 2, 3, 4, 5, 1, 2), singlecount = F)
visitHist(c(1, 2, 2, 3, 4, 5, 1, 2), singlecount = T)

# We want to retrieve the visitation history to each trajectory seperately, so
# let's nest them.
visits <- visits %>% nest(data = -ID)

# Apply the function to each trajectory for each raster. We therefore identify
# the visitation history of each trajectory for each of the different spatial
# resolutions.
visits$History <- pbmclapply(1:nrow(visits)
  , mc.cores            = detectCores() - 1
  , ignore.interactive  = T
  , function(x){
    visitHist(visits$data[[x]]$Raster, singlecount = T)
  })

################################################################################
#### Calculate Network Metrics
################################################################################
# We will create and use networks in order to derive network metrics that we can
# plot spatially. Let's prepare a function that we can use to derive all desired
# network metrics and depict them on a raster.
netMet <- function(
    network   = NULL    # The network based on which the metrics are calculated
  , raster    = NULL    # The raster onto which the metrics are calculated
  , tempfile  = F       # Should the resulting raster go to a temporary file?
  , metrics = c("betweenness", "closeness", "degree")
  ){

    # Calculate the desired network metrics and return them as rasters
    result <- vector(mode = "list", length = 3)
    if ("betweenness" %in% metrics){
      betweenness <- raster
      values(betweenness) <- betweenness(network)
      names(betweenness) <- "betweenness"
      result[[1]] <- betweenness
    }
    if ("closeness" %in% metrics){
      closeness <- raster
      values(closeness) <- closeness(network)
      names(closeness) <- "closeness"
      result[[2]] <- closeness
    }
    if ("degree" %in% metrics){
      degree <- raster
      values(degree) <- degree(network)
      names(degree) <- "degree"
      result[[3]] <- degree
    }

    # Remove NULLs from the lsit
    result <- plyr::compact(result)

    # Put all into a stack
    result <- stack(result)

    # In case the user desires to store a temporary file, do so
    if (tempfile){
      result <- writeRaster(result, tempfile())
    }

    # Return the final stack
    return(result)
}

# Function to combine metrics from multiple graphs
stackMet <- function(
    metrics   = NULL                  # The stack of metrics
  , metric    = "betweenness"         # The metric to be combined
  , tempfile  = F
  , fun       = function(x){mean(x)}  # The function to apply to stack the metric
  ){
    stacked <- stack(lapply(metrics, function(z){z[[metric]]}))
    funned  <- calc(stacked, fun)
    if (tempfile){
      funned <- writeRaster(funned, tempfile())
    }
    return(funned)
}

# Prepare data
visits_all_ww2 <- tibble(
    History  = list(
      do.call(rbind, visits$History) %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections)) %>%
      ungroup() %>%
      mutate(weight = mean(TotalConnections) / TotalConnections)
  )
)

# Create graphs
net_ww2 <- graph_from_data_frame(visits_all_ww2$History[[1]], vertices = vertices)

# Check if weighted
is_weighted(net_ww2)

# Check if connected
is_connected(net_ww2)

# Calculate metrics
res_ww2 <- netMet(network = net_ww2, raster = r)

# Plot the result
plot(sqrt(res_ww2[["betweenness"]]), col = viridis(20), main = "With Weights fun 2")

# Load data for plotting
setwd("03_Data/02_CleanData")
betw <- res_ww2[["betweenness"]]
kaza        <- shapefile("00_General_KAZA_KAZA") %>% as(., "SpatialLines")
africa      <- shapefile("00_General_Africa") %>% as(., "SpatialLines")
africa_crop <- shapefile("00_General_Africa") %>% crop(., kaza)
prot        <- shapefile("02_LandUseTypes_Protected_PeaceParks(1Class)")
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

# Rescale betweenness between 0 and 1
betw <- calc(betw, fun = function(x){
  (x - min(x)) / (max(x) - min(x))
})

# Prepare a plot of the betweenness score only
p1 <- tm_shape(raster(betw)) +
    tm_raster(
      palette     = "black"
    , legend.show = FALSE
  ) +
  tm_shape(prot) +
    tm_polygons(
      col           = "black"
    , border.col    = "black"
  ) +
  tm_shape(sqrt(betw)) +
    tm_raster(
        palette        = magma(50)
      , style          = "cont"
      , title          = "Betweenness Score"
      , labels         = c("Low", "", "High")
      , legend.reverse = T
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
  ) +
  tm_layout(
      bg.color            = "black"
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = "transparent"
    , legend.position     = c("left", "top")
)

# Prepare a plot of the betweenness score and some shapes
p2 <- tm_shape(raster(betw)) +
    tm_raster(
      palette     = "black"
    , legend.show = FALSE
  ) +
  tm_shape(prot) +
    tm_polygons(
      col           = "black"
    , border.col    = "black"
  ) +
  tm_shape(sqrt(betw)) +
    tm_raster(
        palette        = magma(50)
      , style          = "cont"
      , title          = "Betweenness Score"
      , labels         = c("Low", "", "High")
      , legend.reverse = T
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
  ) +
  tm_layout(
      bg.color            = "black"
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = "transparent"
    , legend.position     = c("left", "top")
)

# Prepare a plot of the betweenness score and some shapes
p3 <- tm_shape(raster(betw)) +
    tm_raster(
      palette     = "black"
    , legend.show = FALSE
  ) +
  tm_shape(prot) +
    tm_polygons(
      col           = "black"
    , border.col    = "black"
  ) +
  tm_shape(sqrt(betw)) +
    tm_raster(
        palette        = magma(50)
      , style          = "cont"
      , title          = "Betweenness Score"
      , labels         = c("Low", "", "High")
      , legend.reverse = T
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
  tm_shape(nati_text) +
    tm_text("Name"
      , col       = "white"
      , alpha     = 0.5
      , fontface  = 3
      , size      = 0.5
    ) +
  tm_compass(
      text.color  = "white"
    , color.dark  = "white"
  ) +
  tm_layout(
      bg.color            = "black"
    , legend.text.color   = "white"
    , legend.title.color  = "white"
    , legend.bg.color     = "transparent"
    , legend.position     = c("left", "top")
)

png("test1.png"
  , width = 1080 * 4 / 3.5
  , height = 1080
  , bg = "transparent"
)
p1 + tm_layout(scale = 2.5)
dev.off()

png("test2.png"
  , width = 1080 * 4 / 3.5
  , height = 1080
  , bg = "transparent"
)
p2 + tm_layout(scale = 2.5)
dev.off()

png("test3.png"
  , width = 1080 * 4 / 3.5
  , height = 1080
  , bg = "transparent"
)
p3 + tm_layout(scale = 2.5)
dev.off()
