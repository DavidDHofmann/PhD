################################################################################
#### Animation of Dispersal from Moremi
################################################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(raster)
library(rgdal)
library(davidoff)
library(pbmcapply)
library(spatstat)
library(maptools)
library(viridis)
library(animation)
library(tmap)

################################################################################
#### Prepare Data
################################################################################
# Load simulated dispersal data
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Load protected areas
prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Create SpatialPoints for first location of each trajectory
first <- sims %>%
  dplyr::select("x", "y", "TrackID") %>%
  group_by(TrackID) %>%
  slice(1) %>%
  SpatialPointsDataFrame(
      coords      = cbind(.[["x"]], .[["y"]])
    , proj4string = CRS("+init=epsg:4326")
  )

# Assess from which protected area each trajectory left
first$Origin <- as.character(over(first, prot)$Name)
first <- first@data[, c("TrackID", "Origin")]

# Join information to simulated tracks
sims <- left_join(sims, first, by = "TrackID")
sims$Origin[sims$Area == "Buffer"] <- "Buffer"

# Remove NAs
sims <- sims[!is.na(sims$Origin), ]

# Perpare deisgn
design <- as_tibble(
  expand.grid(
      StepNumber       = seq(1, 2000, by = 1)
    , Origin           = c("Moremi", "Hwange")
    , stringsAsFactors = F
  )
)

# Nest by origin
design <- design %>% nest(Steps = -Origin)
print(design)

# Loop through the design and create desired heatmaps
design$Heatmaps <- lapply(1:nrow(design), function(x){

  # Extract the steps through which we should loop
  steps <- as.vector(as.matrix(design$Steps[[x]]))

  # Subset to correct simulations
  sims_sub <- sims[sims$Origin == design$Origin[x], ]

  # Define extent
  ext <- extent(c(
      min(sims_sub$x)
    , max(sims_sub$x)
    , min(sims_sub$y)
    , max(sims_sub$y)
  )) + c(-1, +1, -1, +1) * metersToDegrees(1000)

  # Prepare reference raster
  r <- raster(ext, res = metersToDegrees(1000))

  # Loop through different number of steps and rasterize trajectories
  heatmaps <- lapply(steps, function(y){
    cat("Create Tracks\n")
    sub <- sims_sub[sims_sub$StepNumber <= y, ]
    sub_track <- sims2tracks(sub)
    cat("Create Heatmap\n")
    heatmap <- rasterizeSpatstat(sub_track, r, mc.cores = 1)
    names(heatmap) <- y
    heatmap <- writeRaster(heatmap, tempfile())
    return(heatmap)
  })
  heatmaps <- stack(heatmaps)

  # Return the heatmaps
  return(heatmaps)

})

# # Identify tracks leaving from moremi or hwange
# sims <- subset(sims, Origin == "Moremi")
#
# # Create trajectories
# sims_tracks <- sims2tracks(sims)
#
# design$Heatmap <- lapply(1:nrow(design), function(x){
#   cat("Create Tracks\n")
#   sub <- sims[sims$StepNumber <= design$StepNumber[x], ]
#   sub_track <- sims2tracks(sub)
#   cat("Create Heatmap\n")
#   heatmap <- rasterizeSpatstat(sub_track, r, mc.cores = 1)
#   heatmap <- writeRaster(heatmap, tempfile())
#   return(heatmap)
# })
#
# # Put all heatmaps together
# heatmaps <- stack(design$Heatmap)

# Store them
writeRaster(heatmaps, "heatmaps_animation.grd")
heatmaps <- stack("03_Data/03_Results/heatmaps_animation.grd")

################################################################################
#### Animate
################################################################################
# Define extent for which we want to plot
extent <- extent(22, 27, -20.5, -17.5) %>% as(., "SpatialPolygons")
crs(extent) <- CRS("+init=epsg:4326")

# # Get extent of kaza
# kaza <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
# africa <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
#
# # Prepare inset map
# test <- tm_shape(kaza) +
#     tm_borders(col = "black") +
#   # tm_shape(africa) +
#   #   tm_borders(col = "gray80") +
#     # tm_text("COUNTRY") +
#   tm_shape(subset(prot, Name == "Moremi")) +
#     tm_polygons(
#         col = "red"
#       , lwd = 0
#     ) +
#   tm_shape(extent) +
#     tm_borders(
#         col = "red"
#       , lty = 2
#     ) +
#   tm_layout(
#       bg.color      = "transparent"
#     , inner.margins = c(0, 0, 0, 0)
#     , outer.margins = c(0, 0, 0, 0)
#     , frame.lwd     = 0
#   )

# Prepare plots for animation
anim <- lapply(1:nlayers(heatmaps), function(n){
  tm_shape(extent) +
      tm_borders() +
    tm_shape(sqrt(heatmaps[[n]])) +
      tm_raster(
          palette            = "magma"
        , style              = "cont"
        , title              = expression("Traversal Frequency"^"0.5")
      ) +
    tm_shape(subset(prot, Name == "Moremi")) +
      tm_borders(
          col = "white"
        , lwd = 1
      ) +
      tm_text(
         "Name"
        , col      = "white"
        , fontface = 3
        , alpha    = 0.8
      ) +
    tm_layout(
        title              = paste0("Steps: ", sprintf("%04d", n))
      , title.color        = "black"
      , title.size         = 0.8
      , title.bg.color     = "white"
      , title.position     = c("center", "top")
      # , legend.position    = c("left", "top")
      # , legend.text.color  = "white"
      # , legend.title.color = "white"
      , legend.show      = F
    ) +
    tm_graticules(
        n.y                 = 5
      , n.x                 = 5
      , labels.inside.frame = FALSE
      , lines               = FALSE
      , ticks               = TRUE
    ) +
    tm_scale_bar(
          position  = c("right", "bottom")
        , text.size = 0.5
        , text.col  = "white"
        , width     = 0.125
    ) +
    tm_compass(
        color.dark  = "white"
      , color.light = "white"
      , text.color  = "white"
      , position    = c("left", "bottom")
  )
})

# Also prepare some leading and trailing plots
leading <- tm_shape(extent) +
      tm_borders() +
    tm_shape(sqrt(heatmaps[[1]])) +
      tm_raster(
          palette            = "black"
        , style              = "cont"
        , title              = expression("Traversal Frequency"^"0.5")
      ) +
    tm_shape(subset(prot, Name == "Moremi")) +
      tm_borders(
          col = "white"
        , lwd = 1
      ) +
      tm_text(
         "Name"
        , col      = "white"
        , fontface = 3
        , alpha    = 0.8
      ) +
    tm_layout(
        title              = paste0("Steps: ", sprintf("%04d", 0))
      , title.color        = "black"
      , title.size         = 0.8
      , title.bg.color     = "white"
      , title.position     = c("center", "top")
      # , legend.position    = c("left", "top")
      # , legend.text.color  = "white"
      # , legend.title.color = "white"
      , legend.show      = F
    ) +
    tm_graticules(
        n.y                 = 5
      , n.x                 = 5
      , labels.inside.frame = FALSE
      , lines               = FALSE
      , ticks               = TRUE
    ) +
    tm_scale_bar(
          position  = c("right", "bottom")
        , text.size = 0.5
        , text.col  = "white"
        , width     = 0.125
    ) +
    tm_compass(
        color.dark  = "white"
      , color.light = "white"
      , text.color  = "white"
      , position    = c("left", "bottom")
  )

# Repeat
leading <- lapply(1:100, function(x){
  leading
})

# Repeat last frame
trailing <- lapply(1:200, function(x){
  anim[[length(anim)]]
})

# Combine all
combined <- c(leading, anim, trailing)

# Prepare the animation
ani.options(interval = 1 / 100, ani.width = 1000, ani.height = 700)
saveVideo({
  for (i in 1:length(combined)){
    suppressWarnings(
      print(combined[[i]])
    )
  }
}, video.name = "Dispersal.mp4")

# # Using tmap
# tmap_animation(
#     combined
#   , filename      = "Dispersal.mp4"
#   , width         = 1200
#   , height        = 600
#   , fps           = 50
# )

################################################################################
#### Creating Time-Kernel
################################################################################
library(spatialEco)
library(adehabitatHR)
water <- readOGR("03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters_GEOFABRIK.shp")
heatmaps <- stack("03_Data/03_Results/heatmaps_animation.grd")
prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
prot <- subset(prot, Desig == "National Park")
# Split heatmaps into packages
test <- splitStack(heatmaps, n = 20)
hrs <- lapply(test, function(z){
  summed <- sum(z)
  pts <- rasterToPoints(summed, spatial = T, fun = function(x){x > 0})
  pts <- as(pts, "SpatialPoints")
  ud <- kernelUD(pts)
  hr <- getverticeshr(ud, percent = 95)
  return(hr)
})
hrs <- do.call(rbind, hrs)
plot(hrs, col = rev(viridis(20)), lwd = 0.1)
plot(water, add = T, col = "white", border = NA)
plot(prot, add = T, lwd = 0.2)
plot(subset(prot, Name == "Moremi")
  , border = "cornflowerblue"
  , col    = colTrans("cornflowerblue", 60)
  , add    = T
  , lwd    = 2
)
text(subset(prot, Name == "Moremi"), "Name", col = "cornflowerblue", font = 3)





plot(hello, col = c(rev(viridis(20)), "black"))

plot(sqrt(heatmaps[[2000]]), col = viridis(100))
plot(hrs[20, ], add = T, border = "red")

test <- as.data.frame(heatmaps, xy = T)

# Identify index when first inundated
test$First <- apply(as.matrix(test[, 3:ncol(test)]), 1, function(x){
  suppressWarnings(
    result <- min(which(x > 0))
  )
  return(result)
})

# Convert to raster
hello <- rasterFromXYZ(xyz = test[, c(1, 2, ncol(test))])
hello <- reclassify(hello, rcl = c(2001, Inf, 2001))
plot(hello, col = c(rev(viridis(20)), "black"))
