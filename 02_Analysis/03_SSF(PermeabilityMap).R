############################################################
#### Predicting a Permeability Surface and Calculating LCPs
############################################################
# Description: In this script we use the result from our model selection to
# predict a permeability surface. On this surface we can then calculate least
# cost paths or least cost corridors

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(ggpubr)
library(raster)
library(data.table)
library(rgeos)
library(tidyverse)
library(glmmTMB)
library(viridis)
library(rasterVis)
library(ggplot2)
library(velox)
library(maptools)
library(spatstat)

# Load custom functions
source("Functions.r")

############################################################
#### Reload Model Results
############################################################
# Read in the model results from the analysis of the 4-hourly fixes
models <- readRDS("03_Data/03_Results/99_ModelSelection.rds")

# Look at the model results again
print(models, n = 100)

# Subset to the best model and identify the coefficients
coeffs <- models$Model[[1]] %>% getCoeffs()

# Plot them
showCoeffs(coeffs[-1, ])

# We want to use the derived model to predict a resistance surface.
# Let's create a simplified dataframe of the model coefficients
pred <- as.data.frame(select(coeffs, Coefficient))

# For easier indexing we assign the covariates as row-names
rownames(pred) <- coeffs$Covariate

# Look at the result
pred

############################################################
#### Load the Covariate Layers
############################################################
# Specify the filenames of the layers we might need for prediction
layers <- c(
      "01_LandCover_Water_Averaged.tif"
    , "01_LandCover_Water_DistanceToWater.tif"
    , "01_LandCover_TreeCover_Averaged.tif"
    , "01_LandCover_NonTreeVegetation_Averaged.tif"
    , "02_LandUseTypes_Protected_PeaceParks(1Class).tif"
    , "04_AnthropogenicFeatures_HumanInfluence(Buffer5000).tif"
    , "04_AnthropogenicFeatures_DistanceToHumans.tif"
    , "04_AnthropogenicFeatures_DistanceToRoads.tif"
  ) %>% paste0("03_Data/02_CleanData/", .) %>%

  # Load each file as rasterlayer into a list
  lapply(., raster) %>%

  # Set some nice layernames
  set_names(c(
      "Water"
    , "DistanceToWater"
    , "Trees"
    , "Shrubs"
    , "Protected"
    , "HumansBuff5000"
    , "DistanceToHumans"
    , "DistanceToRoads"
    )
  )

# We also need to take the sqrt for the "DistanceTo"-layers since we fitted the
# model using the sqrt of these distances
layers[["DistanceToWater"]]   <- sqrt(layers[["DistanceToWater"]])
layers[["DistanceToHumans"]]  <- sqrt(layers[["DistanceToHumans"]])
layers[["DistanceToRoads"]]   <- sqrt(layers[["DistanceToRoads"]])

# Assign layernames again
names <- names(layers)
for (i in 1:length(layers)){
  names(layers[[i]]) <- names[i]
}

############################################################
#### Scale Covariate Layers
############################################################
# Before we scale we create a backup
backup <- layers

# We need to scale the layers using the same scaling parameters that we used to
# scale the covariates during the modelling part
scaling <- read_rds("03_Data/03_Results/99_ScalingSSF.rds")

# Make some nice rownames
rownames(scaling) <- scaling$ColumnName

# Write a function to scale a layer using the table above
scaleLayer <- function(layer, table){
  scale(layer,
      center  = table[names(layer), ]$Center
    , scale   = table[names(layer), ]$Scale
  )
}

# Apply the scaling function to the different layers in the list
layers <- lapply(layers, function(x){
  scaleLayer(x, scaling)
})

# Stack the layers
layers <- stack(layers)

################################################################################
#### NEED TO SCALE THESE VALUES
################################################################################
# Note: We have decided to scale the log_sl_ and cos_ta_ prior to modelling.
# Thus, I'll need to scale them here too.

# We will also need layers for the turning angles and step lengths. Actually,
# they simply inflate the probabilities. For the step length layer we can use
# the average step length (around 2650)
`log(sl_)` <- layers[[1]]
values(`log(sl_)`) <- log(2650)
names(`log(sl_)`) <- "log(sl_)"

# For the turning angles we will use zero
`cos(ta_)` <- layers[[1]]
values(`cos(ta_)`) <- cos(0)
names(`cos(ta_)`) <- "cos(ta_)"

################################################################################
#### NEED TO SCALE THESE VALUES
################################################################################

# Stack all layers into a single rasterstack
layers <- stack(`cos(ta_)`, `log(sl_)`, layers)

# Prepare permeability map using the form "exp(beta0 + beta1 * X1 + ...)"
permeability <- exp(sum(
    pred["cos(ta_)", ]        * layers[["cos.ta_."]]
  , pred["log(sl_)", ]        * layers[["log.sl_."]]
  , pred["Water", ]           * layers[["Water"]]
  , pred["DistanceToWater", ] * layers[["DistanceToWater"]]
  , pred["Shrubs", ]          * layers[["Shrubs"]]
  , pred["HumansBuff5000", ]  * layers[["HumansBuff5000"]]
  , pred["Trees", ]           * layers[["Trees"]]
))

# Remove outliers
upper <- quantile(values(permeability), 0.99, na.rm = TRUE)
lower <- quantile(values(permeability), 0.01, na.rm = TRUE)
values(permeability)[values(permeability) > upper] <- upper
values(permeability)[values(permeability) < lower] <- lower

# # Scale the map values to 0-1
# permeability <- calc(permeability, fun = function(x){
#   (x - min(x)) / (max(x) - min(x))
# })
#
# # Invert the map
# resistance <- calc(permeability, fun = function(x){
#   x * (-1) + max(x)
# })

# Reload the kaza shapefile for plotting
kaza <- shapefile("03_Data/02_CleanData/00_General_KAZA_KAZA")

# Plot the permeability map
png("99_PermeabilityMap.png"
  , width     = 1980
  , height    = 1080
  , bg        = "transparent"
  , pointsize = 30
)
plot(permeability, col = viridis(50))
plot(kaza
  , add = TRUE
  , border = "white"
  , lwd = 5
)
dev.off()

# Store the permeability map to file
writeRaster(permeability
  , "03_Data/03_Results/99_PermeabilityMap.tif"
  , overwrite = TRUE
)

############################################################
#### Compare Permeability (Countries, KAZA, Protected)
############################################################
# Reload permeability map and required shapefiles
perm    <- raster("03_Data/03_Results/99_PermeabilityMap.tif")
kaza    <- shapefile("03_Data/02_CleanData/00_General_KAZA_KAZA")
africa  <- shapefile("03_Data/02_CleanData/00_General_Africa")
prot    <- shapefile("03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(1Class)")

# Aggregate the permeability map to gain efficiency
perm <- aggregate(perm, fact = 4, fun = mean)

# We want to normalize the permeability values
perm <- calc(perm, fun = function(x){
  (x - min(x)) / (max(x) - min(x))
})

# Dissolve the borders
prot <- prot %>% aggregate(., dissolve = TRUE) %>%

  # The geometry is now invalid. We apply a tiny buffer (1m) to make it valid
  gBuffer(., width = 1 / 110 * 0.001) %>%

  # Crop to permeability map
  crop(., perm)

# Polygons for Inside vs. Outside KAZA
count_ins_kaza <- crop(africa, kaza)
count_out_kaza <- erase(crop(africa, extent(perm)), kaza)

# Polygons for Protected vs. Pastoral
count_ins_prot <- crop(africa, prot)
count_out_prot <- erase(crop(africa, extent(perm)), prot)

# Turn the permeability surface to a velox raster
perm <- velox(perm)

# Extract permeability values within vs outside kaza
perm_count_ins_kaza <- perm$extract(count_ins_kaza) %>%
  lapply(., as.data.frame) %>%
  set_names(., count_ins_kaza$COUNTRY) %>%
  rbindlist(., idcol = TRUE) %>%
  mutate(., Group = "Inside KAZA") %>%
  rename(., Permeability = V1, Country = .id)

perm_count_out_kaza <- perm$extract(count_out_kaza) %>%
  lapply(., as.data.frame) %>%
  set_names(., count_out_kaza$COUNTRY) %>%
  rbindlist(., idcol = TRUE) %>%
  mutate(., Group = "Outside KAZA") %>%
  rename(., Permeability = V1, Country = .id)

# Combine the dataframes
count_kaza <- rbind(
    perm_count_ins_kaza
  , perm_count_out_kaza
)

# Remove congo from the dataset
count_kaza <- subset(count_kaza, Country != "Democratic Republic of Congo")

# Extract permeability values within vs outside protected
perm_count_ins_prot <- perm$extract(count_ins_prot) %>%
  lapply(., as.data.frame) %>%
  set_names(., count_ins_prot$COUNTRY) %>%
  rbindlist(., idcol = TRUE) %>%
  mutate(., Group = "Protected") %>%
  rename(., Permeability = V1, Country = .id)

perm_count_out_prot <- perm$extract(count_out_prot) %>%
  lapply(., as.data.frame) %>%
  set_names(., count_out_prot$COUNTRY) %>%
  rbindlist(., idcol = TRUE) %>%
  mutate(., Group = "Pastoral") %>%
  rename(., Permeability = V1, Country = .id)

# Combine the dataframes
count_prot <- rbind(
    perm_count_ins_prot
  , perm_count_out_prot
)

# Remove congo from the dataset
count_prot <- subset(count_prot, Country != "Democratic Republic of Congo")

# Calculate median permeability (rounded to two digits)
count_kaza %>%
  group_by(., Country, Group) %>%
  summarize(.
    , Median  = round(median(Permeability), 2)
    , IQR     = round(IQR(Permeability), 2)
  )
count_kaza %>%
  group_by(., Country) %>%
  summarize(.
    , Median  = round(median(Permeability), 2)
    , IQR     = round(IQR(Permeability), 2)
  )
count_kaza %>%
  group_by(., Group) %>%
  summarize(.
    , Median  = round(median(Permeability), 2)
    , IQR     = round(IQR(Permeability), 2)
  )

count_prot %>%
  group_by(., Country, Group) %>%
  summarize(.
    , Median  = round(median(Permeability), 2)
    , IQR     = round(IQR(Permeability), 2)
  )
count_prot %>%
  group_by(., Country) %>%
  summarize(.
    , Median  = round(median(Permeability), 2)
    , IQR     = round(IQR(Permeability), 2)
  )
count_prot %>%
  group_by(., Group) %>%
  summarize(.
    , Median  = round(median(Permeability), 2)
    , IQR     = round(IQR(Permeability), 2)
  )
round(median(count_prot$Permeability), 2)
round(IQR(count_prot$Permeability), 2)

# Save the data
write_rds(count_kaza, "03_Data/03_Results/99_PermeabilityValues(KAZA).rds")
write_rds(count_prot, "03_Data/03_Results/99_PermeabilityValues(Prot).rds")

# ############################################################
# #### Create Alternative Scenarios (Change Permeability Artificially)
# ############################################################
# # Load shapefile which tells us what to change
# change <- shapefile("99_RemoveResistance")
#
# # Create an alternative covariate layer-stack
# layers <- backup
#
# # We remove the Zambezi river
# layers$Water <- mask(
#     layers$Water
#   , change
#   , inverse = TRUE
#   , updatevalue = 0
# )
#
# # We need to update the distance to water layer
# waterppp <- layers$Water %>%
#
#   # Coerce the inhabited cells to Spatial Points
#   rasterToPoints(., fun = function(x){x == 1}, spatial = TRUE) %>%
#
#   # Transform the points to utm
#   spTransform(., CRS("+init=epsg:32734")) %>%
#
#   # Coerce the object to a ppp object
#   as(., "ppp")
#
# # Load the reference raster to prepare a distance to water
# distance <- raster("00_General_Raster250.tif")
#
# # Now replace the values of the raster with the distances to the nearest water
# # covered cell
# values(distance) <- distance %>%
#
#   # Coerce the layer to a spatial points object
#   as(., "SpatialPoints") %>%
#
#   # Reproject the points to utm
#   spTransform(., CRS("+init=epsg:32734")) %>%
#
#   # Turn the object into a ppp
#   as(., "ppp") %>%
#
#   # Calculate the distance to the nearest water cell
#   nncross(., waterppp) %>%
#
#   # Select only the distance column
#   .[["dist"]]
#
# # Add the layer to the stack
# layers$DistanceToWater <- sqrt(distance)
#
# # Make name proper
# names(layers$DistanceToWater) <- "DistanceToWater"
#
# # Apply the scaling function to the different layers in the list
# layers <- lapply(layers, function(x){
#   scaleLayer(x, scaling)
# })
#
# # Stack the layers
# layers <- stack(layers)
#
# # We will also need layers for the turning angles and step lengths. Actually,
# # they simply inflate the probabilities. For the step length layer we can use
# # the average step length (around 2650)
# `log(sl_)` <- layers[[1]]
# values(`log(sl_)`) <- log(2650)
# names(`log(sl_)`) <- "log(sl_)"
#
# # For the turning angles we will use zero
# `cos(ta_)` <- layers[[1]]
# values(`cos(ta_)`) <- cos(0)
# names(`cos(ta_)`) <- "cos(ta_)"
#
# # Stack all layers into a single rasterstack
# layers <- stack(`cos(ta_)`, `log(sl_)`, layers)
#
# # Prepare permeability map using the form "exp(beta0 + beta1 * X1 + ...)"
# permeability <- exp(sum(
#     pred["cos(ta_)", ]        * layers[["cos.ta_."]]
#   , pred["log(sl_)", ]        * layers[["log.sl_."]]
#   , pred["Water", ]           * layers[["Water"]]
#   , pred["DistanceToWater", ] * layers[["DistanceToWater"]]
#   , pred["Shrubs", ]          * layers[["Shrubs"]]
#   , pred["HumansBuff5000", ]  * layers[["HumansBuff5000"]]
#   , pred["Trees", ]           * layers[["Trees"]]
# ))
#
# # Remove outliers
# upper <- quantile(values(permeability), 0.99, na.rm = TRUE)
# lower <- quantile(values(permeability), 0.01, na.rm = TRUE)
# values(permeability)[values(permeability) > upper] <- upper
# values(permeability)[values(permeability) < lower] <- lower
#
# # Scale the map values to 0-1
# permeability <- calc(permeability, fun = function(x){
#   (x - min(x)) / (max(x) - min(x))
# })
#
# # Invert the map
# resistance <- calc(permeability, fun = function(x){
#   x * (-1) + max(x)
# })
#
# # Store the permeability map to file
# writeRaster(permeability
#   , "99_PermeabilityMap_Alt.tif"
#   , overwrite = TRUE
# )
