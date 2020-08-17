############################################################
#### Algoritm to Model the Flood in the Okavango Delta
############################################################
# Description: Here, we try to parametrize a model that allows us to predict the
# spatial extent of the flood at a given point in time. This is heavily inspired
# by the following post: https://stats.stackexchange.com/questions/244042/trend-
# in-irregular-time-series-data.

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)            # To wrangle code
library(mgcv)                 # For spatial modelling
library(viridis)              # To access nice colors
library(parallel)             # To use multiple cores
library(lubridate)            # To handle dates
library(tictoc)               # To keep track of time
library(pkgcond)              # To suppress warnings
library(raster)               # To manipulate raster data
library(ROCR)                 # To identify threshold
library(rgdal)                # To store shapefiles
library(animation)            # Package to animate
library(spatstat)             # To calculate distances on raster
library(maptools)             # To calculate distances on raster

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load custom functions
source("Functions.r")

# Identify all classified floodmaps
files <- dir(
    path        = "03_Data/02_CleanData/00_Floodmaps/01_Original"
  , pattern     = ".tif$"
  , full.names  = T
)

# Load them
flood <- stack(files)

# Reclassify pixels to flooded (1) not flooded (0) and clouded (NA)
rcl <- data.frame(old = c(0, 127, 255), new = c(1, NA, 0))
flood <- mcreclassify(flood, rcl)

# Calculate how often each pixel was flooded
no_inundated <- calc(flood, function(x){sum(x, na.rm = T)})
names(no_inundated) <- "NoInundated"

# We will use this raster frequently. Thus, let's save it to file
writeRaster(no_inundated
  , filename  = "03_Data/02_CleanData/03_LandscapeFeatures_NoInundated.tif"
  , overwrite = TRUE
)

# Function to create dataframe from stack of floodmaps and covariate
# "no_inundated"
createData <- function(x, y, n = NA){
  if (!is.na(n)){x <- x[[sample(nlayers(x), n)]]}
  dat <- x %>%
    stack(., y) %>%
    as.data.frame(xy = T) %>%
    mutate(PixelID = 1:nrow(.)) %>%
    gather(key = Date, value = Inundated, 3:(ncol(.) - 2)) %>%
    mutate(Date = substr(Date, start = 2, stop = 11)) %>%
    mutate(Date = as.Date(Date, format = "%Y.%m.%d")) %>%
    subset(NoInundated > 0) %>%
    mutate(DayOfYear = yday(Date)) %>%
    mutate(Timestamp = as.numeric(Date))
  return(dat)
}

# Apply the function to the stack of floodmaps and the covariate "no_inundated"
dat <- createData(flood, no_inundated)

############################################################
#### Fit Spatio-Temporal Flood-Model
############################################################
# Subsample and split data into training and testing set
sub <- dat[sample(nrow(dat), 5e6),]
sub <- splitDat(sub, ratio = 0.75)

# Fit GAM model
mod <- bam(Inundated ~
  + s(NoInundated, bs = "tp")
  + s(Timestamp, bs = "tp")
  + s(DayOfYear, bs = "cc")
  + s(x, y, bs = "ds")
  + ti(NoInundated, Timestamp, bs = c("tp", "tp"))
  + ti(NoInundated, DayOfYear, bs = c("tp", "cc"))
  + ti(NoInundated, x, y, d = c(1, 2), bs = c("tp", "ds"))
  + ti(Timestamp, DayOfYear, bs = c("tp", "cc"))
  + ti(Timestamp, x, y, d = c(1, 2), bs = c("tp", "ds"))
  + ti(DayOfYear, x, y, d = c(1, 2), bs = c("cc", "ds"))
  , family    = binomial(link = "logit")
  , data      = sub$Training
  , method    = "fREML"
  , nthreads  = 12
  , discrete  = TRUE
  , gc.level  = 2
)

# # Write the model results to file
# write_rds(mod, "FloodModel.rds")
#
# Look at the model results
summary(mod)
#
# # Visualize them
# plot(mod, pages = 2, shade = T, scale = 0)
#
# # Check diagnostic plots
# gam.check(mod)
#
# # Check for temporal dependence in residuals
# par(mfrow = c(2, 1))
# acf(residuals(mod), lag.max = 200, main = "ACF")
# pacf(residuals(mod), lag.max = 200, main = "pACF")

############################################################
#### Identify Classification Threshold using ROCR
############################################################
# Function to calculate confusion matrix, sensitivity, and specificity. Note
# that x represents predicted values, y represents true values
confMat <- function(x, y, threshold = 0.5){

  # Confusion table
  conf <- table(x, y >= threshold)

  # Performance metrics
  specificity <- conf[1, 1] / sum(conf[1, 1], conf[1, 2])
  sensitivity <- conf[2, 2] / sum(conf[2, 1], conf[2, 2])
  accuracy <- sum(conf[1, 1], conf[2, 2]) / sum(conf[1, ], conf[2, ])

  # Results that are returned
  results <- list(
      Confusion   = conf
    , Specificity = specificity
    , Sensitivity = sensitivity
    , Accuracy    = accuracy
  )
  return(results)
}

# Use the floodmodel to predict the flood extent on the training dataset (to
# find the threshold)
predictTrain <- predict(mod, type = "response")

# Create prediction object. Note that we can only use data that is not NA, as
# the bam removed any NA values automatically
predictROCR <- prediction(
    predictions = predictTrain
  , labels      = sub$Training$Inundated[!is.na(sub$Training$Inundated)]
)

# Create performance object
performROCR <- performance(predictROCR, "tpr", "fpr")

# Plot the results for different thresholds
plot(performROCR
  , colorize          = T
  , print.cutoffs.at  = seq(0, 1, by = 0.1)
  , text.adj          = c(-0.2, 1.7)
)

# It looks like a threshold of 0.5 works perfectly fine. Let's use it to predict
# on the test data
predictTest <- predict(mod, sub$Testing, type = "response")

# Check the confusion matrix
confMat(sub$Testing$Inundated, predictTest, threshold = 0.5)

############################################################
#### Predict Flood
############################################################
# Reload model results
mod <- read_rds("03_Data/03_Results/99_FloodModel.rds")

# Load the covariate layer "no_inundated"
no_inundated <- "03_Data/02_CleanData/03_LandscapeFeatures_NoInundated.tif" %>%
  raster()

# Write a function that predicts the flood extent for a desired date
floodPred <- function(mod, covars, dates = "Desired Date", threshold = 0.5){

  # Prepare new dataframe based on the covariate layer
  dat_new <- covars %>%
    as.data.frame(., xy = TRUE) %>%
    mutate(CellNo = 1:nrow(.)) %>%
    subset(., NoInundated > 0) %>%
    mutate(Date = as.Date(dates)) %>%
    mutate(DayOfYear = yday(Date)) %>%
    mutate(Timestamp = as.integer(Date))

  # Predict likelihood of being inundated and apply threshold
  pred <- predict(mod, dat_new, type = "response")
  pred <- if_else(pred > threshold, 1, 0)

  # Put values into a new raster. Remember that we only predict for cells that
  # have been inundated before
  map <- no_inundated
  values(map)[values(map) > 0] <- pred

  # Write raster to temporary file
  map <- writeRaster(map, rasterTmpFile())

  # Put date of the map as its name
  names(map) <- dates

  # Return the prediction
  return(map)
}

# Create vector of desired dates
dates <- seq(as.Date("2015-01-01"), as.Date("2015-12-31"), by = "1 days")

# Run the prediction to these dates and stack the results
predictions <- mclapply(dates, mc.cores = detectCores() - 1, function(x){
  floodPred(mod = mod, covars = no_inundated, dates = x)
}) %>% stack(.)

# Visualize
plot(predictions[[1]], col = viridis(50))

############################################################
#### Predict Flood Repeatedly
############################################################
# Parameterize floodmodel using random subsets of the total data, then
# predicting flood extent. Let's load the necessary data.

# Function to refit model and predict flood repeatedly
floodPred2 <- function(
    mod
  , covars
  , dat
  , n
  , dates     = "Desired Dates"
  , threshold = 0.5
  , size      = 1e5
  ){

  # Prepare empty list into which the results go
  results <- list()

  # Refit model and predict flood
  for (i in 1:n){
    sub <- dat[sample(nrow(dat), 1e5),]
    mod <- stats::update(mod, data = sub)

    # Apply prediction to these dates
    predictions <- mclapply(dates, mc.cores = detectCores() - 1, function(x){
      floodPred(mod = mod, covars = no_inundated, dates = x)
    }) %>% stack(.)

    # Put results in stack and sum them
    results[[i]] <- predictions

    # Print status
    cat(i, "of", n, "done...\n")
  }

  # Calculate sum over all maps of the same date
  new <- list()
  for (i in 1:nlayers(results[[1]])){
     stacked <- stack(lapply(results, function(x){x[[i]]}))
     stacked <- sum(stacked)
     names(stacked) <- names(results[[1]][[i]])
     new[[i]] <- stacked
  }

  # Stack
  results <- stack(new)

  # Return the results
  return(results)
}

# Select some dates
dates <- seq(as.Date("2015-05-01"), as.Date("2015-05-30"), by = "1 days")

# Run the function
pred <- floodPred2(
    mod       = mod
  , covars    = no_inundated
  , dat       = dat
  , dates     = dates
  , threshold = 0.5
  , n         = 10
  , size      = 1e5
)

# Plot results
plot(pred[[1]])

############################################################
#### Comparison of Predicted and True Flood-Extent
############################################################
# Randomly select some dates for which we have proper floodmaps
dates <- flood %>%
  names(.) %>%
  length(.) %>%
  sample(., 5) %>%
  flood[[.]] %>%
  names(.) %>%
  substr(start = 2, stop = 11) %>%
  as.Date(., format = "%Y.%m.%d")

# Predict flood for these dates
predictions <- mclapply(dates, function(x){
  floodPred(no_inundated, x)
}, mc.cores = n - 1) %>% stack(.)
names(predictions) <- dates

# Plot the results
index <- 1:4
par(mfrow = c(length(index), 2), mar = c(2, 1, 1, 1))
for(i in 1:length(index)){

  # Plot the predicted map
  plot(
      predictions[[index[i]]]
    , main    = paste("Prediction", dates[i])
    , legend  = F
  )

  # Plot the remote sensed (true) map
  plot(
      flood[[names(flood)[names(flood) == names(predictions)[index[i]]]]]
    , main    = paste("True", dates[i])
    , legend  = F
  )
}

############################################################
#### Floodmaps for Dispersal Simulations
############################################################
# Reload flood model
mod <- read_rds("03_Data/03_Results/99_FloodModel.rds")

# Create vector of desired dates
dates <- seq(as.Date("2015-01-01"), as.Date("2015-01-02"), by = "1 days")

# Apply prediction to these dates
predictions <- mclapply(dates, mc.cores = detectCores() - 1, function(x){
  floodPred(mod = mod, covars = no_inundated, dates = x)
}) %>% stack(.)

# Rename the predictions
flood <- predictions

# Load and crop reference raster
r250 <- "03_Data/02_CleanData/00_General_Raster250.tif" %>%
  raster() %>%
  crop(., flood)

# Resample all floodmaps to the reference raster
resampled <- list()
for (i in 1:nlayers(flood)){
  map <- resample(flood[[i]], r250, "ngb")
  map <- writeRaster(map, tempfile())
  resampled[[i]] <- map
  cat(i, "of", length(flood), "done...\n")
}

# Put them back together
flood <- stack(resampled)

# Load globeland and merit layers
globe <- raster("03_Data/02_CleanData/01_LandCoverClasses30_Globeland.tif")
merit <- raster("03_Data/02_CleanData/03_LandscapeFeatures_Rivers_Merit.tif")

# Define the extent for which we will do the simulations (with some buffer)
dat <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv()
ext <- extent(min(dat$x), max(dat$x), min(dat$y), max(dat$y)) +
  1 / 111 * c(-140, +50, -50, +50)
ext <- as(as(ext, "SpatialPolygons"), "SpatialPolygonsDataFrame")
crs(ext) <- CRS("+init=epsg:4326")

# Save it for later
writeOGR(ext
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_DispersalSimulationExtent"
  , driver    = "ESRI Shapefile"
  , overwrite = T
)

# Check the extent
plot(globe)
plot(ext, add = T)
plot(extent(flood), add = T, col = "red")

# Crop the layers to this extent
globe <- crop(globe, ext)
merit <- crop(merit, ext)
flood <- crop(flood, ext)

# We need to remove the globeland waterbodies (Code 40) for the extent of the
# dynamic floodmaps. Let's get a polygon for the extent for which we have
# dynamic floodmaps first
p <- as(extent(predictions[[1]]), "SpatialPolygons")

# Replace the values below the polygon to 0
globe[cellsFromExtent(globe, p)] <- 0

# Before we add the dynamic floodmaps, let's merge the globeland and merit data
globe <- max(globe, merit)

# We need to extend the floodmaps to match the extent of the globeland layer
flood_new <- list()
for (i in 1:nlayers(flood)){
  flood_ext <- extend(flood[[i]], globe, value = NA)
  flood_ext <- writeRaster(flood_ext, tempfile())
  flood_new[[i]] <- flood_ext
  cat(i, "of", nlayers(flood), "done...\n")
}

# Put all back into a stack
flood <- stack(flood_new)

# Finally, we can transfer the water values from the floodmaps to the globeland
# land cover classes
globe <- mask(globe, flood, maskvalue = 1, updatevalue = 1)
names(globe) <- names(flood)

# Store the resulting stack
writeRaster(globe
  , filename  = "03_Data/02_CleanData/03_LandscapeFeatures_PredictedFlood"
  , format    = "raster"
  , overwrite = TRUE
  , options   = c("INTERLEAVE = BAND", "COMPRESS = LZW")
)

############################################################
#### Distance to Water for Dispersal Simulations
############################################################
# Split the stack to make use of mclapply
globe <- splitStack(globe, n = detectCores() - 1)

# Calculate distance to water
distances <- mclapply(globe, mc.cores = detectCores() - 1, function(x){
  distanceTo(x, value = 1)
}) %>% stack() %>% setNames(names(globe))

# Write the resulting stack to file
writeRaster(distances
  , filename  = "03_Data/02_CleanData/03_LandscapeFeatures_PredictedDistanceToWater"
  , format    = "raster"
  , overwrite = TRUE
  , options   = c("INTERLEAVE = BAND", "COMPRESS = LZW")
)

############################################################
#### Animation
############################################################
# Create an animation of the flood
ani.options(interval = .025, ani.width = 1920, ani.height = 1080)
saveVideo({
  for (i in 1:nlayers(predictions)){
    plot(predictions[[i]]
      , col     = c("white", "blue")
      , legend  = FALSE
      , axes    = FALSE
      , box     = FALSE
    )
    text(23.8, -20.5, dates[i], col = "black", cex = 4)
  }
}, video.name = "99_Fluctuations.mp4")

# End cluster
endCluster()
