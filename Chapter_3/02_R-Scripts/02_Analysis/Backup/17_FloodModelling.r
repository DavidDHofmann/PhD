################################################################################
#### Algoritm to Model the Flood in the Okavango Delta
################################################################################
# Description: Here, we try to parametrize a model that allows us to predict the
# spatial extent of the flood at a given point in time. This is heavily inspired
# by the following post: https://stats.stackexchange.com/questions/244042/trend-
# in-irregular-time-series-data. Also check out this article from medium
# https://medium.com/analytics-vidhya/a-guide-to-machine-learning-in-r-for-beginners-part-5-4c00f2366b90

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)            # To wrangle code
library(mgcv)                 # For spatial modelling
library(viridis)              # To access nice colors
library(pbmcapply)            # To use multiple cores and show progress bars
library(lubridate)            # To handle dates
library(tictoc)               # To keep track of time
library(pkgcond)              # To suppress warnings
library(raster)               # To manipulate raster data
library(ROCR)                 # To identify threshold
library(rgdal)                # To store shapefiles
library(animation)            # Package to animate
library(spatstat)             # To calculate distances on raster
library(maptools)             # To calculate distances on raster
library(davidoff)             # Custom functions
library(mgcViz)               # For nice plots

# Set seed for reproducability
set.seed(12345)

# Set the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_1")

# Identify all classified floodmaps
files <- dir(
    path        = "03_Data/02_CleanData/00_Floodmaps/01_Original"
  , pattern     = ".tif$"
  , full.names  = T
)

# Identify dates
filedates <- as.Date(substr(basename(files), start = 1, stop = 10))

# Load them
flood <- stack(files)

# Reclassify pixels to flooded (1) not flooded (0) and clouded (NA) -> In case
# you don't want to use terra::classify, use davidoff::mcreclassify
rcl <- data.frame(old = c(0, 127, 255), new = c(1, NA, 0))
flood <- mcreclassify(flood, rcl)

# Calculate how often each pixel was flooded
no_inundated <- sum(flood, na.rm = T)
names(no_inundated) <- "NoInundated"

# We will use this raster frequently. Thus, let's save it to file
writeRaster(no_inundated
  , filename  = "03_Data/02_CleanData/03_LandscapeFeatures_NoInundated.tif"
  , overwrite = TRUE
)

# Function to translate stack of floodmaps into a dataframe and add a column
# indicating the covariate "no_inundated"
createData <- function(x, y, n = NA){
  if (!is.na(n)){x <- x[[sample(nlyr(x), n)]]}
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

# Clear cache
gc()

################################################################################
#### Fit Spatio-Temporal Flood-Model
################################################################################
# Subsample data and split into training and testing set
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

# Look at the model results
summary(mod)

# Write the model results to file
write_rds(mod, "03_Data/03_Results/99_FloodModel.rds")

# Remove data
rm(dat)
gc()

# Visualize them
plot(mod, pages = 2, shade = T, scale = 0)

# Check diagnostic plots
gam.check(mod)
#
# # Check for temporal dependence in residuals
# par(mfrow = c(2, 1))
# acf(residuals(mod), lag.max = 200, main = "ACF")
# pacf(residuals(mod), lag.max = 200, main = "pACF")

# Some nice visuals
b <- getViz(mod)

# On individual covariates
o <- plot( sm(b, 2) )
o + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

# Summary plot
print(plot(b, select = 1:3), pages = 1)

################################################################################
#### Identify Classification Threshold using ROCR
################################################################################
# Reload model results
mod <- read_rds("03_Data/03_Results/99_FloodModel.rds")

# To predict a floodmap, we'll need to specify a classification threshold. We
# can balance false positives against false negatives to find a meaningful
# threshold. Function to calculate confusion matrix, sensitivity, and
# specificity. Note that x represents predicted values, y represents true values
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
confMat(sub$Testing$Inundated, predictTest, threshold = 0.50)

################################################################################
#### Predict Flood
################################################################################
# Reload model results
mod <- read_rds("03_Data/03_Results/99_FloodModel.rds")

# Load the covariate layer "no_inundated"
no_inundated <- raster("03_Data/02_CleanData/03_LandscapeFeatures_NoInundated.tif")
names(no_inundated) <- "NoInundated"

# Write a function that predicts the flood extent for desired dates
floodPred <- function(mod, covars, dates, threshold = 0.5){

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

  # Clear cache
  gc()

  # Return the prediction
  return(map)
}

# Create vector of desired dates. We will predict the flood for the same range
# of dates for which we have input data
dates <- seq(
    from  = as.Date("2019-01-01")
  , to    = as.Date("2019-12-31")
  , by    = "7 days"
)

# Run the prediction to these dates and stack the results
predictions <- pbmclapply(
    X                   = dates
  , ignore.interactive  = T
  , mc.cores            = detectCores() - 1
  , FUN                 = function(x){
    floodPred(mod = mod, covars = no_inundated, dates = x)
  }) %>% stack(.)

# Visualize
plot(predictions[[300]], col = viridis(50))

################################################################################
#### Comparison of Predicted and True Flood-Extent
################################################################################
# # First we compare the cumulated extent
# dat1 <- data.frame(
#     Date    = as.Date(basename(files))
#   , Extent  = as.vector(cellStats(flood, "sum"))
#   , Type    = "Observed"
# )
# dat2 <- data.frame(
#     Date    = dates
#   , Extent  = as.vector(cellStats(predictions, "sum"))
#   , Type    = "Predicted"
# )
# dat <- rbind(dat1, dat2)
# dat <- subset(dat, Date > "2010-01-01")
#
# # Visualize all
# ggplot(dat, aes(x = Date, y = Extent, col = factor(Type))) +
#   geom_line() +
#   geom_rug(data = subset(dat, Type == "Observed"), sides = "b", length = unit(0.01, "npc")) +
#   scale_color_manual(values = c("gray40", "blue"))
# Identify all classified floodmaps

files <- dir(
    path        = "03_Data/02_CleanData/00_Floodmaps/01_Original"
  , pattern     = ".tif$"
  , full.names  = T
)

# Sample five files
files <- sample(files, 5)

# Identify dates
dates <- as.Date(substr(basename(files), start = 1, stop = 10))

# Load them
flood <- stack(files)

# Reclassify pixels to flooded (1) not flooded (0) and clouded (NA) -> In case
# you don't want to use terra::classify, use davidoff::mcreclassify
rcl <- data.frame(old = c(0, 127, 255), new = c(1, NA, 0))
flood <- reclassify(flood, rcl)

# Predict flood for these dates
predictions2 <- lapply(
    X                   = dates
  # , ignore.interactive  = T
  # , mc.cores            = detectCores() - 1
  , FUN                 = function(x){
    floodPred(mod = mod, covars = no_inundated, dates = x)
  }) %>% stack(.)
names(predictions2) <- dates

# Plot the results
png("test3.png", width = 1980, height = 1080)
i <- 2
par(mfrow = c(1, 2), mar = c(2, 1, 1, 1))
plot(
    predictions2[[i]]
  , main    = paste("Prediction", dates[i])
  , legend  = F
  , col = c("white", "blue")
  , box = F
  , axes = F
)
plot(
    flood[[i]]
  , main    = paste("True", dates[i])
  , legend  = F
  , col = c("white", "blue")
  , box = F
  , axes = F
)
dev.off()

################################################################################
#### Floodmaps for Dispersal Simulations
################################################################################
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

################################################################################
#### Distance to Water for Dispersal Simulations
################################################################################
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

################################################################################
#### Animation
################################################################################
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
