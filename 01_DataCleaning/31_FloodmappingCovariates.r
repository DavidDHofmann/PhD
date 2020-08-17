############################################################
#### Prepare all Covariates Needed to Model the Flood
############################################################
# Description: Here, we prepare covariates that will hopefully allow us to
# predict the spatial extent of the flood at a given point in time. To do so, we
# will create an empty raster of the study area and to each pixel collect
# relevant covariates, such as elevation, xy-coordinates etc.

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)       # To handle raster data
library(tidyverse)    # To wrangle code
library(lubridate)    # To handle dates nicer
library(parallel)     # To use multiple cores
library(elevatr)      # To access elevation data
library(viridis)      # To access nice colors
library(pkgcond)      # To suppress warnings

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load custom functions
source("Functions.r")

# Make use of multicore abilities
beginCluster()

############################################################
#### Prepare Empty Raster
############################################################
# We will load one of the cassified floodmaps as reference raster
r <- "03_Data/02_CleanData/00_Floodmaps/01_Original/2000.04.30.tif" %>%
  raster()
values(r) <- 1:ncell(r)

# Visualize
plot(r, col = viridis(100))

# Add to stack
stacked <- stack(r)
names(stacked) <- "id"

############################################################
#### Covariates: Elevation, Aspect, Slope
############################################################
# Load DEM data
elev <- get_elev_raster(r, z = 8)
elev <- crop(elev, r)

# Reproject the DEM to fit the reference raster
elev <- resample(elev, r, method = "bilinear")

# Also calculate aspect and slope
terrain <- terrain(elev, c("aspect", "slope"), unit = "degrees")

# We can't calculate aspect and slope for the cells at the very edge of the
# raster
terrain[[1]][1, ]
terrain[[1]][, 1]
terrain[[1]][nrow(terrain[[1]]), ]
terrain[[1]][, ncol(terrain[[1]])]

# However, I don't want to trim them and will simply copy edges that are not NaN
terrain[[1]][1, ] <- terrain[[1]][1 + 1, ]
terrain[[1]][nrow(terrain[[1]]), ] <- terrain[[1]][nrow(terrain[[1]]) - 1, ]
terrain[[1]][, 1] <- terrain[[1]][, 1 + 1]
terrain[[1]][, ncol(terrain[[1]])] <- terrain[[1]][, ncol(terrain[[1]]) - 1]

terrain[[2]][1, ] <- terrain[[2]][1 + 1, ]
terrain[[2]][nrow(terrain[[2]]), ] <- terrain[[2]][nrow(terrain[[2]]) - 1, ]
terrain[[2]][, 1] <- terrain[[2]][, 1 + 1]
terrain[[2]][, ncol(terrain[[2]])] <- terrain[[2]][, ncol(terrain[[2]]) - 1]

# Visualize
plot(elev, col = viridis(100))
plot(sqrt(sqrt(terrain$slope)), col = viridis(100))
plot(terrain$aspect, col = viridis(100))

# Give it a nice name
names(elev) <- "elev"

# Add them to the stack
stacked <- stack(stacked, elev, terrain)

############################################################
#### Covariate: Distributary Region
############################################################
# Load polygons of different distributaries
distr <- "03_Data/01_RawData/DAVID/Distributaries.shp" %>%
  shapefile()

# Rasterize the shapefile
distr_r <- rasterize(distr, stacked[[1]])

# Replace NAs with 0s
values(distr_r)[is.na(values(distr_r))] <- 0

# Visualize
plot(distr_r, col = viridis(15))

# Give it a nice name
names(distr_r) <- "distributaries"

# Add covariate to stack
stacked <- stack(stacked, distr_r)

############################################################
#### Covariate: Inundation History
############################################################
# Load all classified floodmaps
flood <- dir(
    path        = "03_Data/02_CleanData/00_Floodmaps/01_Original"
  , pattern     = ".tif$"
  , full.names  = T
) %>% stack()

# Reclassify pixels to flooded (1) not flooded (0) and clouded (NA)
rcl <- data.frame(old = c(0, 127, 255), new = c(1, NA, 0))
flood <- mcreclassify(flood, rcl)

# How often each pixel was inundated
sum <- calc(flood, function(x){sum(x, na.rm = T)})

# Give the layers nicer names
names(sum) <- "no_inundated"
names(flood) <- flood %>%
  names() %>%
  substr(start = 2, stop = 11) %>%
  paste0("inundation_", .)

# We can assign dates to the stack
dates <- names(flood) %>%
  substr(start = 12, stop = 25) %>%
  as.Date(format = "%Y.%m.%d")
flood <- setZ(flood, dates)

# Visualize
plot(sum, col = viridis(100))
plot(flood[[1:4]], col = viridis(100))

# Add summed layer to stack (we will store the inundation history seperately)
stacked <- stack(stacked, sum)

# Write the inundation history to file
writeRaster(
    flood
  , filename  = "03_Data/02_CleanData/03_LandscapeFeatures_InundationHistory"
  , format    = "raster"
  , overwrite = TRUE
  , options   = c("INTERLEAVE = BAND", "COMPRESS = LZW")
)

############################################################
#### Covariate: Pixel Rank
############################################################
# We saw that most peaks occur in the middle of the year (see previous script).
# We can therefore use the beginning of the year to identify pixel ranks. Let's
# prepare a tibble for this containing raster stacks for each year
ranks <- tibble(Year = unique(year(getZ(flood))))

# Now we create stacks containing all floodmaps relating to that year
ranks <- mutate(ranks, Flood = map(Year, function(x){
  index <- which(year(getZ(flood)) == x)
  maps <- flood[[index]]
  return(maps)
}))

# Rather than raster images we need to convert the stacks to dataframes
ranks <- mutate(ranks, Data = map(Flood, function(x){
  as.data.frame(x, xy = T)
}))

# Calculating the pixel ranks is computationally expensive. To reduce the
# workload we will thus only consider pixels that have been inundated at least
# once
index <- values(sum) > 0
sum(index)
ranks <- mutate(ranks, Data = map(Data, function(x){
  x[index, ]
}))

# Keep track of the number of maps we have for each year
ranks <- mutate(ranks, NoMaps = map(Flood, nlayers) %>% do.call(rbind, .))

# Finally we can calculate pixel ranks (on multiple cores)
ranks <- mutate(ranks, PixelRank = mclapply(Data, function(x){
  apply(as.matrix(x[, 3:ncol(x)]), 1, function(y){
    suppress_warnings(min(which(y == 1)))
  })
}, mc.cores = (detectCores()-1)))

# Replace Inf with NA
ranks <- mutate(ranks, PixelRank = map(PixelRank, function(x){
  r <- as.data.frame(x)
  r[r == Inf] <- NA
  return(r)
}))

# Normalize rank scores
ranks <- mutate(ranks, PixelRankNormal = map(PixelRank, function(x){
  (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
}))

# Finally we can access the normalized ranks
ranks_normal <- do.call(cbind, ranks$PixelRankNormal)

# Calculate averages for each pixel
ranks_mean <- rowMeans(ranks_normal, na.rm = T) %>% as.data.frame()
df <- data.frame(ID = 1:ncell(r), Rank = NA)
df$Rank[index] <- ranks_mean[[1]]

# Create a new raster
ranks <- raster(sum)

# And store the pixel rank in there
ranks <- setValues(ranks, df$Rank)

# Visualize result
plot(ranks, col = viridis(50))

# Assign nice name
names(ranks) <- "rank"

# Add to stack
stacked <- stack(stacked, ranks)

# This was the last static covariate. We can now save the stack to file
writeRaster(
    stacked
  , filename  = "03_Data/02_CleanData/03_LandscapeFeatures_FloodmappingCovariates"
  , format    = "raster"
  , overwrite = TRUE
  , options   = c("INTERLEAVE = BAND", "COMPRESS = LZW")
)

# ############################################################
# #### Covariate: Local Rainfall
# ############################################################
# # Load rainfall data
# files <- dir(path = "CHIRPS", pattern = ".tif$", full.names = T)
# rain <- stack(files)
#
# # Crop the raindata to our extent
# rain <- mccrop(rain, stacked)
#
# # We also need to make sure that the resolutions are the same as in the stack
# factor <- round(res(rain) / res(stacked))[1]
# rain <- disaggregate(rain, fact = factor)
#
# # Make sure that summed values of disaggregated raster equal to the original
# # raster
# rain <- rain / factor ** 2
#
# # Resample the data to make the origin the same to our reference raster
# rain <- resample(rain, stacked, method = "bilinear")
#
# # Give the layers nicer names
# names(rain) <- rain %>%
#   names() %>%
#   substr(start = 13, stop = 20) %>%
#   paste0("rain_", .)
#
# # We can assign dates to the stack
# dates <- names(rain) %>%
#   substr(start = 6, stop = 14) %>%
#   paste0(., ".15") %>%
#   as.Date(format = "%Y.%m.%d")
# rain <- setZ(rain, dates)
#
# # Visualize
# plot(rain[[1:4]], col = viridis(100))
#
# # Save rain data to file
# writeRaster(
#     rain
#   , filename = "CleanData/LocalRain"
#   , format = "raster"
#   , overwrite = TRUE
#   , options = c("INTERLEAVE = BAND", "COMPRESS = LZW")
# )
