################################################################################
#### Step Selection Function - Extraction of Covariates
################################################################################
# Description: In this script we will extract all covariates underlying the
# generated steps

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load packages
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(lubridate)    # To handle dates
library(amt)          # To coerce gps fixes to steps
library(raster)       # To handle spatial data
library(rgdal)        # To handle spatial data
library(pbmcapply)    # For multicore abilities

# Load the generated steps
lines <- readOGR("03_Data/02_CleanData/00_General_Dispersers_Popecol(iSSF).shp")

################################################################################
#### Water
################################################################################
# Load the merged water cover dataset
wat <- brick("03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd")

# Also load the dates
dates <- read_csv("03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.csv")

# Now we can extract the percentage cover of Water along each step
extracted <- extrCov(wat, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(dates)

# Let's look at the result
head(extracted)

# We only want to keep the values from the dates that are closest in time to the
# steps
lines$Water <- sapply(1:nrow(lines), function(x){
  index <- which.min(abs(as.Date(lines$t1_[x]) - dates))[1]
  value <- extracted[x, index]
})

# Remove layers and clear garbage
rm(water)
gc()

################################################################################
#### Trees
################################################################################
# Load the treecover map
trees <- brick("03_Data/02_CleanData/01_LandCover_TreeCover_MODIS.tif")

# Read all layers to memory. This substantially decreases extraction times
trees <- readAll(trees)

# Extract the average tree cover along each line
extracted <- extrCov(trees, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(dates)

# Let's look at the result
head(extracted)

# Keep only values closest in date
lines$Trees <- sapply(1:nrow(lines), function(x){
  index <- which.min(abs(as.Date(lines$t1_[x]) - dates))[1]
  value <- extracted[x, index]
})

# Remove layers and clear garbage
rm(Trees)
gc()

################################################################################
#### Shrubs
################################################################################
# Load the shrubcover map
shrubs <- brick("03_Data/02_CleanData/01_LandCover_NonTreeVegetation_MODIS.tif")

# Read all layers to memory. This substantially decreases extraction times
shrubs <- readAll(shrubs)

# Extract the average shrub cover in each ellipse
extracted <- extrCov(shrubs, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(dates)

# Let's look at the result
head(extracted)

# Keep only values closest in date
lines$Shrubs <- sapply(1:nrow(lines), function(x){
  index <- which.min(abs(as.Date(lines$t1_[x]) - dates))[1]
  value <- extracted[x, index]
})

# Remove layers and clear garbage
rm(Shrubs)
gc()

############################################################
#### Extraction of Covariates: Land Use Types
############################################################
# Load the protected areas raster layer
prot <- "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(1Class).tif" %>%
  raster()

# Extract the percentage coverage
extracted <- extrCov(prot, lines)

# Put the extracted values into the dataframe
lines$Protected <- extracted[, 1]

############################################################
#### Extraction of Covariates: Human Influence
############################################################
# Load the human influence layers
humansBase <- "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Base).tif" %>%
  raster()
humansDens <- "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Density).tif" %>%
  raster()
humansBuff5000  <- "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Buffer5000).tif" %>%
  raster()

# Extract the percentage coverage of cells that are inhabited by humans
extracted <- extrCov(humansBase, lines)

# Put the extracted values into the dataframe
lines$HumansBase <- extracted[, 1]

# Extract the average humans density along each line
extracted <- extrCov(humansDens, lines)

# Put extracted values in dataframe
lines$HumansAverage <- extracted[, 1]

# Extract the 5000m buffered average humans density along each line
extracted <- extrCov(humansBuff5000, lines)

# Put extracted values in dataframe
lines$HumansBuff5000 <- extracted[, 1]

############################################################
#### Extraction of Covariates: Social Features
############################################################
# Load the utilisation distributions (uds)
ud <- "03_Data/02_CleanData/05_SocialFeatures_SocialLandscape(UtilisationDistributions).grd" %>%
  stack()

# Prepare dataframe that shows the date for each ud-layer
ud_info <- data.frame(
    LayerNo = 1:nlayers(ud)
  , Date    = ud %>%
      names() %>%
      substr(start = 2, stop = 11) %>%
      as.POSIXct(format = "%Y.%m.%d", tz = "UTC") %>%
      as.Date()
)

# Look at the table
ud_info

# Load the home ranges (hrs)
hr <- "03_Data/02_CleanData/05_SocialFeatures_SocialLandscape(HomeRanges).grd" %>%
  stack()

# Prepare dataframe that shows date for each hr-layer as well as the pack to
# which the home range refers to
hr_info <- data.frame(
    LayerNo = 1:nlayers(hr)
  , Date    = hr %>%
      names() %>%
      substr(start = 2, stop = 11) %>%
      as.POSIXct(format = "%Y.%m.%d", tz = "UTC") %>%
      as.Date()
  , Pack = hr %>%
      names() %>%
      substr(start = 13, stop = 14)
)

# Look at the table
hr_info

# To extract the hr from the correct layer we need to know each dog's
# pre-dispersal pack. Let's load our gps fixes for this
dat <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv()

# Extract the date at which they start to disperse
cut <- dat %>%
  subset(State == "Disperser") %>%
  group_by(DogName) %>%
  summarize(StartDisp = min(Timestamp)) %>%
  mutate(PrePack = NA)

# Identify predispersal pack
for (i in 1:nrow(cut)){
  sub <- subset(dat, DogName == cut$DogName[i] & Timestamp < cut$StartDisp[i])

  # If there are fixes before dispersal, take the last entries pack as pack id
  if (nrow(sub) > 0){
    cut$PrePack[i] <- tail(sub$CurrentPack, 1)

  # If there is no pre-dispersal data, we can't assign the pack
  } else {
    cut$PrePack[i] <- NA
  }
}

# Some packs are missing. Let's add them manually (although most of them will
# likely drop out anyways)
cut$PrePack[cut$DogName == "Lupe"]      <- "KM"
cut$PrePack[cut$DogName == "Rattan"]    <- "DB"
cut$PrePack[cut$DogName == "Scorpion"]  <- "KB"
cut$PrePack[cut$DogName == "Stetson"]   <- "MT"

# Remove the column that indicates the dispersal start
cut <- cut %>% select(-StartDisp)

# Finally, in case of the 24 hour dataset we also include the residents. Let's
# Find their pack identity. For this we first need to identify the individuals
# for which we didn't identify the natal pack yet
names1 <- unique(dat$DogName)
names2 <- unique(cut$DogName)

# Find elements in names1 that are not in names2
names <- names1[!names1 %in% names2]

# Find their pack identities and remove any NAs
cut2 <- dat %>%

  # Subset the data to the individuals for which we dont know the pack yet
  subset(DogName %in% names) %>%

  # Keep only the column that indicates the pack
  select(DogName, CurrentPack) %>%

  # Keep only unique entries
  distinct() %>%

  # Remove missing packs
  na.omit() %>%

  # Make sure we have the same columns as in "cut"
  rename(PrePack = CurrentPack)

# Bind the tables together
cut <- rbind(cut, cut2)

# Look at the table
cut

# Extract UD values
extracted <- extrCov(ud, lines)

# Add the layerdates as column names to the extracted values
names(extracted) <- ud_info$Date

# We only want to keep the values from the dates that are clsoest to the steps
# in time
closest2 <- data.frame(closest1 = NA, closest2 = NA)
for (i in seq_along(lines$t1_)){

  # Identify the date of the step
  date <- as.Date(lines$t1_[i])

  # Get closest date 1 (before)
  closest2[i, 1] <- which(abs(date - ud_info$Date) ==
    min(abs(date - ud_info$Date)))[1]

  # Get closest date 2 (after)
  closest2[i, 2] <- which(abs(date - ud_info$Date) ==
    min(abs(date - ud_info$Date)))[2]

  # Transfer the values of the extraction at the closest date into new columns of
  # the lines dataframe
  lines$OtherPack[i]  <- extracted[i, closest2[i, 1]]
}

# Extract HR values
extracted <- extrCov(hr, lines)

# Add layerdates and packnames as column name
names(extracted) <- paste(hr_info$Date, hr_info$Pack, sep = "_")

# We only want to keep the values from the dates that are clsoest to the steps
# in time. Also, we want only those values that belong to the pack to which the
# disperser belonged before dispersal
closest3 <- data.frame(closest1 = NA, closest2 = NA)
for (i in 1:nrow(lines)){

  # Get some information for the respective row (makes indexing much easier)
  dog   <- lines$id[i]
  pack  <- cut$PrePack[cut$DogName == dog]
  hrs   <- subset(hr_info, Pack == pack)
  date  <- as.Date(lines$t1_[i])

  if (pack %in% hr_info$Pack){

    # Get closest date 1 (before)
    closest3[i, 1] <- hrs$LayerNo[which(
      abs(date - hrs$Date) == min(abs(date - hrs$Date))
    )][1]

    # Get closest date 2 (after)
    closest3[i, 2] <- hrs$LayerNo[which(
      abs(date - hrs$Date) == min(abs(date - hrs$Date))
    )][2]

    # Transfer the values of the extraction at the closest date into new columns
    # of the lines dataframe
    lines$Homerange[i]  <- extracted[i, closest3[i, 1]]
  } else {
    lines$Homerange[i] <- NA
  }

  # We also want to know how many home ranges the step crosses. Let's find all
  # homeranges recorded at the closest date
  hrs <- hr_info$LayerNo[which(
    abs(date - hr_info$Date) == min(abs(date - hr_info$Date))
  )]

  # Identify how many of these home ranges were crossed by the step
  lines$NoHomeranges[i] <- rowSums(extracted[, hrs] > 0)[i]
}

############################################################
#### Extraction of Covariates: Distance to Water
############################################################
# Transform the lines to utm
lines <- spTransform(lines, CRS("+init=epsg:32734"))

# To extract distances we need to coerce the lines to a psp object
linesppp <- lapply(lines@lines, function(z){lapply(z@Lines, as.psp)})
linesppp <- do.call("c", linesppp)

# We also need to create points on the lines to extract average distances. We
# can do so by setting a regular distance (100 meters in this case)
linesppp <- lapply(linesppp, pointsOnLines, eps = 100)

# Now load the layer from which we extract the distance to water
wat <- "03_Data/02_CleanData/01_LandCover_Water(Merged).tif" %>%
  stack()

# Extract dates of floodmaps
wat_dates <- names(wat) %>%
  substr(., start = 2, stop = 11) %>%
  as.Date(., format = "%Y.%m.%d")

# Prepare a list that stores a ppp layer for water (Code 1) for each floodmap
watppp <- suppressMessages(
  lapply(1:nlayers(wat), function(x){
    points <- rasterToPoints(wat[[x]], fun = function(z){z == 1})
    points <- as.data.frame(points)
    coordinates(points) <- c("x", "y")
    crs(points) <- crs(wat)
    points <- spTransform(points, CRS("+init=epsg:32734"))
    points <- as(points, "ppp")
    return(points)
  })
)

# Calculate the average distance to water on the ppp object that is closest in
# date to the actual step
lines$DistanceToWater <- suppressMessages(
  mclapply(1:nrow(lines), mc.cores = detectCores() / 2, function(x){
    index <- which.min(abs(as.Date(lines$t1_[x]) - wat_dates))[1]
    distance <- nncross(linesppp[[x]], watppp[[index]])
    distance <- mean(distance$dist)
    return(distance)
    gc()
  }) %>% do.call(rbind, .)
)

# Transform lines back to WGS84
lines <- spTransform(lines, CRS("+init=epsg:4326"))

############################################################
#### Extraction of Covariates: Distance to Roads
############################################################
# Load the distance to roads raster
DistToRoads <- "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToRoads.tif" %>%
  raster()

# Calculate the average distances to the nearest roads
extracted <- extrCov(DistToRoads, lines)

# Put the extracted values into a new column
lines$DistanceToRoads <- extracted[, 1]

# We also want to know whether a step crossed a road or not. We use the
# roads shapefile for this
roads <- "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.shp" %>%
  shapefile()

# Aggregate the shapefile to one object
roads <- gLineMerge(roads)

# Let's use the gDistance function to determine whether a step crosses a road
# Note that we need to transform the layers to UTM for this
extracted <- gDistance(
    spTransform(lines, CRS("+init=epsg:32734"))
  , spTransform(roads, CRS("+init=epsg:32734"))
  , byid = TRUE
) %>% as.vector()

# Make the variable binary, 1 if there is a road crossing (distance = 0)
extracted <- ifelse(extracted == 0, 1, 0)

# Put the extracted values into the lines dataframe
lines$RoadCrossing <- extracted

############################################################
#### Extraction of Covariates: Distance to Humans
############################################################
# Now get the distances to the next human inhabited cell from the base layer
DistToHumans <- "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToHumans.tif" %>%
  raster()

# Calculate distances to nearest human inhabited cell
extracted <- extrCov(DistToHumans, lines)

# Put the extracted values into a new column
lines$DistanceToHumans <- extracted[, 1]

############################################################
#### Storing the Extracted Values
############################################################
# Reorder the columns
names(lines)
lines@data <- dplyr::select(lines@data, c(
    id
  , State
  , step_id_
  , case_
  , x1_
  , x2_
  , y1_
  , y2_
  , t1_
  , t2_
  , dt_
  , sl_
  , ta_
  , absta_
  , tod_
  , Water
  , DistanceToWater
  , Trees
  , Shrubs
  , Protected
  , HumansBase
  , HumansAverage
  , HumansBuff5000
  , DistanceToHumans
  , DistanceToRoads
  , RoadCrossing
  , OtherPack
  , Homerange
  , NoHomeranges
  )
)

# To store the files we need to coerce the duration column to a numeric
lines$dt_ <- as.numeric(lines$dt_)
lines$DistanceToWater <- as.vector(lines$DistanceToWater)

# Save the object to file. Depending on the selected fix schedule the filename
# will be different
if (selection == 24){
  filename <- "00_General_Dispersers_Popecol(SSF24Hours)"
} else {
  filename <- "00_General_Dispersers_Popecol(SSF4Hours)"
}

# Also save the lines to a spatial lines dataframe
writeOGR(lines
  , "03_Data/02_CleanData"
  , filename
  , driver = "ESRI Shapefile"
  , overwrite = TRUE
)

# Let's also store the data to a regular csv. We can use this file to restore
# the original column names since the ESRI shapefiles will store abbreviated
# names
write.csv(lines@data, paste0("03_Data/02_CleanData/", filename, ".csv"))

# Terminate the cluster
endCluster()
