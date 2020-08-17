############################################################
#### Step Selection Function - Extraction of Covariates
############################################################
# Description: In this script we create regular tracks (either 4 or 24 hours)
# and create random steps for the step selection function. Finally, we also
# extract the desired covariates along each of the steps. Note that the script
# needs to be run twice to get the 4 hourly and 24 hourly extractions. Also be
# warned that the extraction takes a good amount of time to complete.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load packages
library(rgdal)
library(raster)
library(data.table)
library(lubridate)
library(amt)
library(tidyverse)
library(spatstat)
library(rgeos)
library(maptools)
library(velox)
library(parallel)

# Load custom functions
source("Functions.r")

# Set seed
set.seed(1234)

# Make use of multicore abilities
beginCluster()

# Load the merged dataset containing all GPS fixes
data <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv()

# Do you want to extract for the 4 hourly or the 24 hourly fixes?
# selection <- 24
selection <- 4

############################################################
#### Data Cleaning (For 4 Hourly Fixes)
############################################################
# For this part we will use the GPS fixes from dispersers only
tracks <- data %>%

  # Subset to the dispersers only (4 hourly fixes are actually only available
  # during dispersal)
  filter(., State == "Disperser") %>%

  # Ultimately, we will only be using fixes that were collected during the
  # following hours: 3, 7, 15, 19, 23 (rounded to 15 minutes)
  filter(.
    , strftime(round_date(Timestamp, "30 minutes")
      , tz      = "UTC"
      , format  = "%H:%M:%S") %in%
      c("03:00:00", "07:00:00", "15:00:00", "19:00:00", "23:00:00")) %>%

  # We also remove popecols fixes before october 2016 because the fixrate scheme
  # was slightly different before that. We still want to conserve the fixes from
  # Abrahms' dispersers
  filter(.,
    Timestamp >= as.POSIXct('2016-10-01 10:00:00') |
    DogName %in% c("Lupe", "Scorpion", "Stetson")
  ) %>%

  # Using the subsetted GPS fixes we prepare tracks
  make_track(.
    , .x    = x
    , .y    = y
    , .t    = Timestamp
    , id    = DogName
    , crs   = CRS("+init=epsg:4326")
    , State = State
  ) %>%

  # Transform the tracks to utm
  transform_coords(CRS("+init=epsg:32734")) %>%

  # Add a column that indicates the time of the day. I follow Gabrieles paper
  # and night starts when the sun is 18 degree below the horizion
  time_of_day(., solar.dep = 18, include.crepuscule = FALSE) %>%

  # Nest the tracks
  nest(data = -"id") %>%

  # Add column that indicates number of fixes
  mutate(NoFixes = data %>%
    lapply(nrow) %>%
    do.call(rbind, .) %>%
    as.vector()) %>%

  # Keep only those individuals where we have more than 3 fixes during disp.
  subset(NoFixes > 3) %>%

  # Turn to a step representation (look up the amt vignette for details)
  mutate(data = map(data, function(x){
    x %>%

      # The function creates steps from the fixes in each tibble row. The option
      # "keep_cols" allows to keep the time of day (tod) column that we added
      steps(keep_cols = "start") %>%

      # Transform the difftime column to a numeric column to avoid that there are
      # heterogeneous units
      mutate(., dt_ = as.numeric(dt_, units = "hours")) %>%

      # We get rid of any steps that take longer than 4 hours, unless they
      # started at 7:00.
      filter(., strftime(
          round_date(t1_, "30 minutes")
        , tz = "UTC"
        , format = "%H:%M:%S") == "07:00:00" | dt_ <= 4.25) %>%

      # Some steps start at 7 but they take longer than 8 hours.
      # We remove those steps as well
      filter(., dt_ < 8.25) %>%

      # Now that we removed some steps we also need to remove some of the turning
      # angles that we are not able to calculate without them
      (function(x){

        # Loop through all rows (the first row has an NA anyways)
        for (i in 2:nrow(x)){

          # Check if the current step starts where the previous step ended
          if (x$t1_[i] != x$t2_[i-1]){

            # If this is not the case, set the turning angle to NA
            x$ta_[i] <- NA
          }
        }
        # Return the resulting tibble
        return(x)
      })
  })) %>%

  # Unnest the tibble
  unnest() %>%

  # Multiply turning angles with negative one (for some reason the package
  # calculates the turning angles conuterclockwise but we want them clockwise)
  mutate(ta_ = ta_ * (-1)) %>%

  # Add a column indicating the absolute turning angle (important for the
  # ellipses). We can use the function we created above.
  mutate(absta_ = absAngle(.)) %>%

  # We can only work with steps for which we have a turning angle. Let's get rid
  # of any steps where the turning angle is NA
  filter(!is.na(ta_)) %>%

  # We also want a unique step id and an indicator that the steps refer to
  # realized steps
  mutate(
      step_id_  = 1:nrow(.)
    , case_     = 1
  )

# Store the 4 hourly tracks
tracks4hours <- tracks

############################################################
#### Data Cleaning (For 24 Hourly Fixes)
############################################################
# For this part we will use the GPS fixes from dispersers only
tracks <- data %>%

  # Again, we only want to extract the data for dispersers
  filter(., State == "Disperser") %>%

  # We only want to keep one fix per day, preferrably one that is close to 07:00
  # Let's be genereous and set a relatively high threshold. Note a tighter
  # threshold only marginally changes the number of fixes we end up with
  filter(.
    , strftime(round_date(Timestamp, "60 minutes")
      , tz      = "UTC"
      , format  = "%H:%M:%S") %in%
      c("03:00:00", "04:00:00", "05:00:00", "06:00:00", "07:00:00", "08:00:00"
        , "09:00:00", "10:00:00", "11:00:00")) %>%

  # Coerce tibble to dataframe
  as.data.frame() %>%

  # If there are multiple fixes on a day, we want to keep only the fix that is
  # closest to 07:00
  resFix(., hours = 24, start = 7, individual = "DogName") %>%

  # Using the subsetted GPS fixes we prepare tracks
  make_track(.
    , .x = x
    , .y = y
    , .t = Timestamp
    , id = DogName
    , crs = CRS("+init=epsg:4326")
    , State = State
  ) %>%

  # Transform the tracks to utm
  transform_coords(CRS("+init=epsg:32734")) %>%

  # Add a column that indicates the time of the day. Note that we request to get
  # dusk and dawn as well. I follow Gabrieles paper and night starts when the sun
  # is 18 degree below the horizion
  time_of_day(., solar.dep = 18, include.crepuscule = FALSE) %>%

  # Nest the tracks
  nest(data = -"id") %>%

  # Add column that indicates number of fixes
  mutate(NoFixes = data %>%
    lapply(., nrow) %>%
    do.call(rbind, .) %>%
    as.vector(.)) %>%

  # Keep only those individuals where we have more than 3 fixes during disp.
  subset(., NoFixes > 3) %>%

  # Turn to a step representation (look up the amt vignette for details)
  mutate(data = map(data, function(x){
    x %>%

      # The function creates steps from the fixes in each tibble row. The option
      # "keep_cols" allows to keep the time of day (tod) column that we added
      steps(., keep_cols = "start") %>%

      # Transform the difftime column to a numeric column to avoid that there are
      # heterogeneous units
      mutate(., dt_ = as.numeric(dt_, units = "hours")) %>%

      # We get rid of any steps that take longer than 29 hours or shorter than
      # 20 hours
      filter(., dt_ <= 32 & dt_ >= 20) %>%

      # Now that we removed some steps we also need to remove some of the turning
      # angles that we are not able to calculate without them
      (function(x){

        # Loop through all rows (the first row has an NA anyways)
        for (i in 2:nrow(x)){

          # Check if the current step starts where the previous step ended
          if (x$t1_[i] != x$t2_[i-1]){

            # If this is not the case, set the turning angle to NA
            x$ta_[i] <- NA
          }
        }
        # Return the resulting tibble
        return(x)
      })
  })) %>%

  # Unnest the tibble
  unnest(.) %>%

  # Multiply turning angles with negative one (for some reason the package
  # calculates the turning angles conuterclockwise but we want them clockwise)
  mutate(., ta_ = ta_ * (-1)) %>%

  # Add a column indicating the absolute turning angle (important for the
  # ellipses). We can use the function we created
  mutate(., absta_ = absAngle(.)) %>%

  # We can only work with steps for which we have a turning angle. Let's get rid
  # of any steps where the turning angle is NA
  filter(., !is.na(ta_)) %>%

  # We also want a unique step id and an indicator that the steps refer to
  # realized steps
  mutate(., step_id_ = 1:nrow(.), case_ = 1)

# Store the 24 hourly tracks
tracks24hours <- tracks

############################################################
#### Generation of Random Steps
############################################################
# Depending on the fixrate schedule the code will look slightly different
if (selection == 24){
  tracks <- tracks24hours
} else {
  tracks <- tracks4hours
}

# We could use the amt package to create random steps. However, the package
# doesn't allow to sample step lengths from the same distribution for all
# individuals. We will thus write our own code. First of all, we need to make
# all steps larger than zero. Only if all steps are greater than zero we can fit
# a gamma distribution. So let's add one meter to all steps that are currently
# zero (only four steps)
tracks$sl_[tracks$sl_ == 0] <- 1

# Fit a Gamma distribution to the observed step lengths
sl <- fit_distr(tracks$sl_, "gamma")

# Write the distribution to file
saveRDS(sl, "03_Data/03_Results/99_GammaDistribution.rds")

# Define the number of random steps
nsteps <- 24

# Prepare an empty list in which we can store the random steps of each individual
randomSteps <- list()

# Now we use the fitted gamma distribution to sample step lengths. For the
# turning angles we follow Avgar et al. 2016 and use a uniform distribution
for (i in 1:nrow(tracks)){

  # Draw random turning angles
  ta_new <- runif(nsteps, min = -pi, max = pi)

  # Draw random step lengths with the fitted parameters
  sl_new <- rgamma(nsteps
    , shape = sl$params$shape
    , scale = sl$params$scale)

  # Put the step lengths and turning angles into a new dataframe and calculate the
  # new endpoints of each random step. We also indicate that the steps are control
  # steps (i.e. 0)
  randomSteps[[i]] <- tracks[rep(i, nsteps), ] %>%
    mutate(.
      , absta_  = absta_ + (ta_new - ta_)
      , ta_     = ta_new
      , sl_     = sl_new
      , case_   = 0
      , x2_     = x1_ + sin(absta_) * sl_
      , y2_     = y1_ + cos(absta_) * sl_
    )
}

# Collapse the list of dataframes into a single dataframe
randomSteps <- do.call(rbind, randomSteps)

# We need to make sure that the absolute turning angle ranges from 0 to 2 * pi
randomSteps$absta_[randomSteps$absta_ > 2 * pi] <-
  randomSteps$absta_[randomSteps$absta_ > 2 * pi] - 2 * pi
randomSteps$absta_[randomSteps$absta_ < 0] <-
  randomSteps$absta_[randomSteps$absta_ < 0] + 2 * pi

# Merge the dataframes of the observed and random steps
ssf <- rbind(tracks, randomSteps) %>%

  # Sort them according to animal id and case/control
  arrange(step_id_, desc(case_)) %>%

  # Make sure that the case/control column is a logical indicator (coerce 1 to
  # TRUE, 0 to FALSE).
  transform(case_ = as.logical(case_))

# Coerce the steps to spatial lines. We can use the function we defined earlier
# for this
lines <- lineTrack(ssf, CRS("+init=epsg:32734"))

# Transform the lines to WGS84
lines <- spTransform(lines, CRS("+init=epsg:4326"))

# Coerce the duration column to a format that can be stored (numeric rather than
# difftime)
lines$dt_ <- as.numeric(lines$dt_)

# Store the object to a shapefile
if (selection == 24){
  filename <- "00_General_Dispersers_Popecol(Random24Hours)"
} else {
  filename <- "00_General_Dispersers_Popecol(Random4Hours)"
}

# Store the stuff into both output directories
writeOGR(lines
  , dsn       = "03_Data/02_CleanData"
  , layer     = filename
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

############################################################
#### Extraction of Covariates: Water
############################################################
# Load the globeland/ori raster layers into a rasterstack
glo <- "03_Data/02_CleanData/01_LandCover_Water(Merged).tif" %>%
  stack()

# Extract the date of each layer
glo_dates <- names(glo) %>%
  substr(., start = 2, stop = 11) %>%
  as.Date(., format = "%Y.%m.%d")

# Now we can extract the percentage cover of Water along each step
extracted <- extrCov(glo, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(glo_dates)

# Let's look at the result
head(extracted)

# We only want to keep the values from the dates that are closest in time to the
# steps
lines$Water <- sapply(1:nrow(lines), function(x){
  index <- which.min(abs(as.Date(lines$t1_[x]) - glo_dates))[1]
  value <- extracted[x, index]
})

############################################################
#### Extraction of Covariates: Trees
############################################################
# Load the treecover map
trees   <- "03_Data/02_CleanData/01_LandCover_TreeCover_MODIS.tif" %>%
  stack()

# Extract the date of each layer
trees_dates <- names(trees) %>%
  substr(., start = 2, stop = 11) %>%
  as.Date(., format = "%Y.%m.%d")

# Extract the average tree cover along each line
extracted <- extrCov(trees, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(trees_dates)

# Let's look at the result
head(extracted)

# Keep only values closest in date
lines$Trees <- sapply(1:nrow(lines), function(x){
  index <- which.min(abs(as.Date(lines$t1_[x]) - trees_dates))[1]
  value <- extracted[x, index]
})

############################################################
#### Extraction of Covariates: Shrubs
############################################################
# Load the shrubcover map
shrubs  <- "03_Data/02_CleanData/01_LandCover_NonTreeVegetation_MODIS.tif" %>%
  stack()

# Extract the date of each layer
shrubs_dates <- names(shrubs) %>%
  substr(., start = 2, stop = 11) %>%
  as.Date(., format = "%Y.%m.%d")

# Extract the average shrub cover in each ellipse
extracted <- extrCov(shrubs, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(shrubs_dates)

# Let's look at the result
head(extracted)

# Keep only values closest in date
lines$Shrubs <- sapply(1:nrow(lines), function(x){
  index <- which.min(abs(as.Date(lines$t1_[x]) - shrubs_dates))[1]
  value <- extracted[x, index]
})

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
