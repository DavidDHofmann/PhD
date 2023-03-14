############################################################
#### Step Selection Function - Social Landscape Analysis
############################################################
# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Change the working directory
input1 <- "/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/01_SmallExtent"
input2 <- "/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/00_LargeExtent"
setwd(input1)

# Load required packages
library(raster)
library(rgeos)
library(data.table)
library(tidyverse)
library(glmmTMB)
library(lubridate)
library(viridis)

# Load custom functions
source("/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/Functions.r")

############################################################
#### Loading and Cropping Data
############################################################
# Load the 4 hourly fixes (this will take some time). Note that we need to load
# the shapefile, because we want to crop the data to a smaller extent
ssf <- shapefile("00_General_Dispersers_Popecol(SSF4Hours)")

# We will also load the csv because it contains the data with nicer column names
# Additionally, through read_csv we are able to get the columns to the correct
# data types
ssf2 <- "00_General_Dispersers_Popecol(SSF4Hours).csv" %>%

  # We need to manually set some of the column types
  read_csv(., col_types = cols(
        Homerange     = "d"
      , OtherPack     = "d"
      , id            = "f"
      , tod_          = "f"
      , RoadCrossing  = "f"
      , State         = "f")) %>%

  # Remove undesired columns
  select(-c("X1", "x1_", "x2_", "y1_", "y2_", "dt_")) %>%

  # Rename social covariates to have more sensible names
  rename(.,
      OwnHomerange      = Homerange
    , ForeignHomeranges = NoHomeranges
    , ProbEncounter     = OtherPack
  )

# Currently, the ForeignHomeranges still contains the own homerange of the
# wild dog. Let's subtract 1 whenever the wild dog is in his own HR
ssf2$ForeignHomeranges[ssf2$OwnHomerange > 0 & !is.na(ssf2$OwnHomerange)] <-
  ssf2$ForeignHomeranges[ssf2$OwnHomerange > 0 & !is.na(ssf2$OwnHomerange)] - 1

# Check if the two are still overly correlated
cor(ssf2$OwnHomerange, ssf2$ForeignHomeranges, "pairwise.complete.obs")

# Scale Covariates
ssf2 <- ssf2 %>%
  transform(
      Water                 = scale(Water)
    , DistanceToWater       = scale(sqrt(DistanceToWater))
    , Trees                 = scale(Trees)
    , Shrubs                = scale(Shrubs)
    , Protected             = scale(Protected)
    , HumansBase            = scale(HumansBase)
    , HumansAverage         = scale(HumansAverage)
    , HumansBuff5000        = scale(HumansBuff5000)
    , DistanceToHumans      = scale(sqrt(DistanceToHumans))
    , DistanceToRoads       = scale(sqrt(DistanceToRoads))
    , OwnHomerange          = scale(OwnHomerange)
    , ForeignHomeranges     = scale(ForeignHomeranges)
    , ProbEncounter         = scale(ProbEncounter)
  )

# Combine the two datasets
ssf@data <- ssf2

# For the analysis of the social landscape we will be using only those fixes
# that were recorded within a specific bounding box. Let's create the
# corresponding bounding box
bbox <- extent(c(23.25, 24, -19.8, -19.05)) %>% as(., "SpatialPolygons")

# We might also add a crs to the bounding box
crs(bbox) <- CRS("+init=epsg:4326")

# We will only keep those steps that fall within the bounding box
keep <- bbox %>%

  # Check whether a step is contained within the boundary or not
  gContains(ssf, byid = TRUE) %>%

  # Coerce to a dataframe
  as.data.frame() %>%

  # Select only the column with the logicals
  .[["1"]]

# Subset the steps to the ones that are contained within the polygon
ssf <- ssf[keep, ]

# Note that we don't have the full number of random steps for all steps anymore
# (since they don't fall within the bbox). We need to get rid of all steps for
# which we are lacking the desired number of random steps. Let's first identify
# the step_id_ of the steps that need to be removed
remove <- ssf@data %>%

  # Groupd the steps by their step_id_
  group_by(., step_id_) %>%

  # Identify the number of observations for each step_id_
  summarize(., Observations = (length(step_id_))) %>%

  # Subset to the steps with less than 25 observations
  subset(., Observations < 25) %>%

  # Select the column indicating the step_id_ that nees to be removed
  .[["step_id_"]]

# Remove those step ids from the original dataframe
ssf <- subset(ssf, !(step_id_ %in% remove))

# Make sure there are only step_id_ with 25 observations
table(ssf$step_id_) %>%
  as.data.frame(.) %>%
  .[["Freq"]] %>%
  table(.)

# Plot the remaining steps together with the bounding box
plot(ssf)
plot(bbox, add = TRUE, border = "red")

# Now we don't need the spatial information anymore and can get rid of it
ssf <- ssf@data

# Let's check the dates for which we prepared social landscapes
ud_dates <- "05_SocialFeatures_SocialLandscape(UtilisationDistributions).grd" %>%

  # Load the utilisation distributions into a stack
  stack(.) %>%

  # Identify the layernames
  names(.) %>%

  # Extract the dates
  substr(., start = 2, stop = 11) %>%

  # Turn into true dates
  as.Date(., format = "%Y.%m.%d") %>%

  # Check the range
  range(.)

# Let's subset to the fixes that were taken within that timeframe
ssf <- subset(ssf, t1_ >= ud_dates[1] - days(30) & t2_ <= ud_dates[2] + days(30))

# We also need to make sure that we have information about the homerange of the
# individuals
ssf <- subset(ssf, !is.na(OwnHomerange))

# Let's see which individuals we have left for analysis and how many steps there
# are left for each of them
ssf %>% group_by(id) %>% summarize(NoFixes = n())

############################################################
#### Run GLMMTMB Model
############################################################
# We will run our most parsimonious model, but now inlcude the social landscape
# as covariates as well. Let's prepare the model call of the standard model
models <- readRDS("99_ModelSelection.rds")

# Get the formula of the best model
form <- models$Formula[[1]]

# Now update the formula and include the social landscape covariates
form <- update(form, ~ .
  + ProbEncounter
  + OwnHomerange
  + ForeignHomeranges
  + ProbEncounter:OwnHomerange
  + ProbEncounter:ForeignHomeranges
  + (0 + ProbEncounter|id)
  + (0 + OwnHomerange|id)
  + (0 + ForeignHomeranges|id)
)

# And rund the model
mod <- glmm_clogit(form, ssf)

# Look at the summary
summary(mod)

# Plot the coefficients
showCoeffs(getCoeffs(mod)[-1, ])

# Store the result for later
setwd(input2)
write_rds(mod, "99_SocialModel.rds")
setwd(input1)
write_rds(mod, "99_SocialModel.rds")
setwd(input1)

# Reload the results
models <- read_rds("99_MovementModel.rds")

############################################################
#### Plot the Interactions
############################################################
# Plot the interactions and store the plots
png("99_SocialInteractions(ProbencounterOwnhomerange).png"
  , width = 1080
  , height = 1080
  , bg = "transparent"
  , pointsize = 30
)
visInt(mod, "ProbEncounter", "OwnHomerange", label.color = "white")
dev.off()

png("99_SocialInteractions(ProbencounterForeignhomeranges).png"
  , width = 1080
  , height = 1080
  , bg = "transparent"
  , pointsize = 30
)
visInt(mod, "ProbEncounter", "ForeignHomeranges", label.color = "white")
dev.off()
