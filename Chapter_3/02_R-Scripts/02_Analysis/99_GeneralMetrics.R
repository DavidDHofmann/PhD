################################################################################
#### Some General Metrics to Include in the Manuscript
################################################################################
# Description: Calculation of some general metrics that will be included in the
# manuscript

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load packages
library(tidyverse)    # For data wrangling
library(lubridate)    # To handle dates
library(terra)        # To handle spatial data

# Custom functions
source("02_R-Scripts/00_Functions.R")

# Store the metrics into a named vector
metrics <- c()

# Load (subsampled) dispersal data
dat <- read_csv("03_Data/02_CleanData/DispersersSubsampled.csv")
dat <- subset(dat, State == "Disperser")

# Identify all dogs
names <- dat %>%
  pull(ID) %>%
  unique()

# Load the step selection data
ssf <- "03_Data/02_CleanData/SSF.csv" %>%
  read_csv() %>%
  subset(case == 1 & ID %in% names)

# Average step-length
mean(ssf$sl)

# Also load the non-subsampled data (to compute acquisition rate)
all <- read_csv("03_Data/02_CleanData/Dispersers.csv")
all <- subset(all, State == "Disperser" & Source == "Popecol")

################################################################################
#### Number of Coalitions in Original Data
################################################################################
# Number of males
metrics[1] <- dat %>%
  subset(Sex == "M") %>%
  select(ID, Sex) %>%
  distinct() %>%
  nrow()

# Number of females
metrics[2] <- dat %>%
  subset(Sex == "F") %>%
  select(ID, Sex) %>%
  distinct() %>%
  nrow()

# Number in total
metrics[3] <- dat %>%
  select(ID, Sex) %>%
  distinct() %>%
  nrow()

# Put names
names(metrics)[1:3] <- c("CollarsMales", "CollarsFemales", "CollarsTotal")

################################################################################
#### Number of GPS datapoints
################################################################################
# In total
metrics[4] <- dat %>%
  nrow() %>%
  round() %>%
  format(big.mark = ",")

# Mean per coalition +- sd
metrics[5] <- dat %>%
  count(ID) %>%
  select(n) %>%
  summarize(mean = round(mean(n)), sd = round(sd(n))) %>%
  mutate(comb = paste0(mean, " $\\pm$ ", sd)) %>%
  pull(comb)

# Put names
names(metrics)[4:5] <- c("FixesTotal", "FixesMeanSD")

################################################################################
#### Expanse of the Study Area
################################################################################
# Extent of the study area
exte <- ext(c(22, 27, -20.7, -17.9))

# Compute expanse
metrics[6] <- exte %>%
  as.polygons(crs = "+init=epsg:4326") %>%
  expanse(unit = "km") %>%
  round(-4) %>%
  format(big.mark = ",")
names(metrics)[6] <- c("SizeStudyArea")

################################################################################
#### Percentage Protected
################################################################################
# Load protected areas
prot <- vect("03_Data/02_CleanData/Protected.gpkg")
prot <- crop(prot, exte)
back <- rast(prot, res = 1 / 111)

# Rasterize them
prot <- rasterize(prot, back, background = 0)

# Compute fraction
metrics[7] <- prot %>%
  expanse(byValue = T, unit = "km") %>%
  as.data.frame() %>%
  mutate(Percentage = area / sum(area) * 100) %>%
  subset(value == 1) %>%
  pull(Percentage) %>%
  round()
names(metrics)[7] <- c("PercentageProtected")

################################################################################
#### Center of Study Area
################################################################################
# Compute center point
x <- mean(c(exte[1], exte[2]))
y <- mean(c(exte[3], exte[4]))

# Convert to degrees and minutes
x_degs <- abs(floor(x))
x_mins <- round(x%%1 * 60)
y_degs <- abs(floor(y))
y_mins <- round(y%%1 * 60)

# Print nicely
x  <- paste0(x_degs, "\\degree", x_mins, "'E")
y  <- paste0(y_degs, "\\degree", y_mins, "'S")
xy <- paste(x, y)
metrics[8] <- xy
names(metrics)[8] <- c("StudyAreaCenter")

################################################################################
#### Seasonal Covariates
################################################################################
# Reload the data
load("03_Data/03_Results/SeasonalCovariates.rds")

# Compute metrics of interest for precipitation
metrics[9:10] <- means_prec %>%
  arrange(Precipitation) %>%
  pull(Precipitation) %>%
  range() %>%
  round(-1)
metrics[11] <- means_prec %>%
  pull(Precipitation) %>%
  sum() %>%
  round(-1)
names(metrics)[9:11] <- c("MinimumPrecipitation", "MaximumPrecipitation", "TotalPrecipitation")

# Compute metrics of interest for temperature
metrics[12:13] <- means_temp %>%
  arrange(Temperature) %>%
  pull(Temperature) %>%
  range() %>%
  round()
names(metrics)[12:13] <- c("MinimumTemperature", "MaximumTemperature")

################################################################################
#### Acquisition Rate
################################################################################
# Compute number of successfull fixes by collar and Dog
metrics[14] <- all %>%
  dplyr::select(ID, CollarID, x) %>%
  mutate(Valid = if_else(is.na(x), "invalid", "valid")) %>%
  count(ID, CollarID, Valid) %>%
  pivot_wider(names_from = Valid, values_from = n, values_fill = 0) %>%
  mutate(AcquisitionRate = valid / (valid + invalid) * 100) %>%
  ungroup() %>%
  summarize(
      MeanAcquisitionRate = round(mean(AcquisitionRate), 0)
    , SDAcquisitionRate   = round(sd(AcquisitionRate))
  ) %>%
  as.character() %>%
  paste(., collapse = " $\\pm$ ")
names(metrics)[14] <- c("AcquisitionRate")

################################################################################
#### Number of Coalitions in Step Selection Data
################################################################################
# Number of males
metrics[15] <- ssf %>%
  subset(Sex == "M") %>%
  select(ID, Sex) %>%
  distinct() %>%
  nrow()

# Number of females
metrics[16] <- ssf %>%
  subset(Sex == "F") %>%
  select(ID, Sex) %>%
  distinct() %>%
  nrow()

# Number in total
metrics[17] <- ssf %>%
  select(ID, Sex) %>%
  distinct() %>%
  nrow()

# Put names
names(metrics)[15:17] <- c("SSFMales", "SSFFemales", "SSFTotal")

################################################################################
#### Number of Steps
################################################################################
# Load step-selection data
metrics[18] <- ssf %>%
  nrow() %>%
  format(big.mark = ",")
metrics[19] <-ssf %>%
  count(ID) %>%
  summarize(
      MeanNumberSteps = round(mean(n), 0)
    , SDNumberSteps   = round(sd(n))
  ) %>%
  as.character() %>%
  paste(., collapse = " $\\pm$ ")
names(metrics)[18:19] <- c("StepsTotal", "StepsMeanSD")

################################################################################
#### Number of steps by Season
################################################################################
# Derive the season for each GPS fix
ssf$Season1 <- getRainySeason(ssf$Timestamp)
ssf$Season2 <- getHerbivoreSeason(ssf$Timestamp)

# Get number of fixes in different seasons
metrics[20:21] <- ssf %>%
  count(Season1) %>%
  pull(n) %>%
  format(big.mark = ",")
metrics[22:23] <- ssf %>%
  count(Season2) %>%
  pull(n) %>%
  format(big.mark = ",")

# Put names
names(metrics)[20:21] <- c("NumberStepsDry", "NumberStepsWet")
names(metrics)[22:23] <- c("NumberStepsConcentrated", "NumberStepsDispersed")

# Get percentage of fixes in different seasons
metrics[24:25] <- ssf %>%
  count(Season1) %>%
  mutate(Freq = n / sum(n) * 100) %>%
  pull(Freq) %>%
  round()
metrics[26:27] <- ssf %>%
  count(Season2) %>%
  mutate(Freq = n / sum(n) * 100) %>%
  pull(Freq) %>%
  round()

# Put names
names(metrics)[24:25] <- c("PercentageStepsDry", "PercentageStepsWet")
names(metrics)[26:27] <- c("PercentageStepsConcentrated", "PercentageStepsDispersed")

################################################################################
#### Training Polygons
################################################################################
pols <- vect("03_Data/02_CleanData/TrainingClasses.gpkg")

# Number of polygons by class
class <- "03_Data/02_CleanData/TrainingClasses.gpkg" %>%
  vect() %>%
  values() %>%
  count(Class)

# Put into mterics
metrics[28:30] <- class$n
metrics[31]    <- sum(class$n)

# Give nice names
names(metrics)[28:30] <- paste0("Training", class$Class)
names(metrics)[31] <- "TrainingTotal"

# Also extract the dates at which those polygons were generated.
metrics[32] <- pols %>%
  values() %>%
  pull(Date) %>%
  unique() %>%
  length()
names(metrics)[32] <- "TrainingDates"

# Compute the relative size of wet-pans compared to other training classes
pols$Area <- expanse(pols)
metrics[33] <- pols %>%
  values() %>%
  group_by(Class) %>%
  summarize(Area = sum(Area), .groups = "drop") %>%
  mutate(AreaPercentage = round(Area / sum(Area) * 100, 2)) %>%
  subset(Class == "Wetpan") %>%
  pull(AreaPercentage)
names(metrics)[33] <- "PercentageTrainingWetpan"

################################################################################
#### Range of Dates of GPS timestamps
################################################################################
metrics[34:35] <- year(range(all$Timestamp))
names(metrics)[34:35] <- c("GPSFromYear", "GPSToYear")

################################################################################
#### Step Distribution
################################################################################
metrics[36:37] <- "03_Data/03_Results/StepLengthDistribution.rds" %>%
  read_rds() %>%
  unlist() %>%
  round(2)
names(metrics)[36:37] <- c("GammaShape", "GammaScale")

################################################################################
#### Sentinel 2 Tiles
################################################################################
metrics[38:39] <- "/media/david/Elements/Todownload.rds" %>%
  read_rds() %>%
  count(level) %>%
  pull(n) %>%
  prettyNum(big.mark = ",")
metrics[40] <- "/media/david/Elements/Todownload.rds" %>%
  read_rds() %>%
  nrow() %>%
  prettyNum(big.mark = ",")
names(metrics)[38:39] <- c("NumberTiles1C", "NumberTiles2A")
names(metrics)[40] <- c("NumberTilesTotal")

################################################################################
#### Store Metrics
################################################################################
# Go through each metric and store it as tex
dir.create("04_Manuscript/GeneralMetrics", showWarnings = F)
for (i in 1:length(metrics)) {
  name <- paste0("04_Manuscript/GeneralMetrics/", names(metrics)[i], ".tex")
  writeLines(metrics[[i]], name)
}
