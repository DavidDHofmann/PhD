################################################################################
#### Computing the Flood Extent
################################################################################
# Description: Computing the flood extent of the Okavango Delta at minimum and
# minimum flood

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Change the working directory
library(terra)    # To handle spatial data
library(plyr)     # To wrangle

# Load the watermaps
dat <- rast("03_Data/02_CleanData/WaterCover.tif")
dat <- dat[[c("min", "max")]]

# We only want to focus on the Okavango delta, thus, let's load a shapefile of
# it
oka <- vect("03_Data/02_CleanData/MajorWaters.shp")
oka <- oka[oka$name == "Okavango Delta", ]

# Buffer it slightly
oka <- buffer(oka, width = 10000)
oka <- fillHoles(oka)

# Let's take a look at the overlap
plot(dat[["max"]])
plot(oka, add = T)

# This looks good. We now mask anything outside that polygon
dat <- mask(dat, oka)
dat <- trim(dat)
dat <- subst(dat, 0, NA)

# Now let's compute the total area covered by the flood in the two cases
areas <- expanse(dat, unit = "km")
plot(dat, col = "cornflowerblue", main = paste0(names(dat), ": ", round(areas$area), " km2"))

# Store the numbers to file
areas_round <- round_any(areas$area, 500)
areas_round <- format(areas_round, big.mark = ",")
areas_round

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/06_ExtentFlood.rds")

# Print to terminal
cat("Done :)\n")
