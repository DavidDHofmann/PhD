################################################################################
#### Computing the Flood Extent
################################################################################
# Description: Computing the flood extent of the Okavango Delta at minimum and
# minimum flood

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
setwd(wd)

# Change the working directory
library(terra)
library(plyr)

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
plot(dat, col = "cornflowerblue", main = paste0(names(dat), ": ", round(areas), " km2"))

# Store the numbers to file
areas_round <- round_any(areas, 500)
writeLines(as.character(min(areas_round)), "04_Manuscript/99_FloodExtentMinimum.tex")
writeLines(as.character(max(areas_round)), "04_Manuscript/99_FloodExtentMaximum.tex")
