################################################################################
#### Sentinel Data Processing L1C to L2A
################################################################################
# Description: Processing the Sentinel 2 L1C data to achieve L2A products.

# Clear R's brain
rm(list = ls())

# Load required packages
library(sen2r)     # To automate the correction of Sentinel 2 data
library(parallel)  # To check for the number of cores

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Specify the directories to the sentinel files
dir_l1c <- "/media/david/Elements/L1C"
dir_l2a <- "/media/david/Elements/L2A"
dir_l1c <- "/home/david/Schreibtisch/Sentinel/L1C"
dir_l2a <- "/home/david/Schreibtisch/Sentinel/L2A"

# Check out all data that needs to be adjusted
files <- dir(
    path         = dir_l1c
  , include.dirs = T
  , full.names   = T
  , pattern      = ".SAFE$"
)

# Run the correction
sen2cor(files
  , outdir    = dir_l2a
  , parallel  = detectCores() - 1
  , overwrite = F
  , use_dem   = F
)

################################################################################
#### Testing
################################################################################
################################################################################
#### WRITE CUSTOM FUNCTIONS!
################################################################################
plotSCL <- function(x) {


}



setwd("/home/david/Schreibtisch/Sentinel")

file <- "/home/david/Schreibtisch/Sentinel/L2A/S2B_MSIL2A_20210929T080719_N0301_R078_T34KGC_20210929T105251.SAFE"
file <- s2_translate(file)
mask <- "/home/david/Schreibtisch/Sentinel/L2A/S2B_MSIL2A_20210929T080719_N0301_R078_T34KGC_20210929T105251.SAFE"
mask <- s2_translate(mask, prod_type = "SCL")



rf <- rast(file)
rm <- rast(mask)

# Need to disaggregate the mask to match the other file
rm <- disagg(rm, fact = 2)
compareGeom(rf, rm)

# Mask has the following classes (check here: https://sen2r.ranghetti.info/articles/outstructure#accessory-layers)
classes <- c("NoData", "Saturated or Defective", "Dark Area Pixels", "Cloud Shadows", "Vegetation", "Bare Soils", "Water", "Clouds Low Probability / Unclassified", "Clouds Medium Probability", "Clouds High Probability", "Cirrus", "Snow / Ice")
length(classes)
levels(rm) <- classes
plot(rm)


final <- s2_mask(file, mask, mask_type = "scl_0_3_8_9")


r <- rast(ncols = 10, nrows = 10)
values(r) <- sample(3, ncell(r), replace = T)
cls <- c("Bare", "Water", "Pan", "Vegetation")
r <- r - 1
levels(r) <- cls
plot(r)

set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE)
is.factor(r)
cls <- c("forest", "water", "urban")
plot(r)
# make the raster start at zero
x <- r - 1
levels(x) <- cls
names(x) <- "land cover"
is.factor(x)
x
plot(x, col=c("green", "blue", "light gray"))
text(x, digits=3, cex=.75, halo=TRUE)
# raster starts at 3
x <- r + 2
is.factor(x)
# approach 1
levels(x) <- c("", "", "", "forest", "water", "urban")
# approach 2, also showing the use of two categories
d <- data.frame(id=3:5, cover=cls, letters=letters[1:3], value=10:12)
levels(x) <- d
x
## switch categories
cats(x, 1)
# get current index
activeCat(x)
# set index
activeCat(x) <- 3
plot(x, col=c("green", "blue", "light gray"))
text(x, digits=3, cex=.75, halo=TRUE)
r <- as.numeric(x)
r
activeCat(x) <- 2
p <- as.polygons(x)
plot(p, "letters", col=c("green", "blue", "light gray"))














# Check out all downloaded files
files <- dir(path = "L1C", include.dirs = T, full.names = T)

# Read metadata of the files
meta <- safe_getMetadata(files)

# Run sen2cor
sen2cor(files, outdir = "L2A", use_dem = T, parallel = detectCores() - 1, overwrite = F)

# Detect L2A files
files <- dir(path = "L2A", include.dirs = T, full.names = T, pattern = "SAFE")
print(files)

# Translate
lapply(files, s2_translate)

# Apply cloud mask to scl classes (check here
# https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR)
files_vrt <- dir(path = "L2A", pattern = ".vrt", full.names = T)
s2_mask(infiles = files_vrt[1], maskfiles = files_vrt[1], mask_type = "scl_0_3_8_9")

library(sf)
foot <- st_as_sfc(meta$footprint)
foot <- st_as_sf(foot)
plot(foot)
