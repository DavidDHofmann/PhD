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

# Check out all data that needs to be adjusted (all from L1C)
files <- dir(
    path         = dir_l1c
  , include.dirs = T
  , full.names   = T
  , pattern      = ".SAFE$"
)

# Put files into groups of 15
group <- rep(1:ceiling(length(files) / 15), length.out = length(files))

# Function to determine the name of the corrected product
correctedName <- function(x) {
  corr <- gsub(basename(x), pattern = "MSIL1C", replacement = "MSIL2A")
  return(corr)
}

# Go through the groups, move the files to the computer and run the correction
lapply(group, function(x) {

  # Copy the files
  files_sub <- files[group == x]
  files_new <- file.path(tempdir(), basename(files_sub))
  file.copy(files_sub, tempdir(), recursive = T)

  # Run the correction
  sen2cor(files_new
    , outdir    = dir_l2a
    , parallel  = T
    , overwrite = F
    , use_dem   = F
  )

  # Make sure that all files have been correctly converted
  success <- all(file.exists(file.path(dir_l2a, correctedName(files_sub))))

  # If all files have been correctly converted, remove the originals as well as
  # the copies in the tempdir
  unlink(files_new, recursive = T)
  unlink(files_sub, recursive = T)

})

################################################################################
#### Testing
################################################################################
################################################################################
#### WRITE CUSTOM FUNCTIONS!
################################################################################
plotSCL <- function(x) {


}



setwd("/home/david/Schreibtisch/Sentinel")

files <- sample(dir(path = "L2A", include.dirs = T, full.names = T), size = 1)
file <- s2_translate(files)
mask <- files
mask <- s2_translate(mask, prod_type = "SCL")

library(terra)
rf <- rast(file)
rm <- rast(mask)
freq(rm)

# Need to disaggregate the mask to match the other file
rm <- disagg(rm, fact = 2)
compareGeom(rf, rm)

# Mask has the following classes (check here: https://sen2r.ranghetti.info/articles/outstructure#accessory-layers)
classes <- c("NoData", "Saturated or Defective", "Dark Area Pixels", "Cloud Shadows", "Vegetation", "Bare Soils", "Water", "Clouds Low Probability / Unclassified", "Clouds Medium Probability", "Clouds High Probability", "Cirrus", "Snow / Ice")
length(classes)
levels(rm) <- classes
plot(rm)
freq(rm)
test <- mask(rf, rm, maskvalue = c(0, 3, 8, 9, 10), updatevalue = NA)
