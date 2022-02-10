################################################################################
#### Pan Mapping Training: Comparison between Landsat 8 and Sentinel 2
################################################################################
# Description: Here, we'll download Landsat 8 and Sentinel 2 imagery and train a
# classifier that is supposed to detect wet pans. We will than compare the
# performance of landsat versus sentinel. We will later use the trained
# classifier together with the chosen product (either Sentinel 2 or Landsat 8)
# to predict the locations of pans using imagery that is updated every month or
# so.

# Clear R's brain
rm(list = ls())

# Load required packages
library(rgdal)         # To handle spatial data
library(sf)            # To handle spatial data
library(raster)        # To handle spatial data
library(terra)         # To handle spatial data
library(lubridate)     # To handle dates
library(tidyverse)     # For data wrangling
library(parallel)      # To identify the number of cores
library(rgee)          # To download Landsat 8 data
library(sen2r)         # To download Sentinel 2 data
library(ggridges)      # For ridgelines in ggplot
library(rpart)         # To train a classifier
library(randomForest)  # To train a classifier
library(rpart.plot)    # To visualize classifier
library(caret)         # To assess variable importance
library(ggpubr)        # To arrange multiple ggplots
library(ggstance)      # To plot pointranges

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Create necessary directories
dir.create("03_Data/01_RawData/LANDSAT", showWarnings = F)
dir.create("03_Data/01_RawData/SENTINEL", showWarnings = F)
dir.create("03_Data/01_RawData/SENTINEL/L1C", showWarnings = F)
dir.create("03_Data/01_RawData/SENTINEL/L2A", showWarnings = F)
outdir_l1c <- "03_Data/01_RawData/SENTINEL/L1C"
outdir_l2a <- "03_Data/01_RawData/SENTINEL/L2A"

# Function to compute the normalized difference (nd) index of two bands
nd <- function(img, band_x, band_y) {
  x <- img[[band_x]]
  y <- img[[band_y]]
  nd <- (x - y) / (x + y)
  return(nd)
}

################################################################################
#### Prepare Training Data
################################################################################
# Identify ground truth data
files <- dir(
    path       = "03_Data/01_RawData/GABRIELE"
  , pattern    = ".kml$"
  , full.names = T
)

# Go through the files and load the different land cover training data
classes <- lapply(files, function(x) {

  # List available land cover categories and keep only the ones of interest
  cats <- ogrListLayers(x)
  cats <- cats[cats %in% c("Dry_Pans", "Wet_Pans", "Pools", "Dry_Land")]
  date <- ymd(basename(x))

  # Go through the categories and load data on the different categories
  cats <- lapply(cats, function(y) {
    cat <- readOGR(x, layer = y, D3_if_2D3D_points = T)
    cat@data <- data.frame(Date = rep(date, nrow(cat)), Class = rep(y, nrow(cat)))
    return(cat)
  })

  # Put them together
  cats <- do.call(rbind, cats)

  # Return the result
  return(cats)

}) %>% do.call(rbind, .)

# Rename the classes nicely
classes$Class[classes$Class == "Dry_Land"] <- "Dryland"
classes$Class[classes$Class == "Dry_Pans"] <- "Drypan"
classes$Class[classes$Class == "Wet_Pans"] <- "Wetpan"
classes$Class[classes$Class == "Pools"]    <- "Water"

# Drypans are also dryland
classes$Class[classes$Class == "Drypan"] <- "Dryland"

# Check out the resulting dataframe
as.data.frame(classes)
table(classes$Class)

# Visualize the data
ggplot(st_as_sf(classes), aes(col = Class, fill = Class)) +
  geom_sf(alpha = 0.5) +
  facet_wrap("Date") +
  theme_minimal()

# Check the unique dates for which we have training data. We will have to
# download satellite data for those dates as well
dates <- unique(classes$Date)

################################################################################
#### Download Landsat 8 Data
################################################################################
# Specify correct python environment for rgee
ee_install_set_pyenv(
    py_path = "/home/david/miniconda3/envs/rgee/bin/python"
  , py_env  = "rgee"
)

# Make sure we have all installed for rgee
ee_Initialize()

# Function to mask cloudy pixels (check out the pixel descriptions here:
# https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_SR)
# Code for the mask is based on this example:
# https://r-spatial.github.io/rgee/articles/rgee03.html
maskLandsat <- function(image) {
  mask_shadow <- bitwShiftL(1, 3)
  mask_clouds <- bitwShiftL(1, 5)
  qa <- image$select("pixel_qa")
  mask <- qa$bitwiseAnd(mask_shadow)$eq(0)$And(qa$bitwiseAnd(mask_clouds)$eq(0))
  image$updateMask(mask)$
    divide(10000)$
    select("B[0-9]*")$
    copyProperties(image, list("system:time_start"))
}

# Loop through the dates and download landsat data
for (i in 1:length(dates)) {

  # Subset to the classes of interest
  classes_sub <- subset(classes, Date == dates[i])

  # Define the area of interest
  aoi <- ee$Geometry$Polygon(
    list(
        c(extent(classes_sub)[1], extent(classes_sub)[3])
      , c(extent(classes_sub)[1], extent(classes_sub)[4])
      , c(extent(classes_sub)[2], extent(classes_sub)[4])
      , c(extent(classes_sub)[2], extent(classes_sub)[3])
    )
  )

  # Specify filename
  filename <- paste0("03_Data/01_RawData/LANDSAT/", dates[i], ".tif")

  # If the file doesn't exist, download it
  if (!file.exists(filename)) {

    # Query satellite data
    collection <- ee$ImageCollection("LANDSAT/LC08/C01/T1_SR")$
      filterDate(as.character(dates[i] - days(10)), as.character(dates[i] + days(10)))$
      filterBounds(aoi)$
      map(maskLandsat)$
      median()$
      clip(aoi)

    # Check it
    # ee_print(collection)

    # Download
    tempfile <- tempfile(fileext = ".tif")
    ee_as_raster(collection, region = aoi, dsn = tempfile, scale = 30)

    # Load the data
    dat <- rast(tempfile)

    # Compute indices
    ndvi <- nd(dat, "B5", "B4")
    ndwi <- nd(dat, "B3", "B5")
    ndmi <- nd(dat, "B5", "B6")
    ndsi <- nd(dat, "B3", "B6")
    best <- nd(dat, "B3", "B7")

    # Put all together into a single stack
    all <- c(dat, ndvi, ndwi, ndmi, ndsi, best)
    names(all) <- c(names(dat), "ndvi", "ndwi", "ndmi", "ndsi", "best")

    # Store the raster to file again
    writeRaster(all, filename, overwrite = T)
  }
}

################################################################################
#### Download Sentinel 2 Data
################################################################################
# Login to scihub
write_scihub_login("dodx9", "Scihubbuster69_")

# Identify the files that we need to download
todownload <- lapply(dates, function(x) {

  # Subset to the classes of interest
  classes_sub <- subset(classes, Date == x)

  # Define a polygon that indicates the area for which we need to download sentinel
  # 2 data
  aoi <- extent(classes_sub)
  aoi <- as(aoi, "SpatialPolygons")
  crs(aoi) <- "+init=epsg:4326"

  # Find products for the range of dates
  todownload <- s2_list(
      spatial_extent = st_as_sf(aoi)
    , time_interval  = c(x - days(10), x + days(10))
    , level          = "auto"
    , server         = c("scihub", "gcloud")
  )

  # Check the metadata and create a regular date column
  meta <- as.data.frame(todownload)
  meta <- mutate(meta, Date = as.Date(sensing_datetime))

  # Remove data with high cloud cover
  index <- which(meta$clouds < 5)
  todownload <- todownload[index]
  meta <- meta[index, ]

  # Find the image(s) that is/are closest in date to the input date
  index <- which(abs(meta$Date - x) == min(abs(meta$Date - x)))
  todownload <- todownload[index]
  meta <- meta[index, ]

  # Put everything into a tibble
  todownload  <- tibble(todownload, meta)
  todownload $footprint <- NULL

  # Assign the date from the training data
  todownload $Date <- x

  # Return the respective files
  return(todownload )

}) %>% do.call(rbind, .)

# Specify output directory depending on the product level
todownload$outdir <- ifelse(todownload$level == "1C", outdir_l1c, outdir_l2a)

# Store the download table to file
write_rds(todownload, "03_Data/01_RawData/SENTINEL/Todownload.rds")
todownload <- read_rds("03_Data/01_RawData/SENTINEL/Todownload.rds")

# Subset to the files that we haven't downloaded yet
exists <- dir.exists(paste0(todownload$outdir, "/", todownload$name))
todownload_sub <- subset(todownload, !exists)

# Download the files
if (nrow(todownload_sub) > 0) {

  # Order the products if necessary (need to wait after this)
  s2_order(todownload_sub$todownload)

  # Keep only the ones available
  online <- safe_is_online(todownload_sub$todownload)
  todownload_sub <- todownload_sub[online, ]

  # Download the products
  lapply(1:nrow(todownload_sub), function(x) {
    tryCatch({
      s2_download(todownload_sub$todownload[x], outdir = todownload_sub$outdir[x])
    }, error = function(e) {return(e)})
  })
}

# Run correction on the L1C products
files <- dir(path = outdir_l1c, include.dirs = T, full.names = T)
sen2cor(files
  , outdir    = outdir_l2a
  , parallel  = T
  , overwrite = F
  , use_dem   = F
)

# Determine final output directory and filenames
todownload$outdir_final <- outdir_l2a
todownload$name_final <- gsub(todownload$name
  , pattern     = "MSIL1C"
  , replacement = "MSIL2A"
)

# Loop through the dates and process the respective images
maps <- lapply(dates, function(x) {

  # Prepare filename for final image
  filename <- paste0("03_Data/01_RawData/SENTINEL/", x, ".tif")

  # Only need to do this if the file doesn't exist yet
  if (!file.exists(filename)) {

    # Specify area of interest
    aoi <- extent(subset(classes, Date == x))
    aoi <- as(aoi, "SpatialPolygons")
    crs(aoi) <- "+init=epsg:4326"

    # Create filepath to the respective files
    stitch <- subset(todownload, Date == x)
    stitch <- file.path(wd, stitch$outdir_final, stitch$name_final)

    # Loop through them and prepare masked image
    stitch <- lapply(stitch, function(y) {
      vrt <- s2_translate(y, format = "VRT", outdir = tempdir())
      mask <- s2_translate(y, prod_type = "SCL", outdir = tempdir())

      # Load them
      rf <- rast(vrt)
      rm <- rast(mask)
      if (nrow(rf) != nrow(rm)) {
        rm <- disagg(rm, fact = 2)
      }

      # Assign classes to mask
      cat <- c("NoData", "Saturated or Defective", "Dark Area Pixels", "Cloud Shadows", "Vegetation", "Bare Soils", "Water", "Clouds Low Probability / Unclassified", "Clouds Medium Probability", "Clouds High Probability", "Cirrus", "Snow / Ice")
      levels(rm) <- cat

      # Mask pixels (no-data, cloud shadows, clouds medium, clouds high, cirrus)
      rf <- mask(rf, rm, maskvalue = c(0, 3, 8, 9, 10), updatevalue = NA)

      # Clip to our areas of interest
      dat <- crop(rf, project(vect(aoi), crs(rf)), snap = "out")

      # Return the file
      return(dat)
    })

    # Stitch them together
    stitch <- sprc(stitch)
    dat <- mosaic(stitch, fun = "median")

    # Rename the bands (check here:
    # https://sen2r.ranghetti.info/articles/outstructure#:~:text=Band%201%20%E2%80%93%20Aerosol%20(443%20nm,4%20%E2%80%93%20Red%20(665%20nm))
    names(dat) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B11", "B12")

    # Compute indices of interest
    ndvi <- nd(dat, "B8", "B4")
    ndwi <- nd(dat, "B3", "B8")
    ndmi <- nd(dat, "B8", "B11")
    ndsi <- nd(dat, "B3", "B11")
    best <- nd(dat, "B3", "B12")

    # Put all together into a single stack
    all <- c(dat, ndvi, ndwi, ndmi, ndsi, best)
    names(all) <- c(names(dat), "ndvi", "ndwi", "ndmi", "ndsi", "best")

    # Store the resulting raster to file
    writeRaster(all, filename, overwrite = T)
  }
})

# Once all files have been curated, remove the big folders
unlink(outdir_l1c, recursive = T)
unlink(outdir_l2a, recursive = T)

################################################################################
#### Training the Classifier
################################################################################
# Convert training classes to terra vect
class <- vect(classes)
class$ID <- 1:nrow(class)

# Design matrix that we want to go through
design <- expand_grid(
    Dates     = unique(class$Date)
  , Satellite = c("LANDSAT", "SENTINEL")
)

# Go through the design and extract reflectances
design$extracted <- lapply(1:nrow(design), function(x) {

  # Extract information of the specific iteration
  date <- design$Dates[x]
  sate <- design$Satellite[x]

  # Subset to training data of that specific date
  class_sub <- class[class$Date == date]

  # We need to assign new IDs because during the extraction, only the index of
  # the class will be returned.
  class_sub$ID <- 1:nrow(class_sub)

  # Load the satellite image that matches the date of the training data
  img <- rast(paste0("03_Data/01_RawData/", sate, "/", date, ".tif"))

  # Extract reflectances below the training classes
  extracted <- terra::extract(img, project(class_sub, crs(img)))

  # Remove any NAs
  extracted <- na.omit(extracted)

  # Join by id to add back the information about the date and land cover class
  # Order columns a bit nicer and return the result
  extracted <- left_join(extracted, values(class_sub), by = "ID")
  extracted <- dplyr::select(extracted, ID, Date, Class, everything())
  return(extracted)

})

# Keep 1000 sample points per category
design$Samples <- lapply(design$extracted, function(x) {
  samples <- x %>%
    group_by(Date, Class) %>%
    sample_n(size = 1000, replace = T) %>%
    ungroup()
  return(samples)
})

# Split dataset by satellite
samp_land <- design %>%
  subset(Satellite == "LANDSAT") %>%
  dplyr::select(c(Samples)) %>%
  unnest(Samples)
samp_sent <- design %>%
  subset(Satellite == "SENTINEL") %>%
  dplyr::select(c(Samples)) %>%
  unnest(Samples)

# Remove ID and Date, we won't need them anymore
samp_land <- select(ungroup(samp_land), -c(ID, Date))
samp_sent <- select(ungroup(samp_sent), -c(ID, Date))

# Create histograms of reflectances
samp_land %>%
  pivot_longer(B1:best, names_to = "Band", values_to = "Reflectance") %>%
  ggplot(aes(x = Reflectance, y = Class, col = Class, fill = Class)) +
    geom_density_ridges() +
    facet_wrap("Band", scales = "free") +
    theme_minimal() +
    theme(
        axis.text.x = element_blank()
      , axis.text.y = element_blank()
    ) +
    scale_fill_viridis_d(alpha = 0.5) +
    scale_color_viridis_d() +
    ggtitle("Landsat Reflectance Profiles")
samp_sent %>%
  pivot_longer(B1:best, names_to = "Band", values_to = "Reflectance") %>%
  ggplot(aes(x = Reflectance, y = Class, col = Class, fill = Class)) +
    geom_density_ridges() +
    facet_wrap("Band", scales = "free") +
    theme_minimal() +
    theme(
        axis.text.x = element_blank()
      , axis.text.y = element_blank()
    ) +
    scale_fill_viridis_d(alpha = 0.5) +
    scale_color_viridis_d() +
    ggtitle("Sentinel Reflectance Profiles")

# Train cart classifiers
mod_cart_land <- rpart(as.factor(Class) ~ ., data = samp_land, method = "class")
mod_cart_sent <- rpart(as.factor(Class) ~ ., data = samp_sent, method = "class")

# Visualize decision trees
par(mfrow = c(2, 1))
prp(mod_cart_land, main = "Landsat Classification Tree")
prp(mod_cart_sent, main = "Sentinel Classification Tree")

# Train random forest classifiers
mod_rand_land <- randomForest(as.factor(Class) ~ ., data = samp_land)
mod_rand_sent <- randomForest(as.factor(Class) ~ ., data = samp_sent)

# Check out the results
print(mod_rand_land)
print(mod_rand_sent)

################################################################################
#### Valiation
################################################################################
# Prepare a tibble that helps us to keep track of what we're doing
validation <- tibble(
    Satellite = c("Landsat", "Landsat", "Sentinel", "Sentinel")
  , Data      = list(samp_land, samp_land, samp_sent, samp_sent)
  , Model     = c("Cart", "RandomForest", "Cart", "RandomForest")
  , ModelObject = list(mod_cart_land, mod_rand_land, mod_cart_sent, mod_rand_sent)
)

# Compute variable importance
validation$Varimp <- lapply(validation$ModelObject, function(x) {
  importance <- varImp(x)
  importance <- as.data.frame(importance)
  importance$Band <- rownames(importance)
  return(importance)
})

# Go through the different configurations and run the validation
validation$Validation <- lapply(1:nrow(validation), function(x) {

  # Make bookkeeping a bit easier
  dat <- validation$Data[[x]]
  mod <- validation$Model[x]

  # Group the data and run k-fold cross validation
  k <- 5
  group <- sample(rep(1:k, each = round(nrow(dat) / k)))[1:nrow(dat)]
  valid <- lapply(1:k, function(y) {

    # Split testing and training data
    dat_trai <- dat[group != y, ]
    dat_test <- dat[group == y, ]

    # Fit model
    if (mod == "Cart") {
        model <- rpart(as.factor(Class) ~ ., data = dat_trai)
      } else {
        model <- randomForest(as.factor(Class) ~ ., data = dat_trai)
    }

    # Make prediction
    pred <- predict(model, dat_test, type = "class")

    # Put all into a dataframe
    result <- data.frame(
        Truth = as.factor(dat_test$Class)
      , Predi = pred
    )

    # Return it
    return(result)
  }) %>% do.call(rbind, .)
  return(valid)
})

# Function for accuracy assessment
accuracy <- function(x, y) {
  conf <- table(x, y, dnn = c("observed", "predicted"))
  specificity <- conf[1, 1] / sum(conf[1, 1], conf[1, 2])
  sensitivity <- conf[2, 2] / sum(conf[2, 1], conf[2, 2])
  accuracy <- sum(conf[1, 1], conf[2, 2]) / sum(conf[1, ], conf[2, ])
  res <- list(
      Confusion   = conf
    , Specificity = specificity
    , Sensitivity = sensitivity
    , Accuracy    = accuracy
  )
  return(res)
}

# Run the function on the predictions
validation$Confmat <- lapply(validation$Validation, function(x) {
  accuracy(x$Truth, x$Predi)
})

# Extract the information and put it into columns
validation <- mutate(validation
  , Confusion   = lapply(Confmat, function(x) {x$Confusion})
  , Specificity = sapply(Confmat, function(x) {x$Specificity})
  , Sensitivity = sapply(Confmat, function(x) {x$Sensitivity})
  , Accuracy    = sapply(Confmat, function(x) {x$Accuracy})
)
validation$Confmat <- NULL

# Check out the confusion matrices
validation$Confusion

# Convert them to dataframes
validation$Confusion <- lapply(validation$Confusion, as.data.frame)

# Store the validation object to file
write_rds(validation, "03_Data/03_Results/99_PanMapping.rds")

# Also store the pan mapping classes to file
writeVector(vect(classes)
  , "03_Data/02_CleanData/00_General_TrainingClasses.shp"
  , overwrite = T
)

################################################################################
#### CONTINUE HERE!!!
################################################################################
# Make a prediction
pred_land <- predict(land, mod_rand_land, na.rm = T)
pred_sent <- predict(sent, mod_rand_sent, na.rm = T)

# # Keep only best classes
# pred_land <- which.max(pred_land)
# pred_sent <- which.max(pred_sent)

# Remove dryland
pred_land <- subst(pred_land, 1, NA)
pred_sent <- subst(pred_sent, 1, NA)

# Store the predictions
writeRaster(pred_land, "Prediction_Landsat.tif", overwrite = T)
writeRaster(pred_sent, "Prediction_Sentinel.tif", overwrite = T)

################################################################################
#### Visualizations
################################################################################
