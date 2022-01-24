################################################################################
#### Sentinel 2 Data Download for Pan Mapping
################################################################################
# Description: Here, we'll download Sentinel 2 imagery and train a classifier
# that is intended to detect wet pans. We will later use the trained classifier
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
wd <- "/home/david/Schreibtisch/PanMapping2"
setwd(wd)

# Create necessary directories
dir.create("Landsat", showWarnings = F)
dir.create("Sentinel", showWarnings = F)
dir.create("Sentinel/L1C", showWarnings = F)
dir.create("Sentinel/L2A", showWarnings = F)
outdir_l1c <- "Sentinel/L1C"
outdir_l2a <- "Sentinel/L2A"

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
    path       = "TrainingClasses"
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
  coord_sf() +
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

# Define the area of interest
aoi <- ee$Geometry$Polygon(
  list(
      c(extent(classes)[1], extent(classes)[3])
    , c(extent(classes)[1], extent(classes)[4])
    , c(extent(classes)[2], extent(classes)[4])
    , c(extent(classes)[2], extent(classes)[3])
  )
)

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

  # Specify filename
  filename <- paste0(wd, "/Landsat/", dates[i], ".tif")

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
    ee_print(collection)

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

# Define a polygon that indicates the are for which we need to download sentinel
# 2 data
aoi <- extent(classes)
aoi <- as(aoi, "SpatialPolygons")
crs(aoi) <- "+init=epsg:4326"

# Check the overlap with sentinel 2 tiles
tile <- tiles_intersects(st_as_sf(aoi), out_format = "sf")
tile <- as(tile, "Spatial")

# Visualize it
plot(tile)
plot(aoi, add = T, border = "red", lty = 2)
plot(classes, add = T)

# Go through the dates and download the sentinel data
for (i in 1:length(dates)) {

  # Specify the final filename
  filename <- paste0(wd, "/Sentinel/", dates[i], ".tif")

  # If the file does not exist yet, download and process it
  if (!file.exists(filename)) {

    # Find products for the range of dates
    available <- s2_list(
        spatial_extent = st_as_sf(aoi)
      , time_interval  = c(dates[i] - days(10), dates[i] + days(10))
      , level          = "auto"
    )

    # Check the metadata and create a regular date column
    meta <- as.data.frame(available)
    meta$Date <- as.Date(meta$sensing_datetime)

    # Find the image(s) that is/are closest in date to the input date
    index <- which(abs(meta$Date - dates[i]) == min(abs(meta$Date - dates[i])))

    # Order the product(s)
    s2_order(available[index])

    # Depending on the product level, specify the correct output directory
    outdir <- ifelse(meta$level[index] == "1C", outdir_l1c, outdir_l2a)

    # Download the specific product
    success <- tryCatch({
      s2_download(available[index], outdir = outdir)
    }, error = function(e) {return(e)})

    # Process the downloaded products and prepare L2C
    files <- dir(path = outdir_l1c, include.dirs = T, full.names = T)
    sen2cor(files
      , outdir    = outdir_l2a
      , parallel  = detectCores() - 1
      , overwrite = F
    )

    # Create a VRT
    files <- dir(path = outdir_l2a, include.dirs = T, full.names = T)
    vrt <- s2_translate(files, format = "VRT", outdir = tempdir())

    # Load it and clip it to our area of interest
    dat <- rast(vrt)
    dat <- crop(dat, project(vect(aoi), crs(dat)), snap = "out")

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
    writeRaster(all, paste0("Sentinel/", dates, ".tif"), overwrite = T)
  }
}


################################################################################
#### Training the Classifier
################################################################################
# Load landsat 8 data
land <- rast("/home/david/Schreibtisch/PanMapping2/Landsat/2018-08-18.tif")
sent <- rast("/home/david/Schreibtisch/PanMapping2/Sentinel/2018-08-18.tif")

# Convert training classes to terra vect
class <- vect(classes)
class$ID <- 1:nrow(class)

# Extract reflectances below the random points
extracted_land <- terra::extract(land, project(class, crs(land)))
extracted_sent <- terra::extract(sent, project(class, crs(sent)))

# Check for NAs
apply(extracted_land[, 2:15], 2, function(x) {
  sum(is.na(x))
})
apply(extracted_sent[, 2:17], 2, function(x) {
  sum(is.na(x))
})

# Remove them
extracted_sent <- na.omit(extracted_sent)
extracted_land <- na.omit(extracted_land)

# Join by id to add back the information about the date and land cover class
extracted_land <- left_join(extracted_land, values(class), by = "ID")
extracted_sent <- left_join(extracted_sent, values(class), by = "ID")

# Order columns a bit nicer
extracted_land <- dplyr::select(extracted_land, ID, Date, Class, everything())
extracted_sent <- dplyr::select(extracted_sent, ID, Date, Class, everything())

# For each class, keep only 1000 points
samp_land <- extracted_land %>% group_by(Date, Class) %>% sample_n(size = 1000, replace = T)
samp_sent <- extracted_sent %>% group_by(Date, Class) %>% sample_n(size = 1000, replace = T)

# Let's only keep two classes
# samp_land$Class <- ifelse(samp_land$Class %in% c("Wetpan", "Water"), "Water", "DryLand")
# samp_sent$Class <- ifelse(samp_sent$Class %in% c("Wetpan", "Water"), "Water", "DryLand")
table(samp_land$Class)
table(samp_sent$Class)

# Create histograms of reflectances
p1 <- samp_land %>%
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
p2 <- samp_sent %>%
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

# Remove ID and Date
samp_land <- select(ungroup(samp_land), -c(ID, Date))
samp_sent <- select(ungroup(samp_sent), -c(ID, Date))

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

# Assess and plot variable importance
p3 <- varImp(mod_rand_land) %>%
  as.data.frame() %>%
  mutate(Band = rownames(.)) %>%
  ggplot(aes(x = reorder(Band, Overall), y = Overall)) +
    geom_segment(aes(x = reorder(Band, Overall), xend = reorder(Band, Overall), y = 0, yend = Overall), col = "gray50") +
    geom_point(col = "orange") +
    theme_minimal() +
    coord_flip() +
    xlab("Importance") +
    ylab("Band / Index") +
    ggtitle("Landsat 8 Variable Importance")
p4 <- varImp(mod_rand_sent) %>%
  as.data.frame() %>%
  mutate(Band = rownames(.)) %>%
  ggplot(aes(x = reorder(Band, Overall), y = Overall)) +
    geom_segment(aes(x = reorder(Band, Overall), xend = reorder(Band, Overall), y = 0, yend = Overall), col = "gray50") +
    geom_point(col = "orange") +
    theme_minimal() +
    coord_flip() +
    xlab("Importance") +
    ylab("Band / Index") +
    ggtitle("Sentinel 2 Variable Importance")

# Put legend to bottom
p1 <- p1 + theme(legend.position = "bottom")
p2 <- p2 + theme(legend.position = "bottom")

# Put the plots together
p <- ggarrange(p1, p2, p3, p4, labels = "auto")
ggsave(
    plot     = p
  , filename = "/home/david/ownCloud/University/15. PhD/Chapter_2/04_Manuscript/99_Reflectances.png"
  , bg       = "white"
  , scale    = 1.2
)

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

# Visualize variable importance
p5 <- validation %>%
  select(Satellite, Model, Varimp) %>%
  unnest(Varimp) %>%
  ggplot(aes(y = reorder(Band, Overall), x = Overall, xmin = 0, xmax = Overall, colour = Model)) +
    geom_pointrangeh(position = position_dodgev(height = 0.5), fatten = 1.5) +
    theme_minimal() +
    xlab("Importance") +
    ylab("Band / Index") +
    facet_wrap("Satellite", scales = "free") +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("cornflowerblue", "orange"))

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

# Visualize them
p6 <- validation %>% select(Satellite, Model, Confusion) %>%
  unnest(Confusion) %>%
  ggplot(aes(x = observed, y = predicted, fill = Freq)) +
    geom_tile(colour = "black") +
    geom_text(aes(label = Freq), col = "black", size = 2) +
    facet_wrap(~ Satellite + Model, ncol = 4) +
    scale_fill_gradient(low = "white", high = adjustcolor("orange", alpha.f = 0.3)) +
    theme_minimal() +
    theme(legend.position = "none", panel.grid = element_blank()) +
    coord_equal() +
    xlab("Observed") +
    ylab("Predicted")
p7 <- validation %>%
  select(Satellite, Model, Specificity, Sensitivity, Accuracy) %>%
  pivot_longer(3:5, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(y = Metric, x = Value, xmin = 0, xmax = Value)) +
    geom_vline(xintercept = 1, lty = 2, col = "gray50", lwd = 0.5) +
    geom_pointrangeh(position = position_dodgev(height = 0.5), fatten = 1.5, col = "orange") +
    theme_minimal() +
    xlab("") +
    ylab("Metric") +
    facet_wrap(~ Satellite + Model, ncol = 4) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("cornflowerblue", "orange"))

# Put the plots together
p <- ggarrange(p6, p7, nrow = 2)

# Store them
ggsave(
    plot     = p
  , filename = "/home/david/ownCloud/University/15. PhD/Chapter_2/04_Manuscript/99_ClassificationValidation.png"
  , bg       = "white"
)

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
# Get satellite image of study area
sat <- read_osm(st_as_sf(study), type = "bin", zoom = 11)
sat <- as(sat, "Raster")
sat <- rast(sat)
sat <- project(sat, crs(class))
study <- project(study, crs(class))

# Prepare colors for plot
cols <- c("orange", "cornflowerblue")
cols <- cols[as.factor(class$Class)]

# Plot of study area
plotRGB(sat)
plot(study, add = T, border = "black", lwd = 4, col = NA)
plot(class, add = T, col = cols, border = cols)
text(class, "Class", pos = 3, halo = F, cex = 0.9, col = "white")

# Show reflectances below different polygons
dat %>%
  pivot_longer(2:ncol(.), names_to = "Band", values_to = "Reflectance") %>%
  group_by(Class, Band) %>%
  ggplot(aes(x = Class, y = Reflectance, col = Class)) +
    geom_point(alpha = 0.4) +
    scale_color_manual(values = c("orange", "cornflowerblue")) +
    facet_wrap("Band", scales = "free") +
    theme_minimal()

# Show classification tree
prp(cartmodel, box.palette = c("orange", "cornflowerblue"))

# Load predictions
pred1 <- rast("Sentinel/Prediction_2021-04-05.tif")
pred2 <- rast("Sentinel/Prediction_2021-10-12.tif")

# Visualize predictions
plotRGB(sat, main = "Prediction for 2021-04-05", mar = c(1, 1, 1, 1))
plot(pred1, col = "white", frame = F, axes = F, add = T)
plotRGB(sat, main = "Prediction for 2021-10-12", mar = c(1, 1, 1, 1))
plot(pred2, col = "white", frame = F, axes = F, add = T)
