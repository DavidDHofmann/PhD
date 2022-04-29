################################################################################
#### Preparation of ORI Floodmaps
################################################################################
# Description: Preparation of the ORI floodmaps, making sure they match our
# reference raster. These floodmaps have been sent to me by ORI. We will then
# complement those maps by classifying our own floodmaps from MODIS satellite
# imagery.

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# load packages
library(tidyverse)    # For data wrangling
library(terra)        # For fast raster manipulation
library(lubridate)    # To handle dates
library(floodmapr)    # For floodmapping
library(pbmcapply)    # For multicore use

# Specify username and password to Earth Data account
username <- "DoDx9"
password <- "EarthData99"

################################################################################
#### Reclassify and Resample ORI Floodmaps
################################################################################
# Reclassify and resample ORI floodmaps
files <- dir(
    path        = "03_Data/01_RawData/ORI/01_Maps"
  , pattern     = "*.tif$"
  , full.names  = T
)

# Define filenames for resampled files
dir.create("03_Data/02_CleanData/00_Floodmaps", showWarnings = F)
dir.create("03_Data/02_CleanData/00_Floodmaps/01_Original", showWarnings = F)
newfiles <- paste0("03_Data/02_CleanData/00_Floodmaps/01_Original/", basename(files))

# Prepare reclassification table
rcl <- data.frame(old = c(0, 127, 255), new = c(1, 2, 0))

# Load reference floodmap
r <- rast("03_Data/01_RawData/ORI/FloodmapReference.tif")

# Loop through the files, resample them to the reference raster, reclassify
# values, and store them. Note that I will trim the raster to save space.
files <- files[!file.exists(newfiles)]
newfiles <- newfiles[!file.exists(newfiles)]
if (length(files) > 0) {
  cat("Preparing all layers provided by ORI...\n")
  for (i in 1:length(files)) {
    ori <- rast(files[i])[[1]]
    ori <- resample(ori, r, "near")
    ori <- classify(ori, rcl)
    names(ori) <- ymd(basename(files[[i]]))
    writeRaster(ori, newfiles[i], overwrite = T)
    cat(i, "out of", length(files), "done...\n")
  }
}

# Remove the "aux.xml" files that were created
file.remove(dir(
    path        = "03_Data/02_CleanData/00_Floodmaps"
  , pattern     = ".aux.xml$"
  , full.names  = T
))

################################################################################
#### Determine Missing Floodmaps
################################################################################
# Identify all ORI files
files <- dir(
    path        = "03_Data/01_RawData/ORI/01_Maps"
  , pattern     = "*.tif$"
  , full.names  = T
)

# Extract floodmap date from filenames
ori_dates <- files %>%
  basename() %>%
  ymd()

# Given the MODIS satellite schedule, floodmaps could be generated every 8-th
# day, starting from day one of each year. Generate a vector that shows those
# dates.
pot_dates <- lapply(2000:2021, function(x) {
  seq(as.Date(paste0(x, "-01-01")), as.Date(paste0(x, "-12-31")), by = "8 days")
}) %>% do.call(c, .)

# Make sure that all floodmap dates are actually contained within those
# potential dates
all(ori_dates %in% pot_dates)

# Check for which of those dates we actually have a floodmap. Those dates for
# which we miss a floodmap, we could try to get one ourselves
pot_dates <- data.frame(
    Date             = pot_dates
  , Floodmap         = pot_dates %in% ori_dates
  , stringsAsFactors = F
)

# Let's try to better understand for what dates we are missing floodmaps
table(pot_dates$Floodmap, year(pot_dates$Date))
table(pot_dates$Floodmap, month(pot_dates$Date))

# We can make two observations. First, we miss data since end of 2018, which
# makes sense, because ORI has only provided me with data until then. Second, we
# find that many maps are missing in the rainy season, which again makes sense,
# as bimodality of the MODIS satellite maps is often not given during those
# months. We can actually plot the latter quite nicely.
pot_dates %>% mutate(Month = month(Date)) %>% count(Month, Floodmap) %>% subset(Floodmap) %>%
  ggplot(aes(x = Month, y = n)) +
    geom_line() +
    geom_point() +
    theme_minimal()

# Load dispersal data and identify first and last date. We can use those to
# determine for timerange we should try to download additional floodmaps
dis_dates <- "03_Data/02_CleanData/00_General_Dispersers.csv" %>%
  read_csv(show_col_types = F) %>%
  subset(State == "Disperser") %>%
  pull(Timestamp) %>%
  as.Date() %>%
  range()

# Subset to those dates (and expand the range slightly)
todownload <- subset(pot_dates
  , Date >= dis_dates[1] - days(8)
  & Date <= dis_dates[2] + days(8)
)

# Obviously we don't need to download data for those dates for which we already
# have a floodmap
todownload <- subset(todownload, !Floodmap) %>% pull(Date)

# Double check if the files have already been downloaded
exists <- file.exists(paste0("03_Data/01_RawData/MODIS/MCD43A4/", todownload, ".tif"))
todownload <- todownload[!exists]

# How many do we still need to download?
length(todownload)

################################################################################
#### Download Modis Maps
################################################################################
# Download the maps
downloaded <- c()
if (length(todownload) > 0) {
    for (i in 1:length(todownload)){

      # Download i-th file
      tryCatch(modis_download(
          dates     = todownload[i]
        , outdir    = "03_Data/01_RawData/MODIS/MCD43A4"
        , tmpdir    = tempdir()
        , username  = username
        , password  = password
        , overwrite = F
      ), error = function(e) {cat("There was an error with", i, "\n")})

      # Print update
      cat(i, "out of", length(todownload), "downloaded\n")
    }
  } else {
    cat("No more MODIS images need to be downloaded...\n")
}

################################################################################
#### Classify Them
################################################################################
# We will run the classification algorithm twice. First, using the static
# polygon for watermaps (which is the default in the floodmapr package).
# Afterwards, we'll run the algorithm again, this time using dynamic polygons
# for watermaps.
for (i in 1:2) {

  # Specify if watermask should be dynamic or not
  dynamic <- ifelse(i == 1, F, T)
  if (dynamic) {
      cat("Classifying floodmaps from MODIS images using dynamic water polygons...\n")
    } else {
      cat("Classifying floodmaps from MODIS images using static water polygons...\n")
  }

  # Load reference raster (this time we'll use another floodmap as reference)
  r <- rast("03_Data/01_RawData/ORI/FloodmapReference.tif")

  # Retrieve downloaded files
  downloaded <- dir(
      path        = "03_Data/01_RawData/MODIS/MCD43A4"
    , pattern     = ".tif$"
    , full.names  = T
  )

  # Also check all already classified maps
  files <- dir(
      path        = "03_Data/02_CleanData/00_Floodmaps/01_Original"
    , pattern     = ".tif$"
    , full.names  = T
  )

  # Subset to files that have not been classified yet
  downloaded <- downloaded[!basename(downloaded) %in% basename(files)]

  # Retrieve dates from filenames
  dates <- downloaded %>%
    basename() %>%
    ymd()

  # Loop through each of the downloaded images and classify them
  if (length(downloaded) > 0) {
      for (i in 1:length(downloaded)) {

        # Load respective file
        loaded <- modis_load(downloaded[i])

        # Identify filenames of all already classified watermaps
        filenames <- dir(
            path        = "03_Data/02_CleanData/00_Floodmaps/01_Original"
          , pattern     = ".tif$"
          , full.names  = T
        )

        # Identify dates of all already classified watermaps
        filedates <- filenames %>%
          basename() %>%
          ymd()

        # Create a dynamic watermask
        if (dynamic) {
          watermask <- modis_watermask(
              date      = dates[i]
            , filenames = filenames
            , filedates = filedates
            , years     = 5
            , threshold = 0.975
          )
        }

        # Check for bimodality in the modis image
        if (dynamic) {
          is_bimodal <- modis_bimodal(
              x         = loaded
            , watermask = watermask
            , drymask   = NULL
          )
        } else {
          is_bimodal <- modis_bimodal(
              x         = loaded
            , watermask = NULL
            , drymask   = NULL
          )
        }

        # In case the file is not bimodal, we can't classify it and need to skip
        if (!is_bimodal){
          cat(paste0("Image ", dates[i], " is not bimodal. Can't classify watermap. Skipping to next image...\n"))
          next
        } else {
          cat(paste0("Image ", dates[i], " is bimodal. Will be classified now...\n"))
        }

        # Classify modis image. We provide a dynamic watermask here
        if (dynamic) {
            classified <- modis_classify(
                x                 = loaded
              , watermask         = watermask
              , drymask           = NULL
              , ignore.bimodality = T
            )
          } else {
            classified <- modis_classify(
                x                 = loaded
              , watermask         = NULL
              , drymask           = NULL
              , ignore.bimodality = T
            )
        }

        # Print update
        cat(paste0("Image ", dates[i], " Classified. Resampling and storing now...\n"))

        # We also want to resample the classified image, so that it matches the
        # reference raster. We'll use the nearest neighbor method here
        classified <- terra::resample(
            x       = classified
          , y       = r
          , method  = "near"
        )

        # Specify final filename
        name <- dates[i]
        name <- paste0("03_Data/02_CleanData/00_Floodmaps/01_Original/", name, ".tif")
        names(classified) <- dates[i]

        # Store raster
        writeRaster(classified, name, overwrite = T)

        # Print update
        cat("Image", i, "out of", length(downloaded), "finished...\n")
      }
    } else {
    cat("No more maps to classify...\n")
  }
}

# Remove any .aux.xml files
file.remove(
    dir(
      path       = "03_Data/02_CleanData/00_Floodmaps/01_Original"
    , pattern    = ".aux.xml$"
    , full.names = T
  )
)

################################################################################
#### Resample Maps to Study Area
################################################################################
# Identify all floodmaps
files <- dir(
    path       = "03_Data/02_CleanData/00_Floodmaps/01_Original"
  , pattern    = ".tif$"
  , full.names = T
)

# How many are there?
length(files)

# Only keep the ones that are not resampled yet
files <- files[!basename(files) %in% dir(path = "03_Data/02_CleanData/00_Floodmaps/02_Resampled")]

# Prepare new filenames
newnames <- paste0("03_Data/02_CleanData/00_Floodmaps/02_Resampled/", basename(files))

# Resample floodmaps
dir.create("03_Data/02_CleanData/00_Floodmaps/02_Resampled", showWarnings = F)
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")
if (length(files) > 0) {
  cat("Resampling floodmaps...\n")
  invisible(pbmclapply(1:length(files), ignore.interactive = T, mc.cores = 1, function(x) {
    flood <- rast(files[x])
    flood <- resample(flood, r, method = "near")
    flood <- trim(flood)
    writeRaster(flood, filename = newnames[x])
  }))
}

# ################################################################################
# #### Visualizations
# ################################################################################
# # Let's visualize a random subset of the classified floodmaps
# sample <- sample(newnames, 6)
# flood <- rast(sample)
# names(flood) <- ymd(basename(sample))
#
# # Plot
# flood %>%
#   as.data.frame(xy = T) %>%
#   pivot_longer(cols = 3:8, names_to = "Date", values_to = "Value") %>%
#   mutate(Value = as.factor(Value)) %>%
#   mutate(Value = fct_recode(Value, "Dryland" = "0", "Water" = "1", "Clouds" = "2")) %>%
#   ggplot(aes(x = x, y = y, fill = as.factor(Value))) +
#     geom_raster() +
#     facet_wrap("Date", labeller = labeller(Date = function(x){ymd(substr(x, 2, 11))})) +
#     theme_minimal() +
#     coord_sf() +
#     scale_fill_manual(values = c("white", "cornflowerblue", "gray"), name = "Category")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/04_Floodmaps_SessionInfo.rds")
cat("Done :)\n")
