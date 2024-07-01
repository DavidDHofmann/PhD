################################################################################
#### Covariate Extraction
################################################################################
# Description: Use the step selection data generated in the previous script and
# extract covariate values below each observed and random step. To speed up
# extractions, we will extract covariate values along interpolated points
# instead of along lines. For dynamic covariates, we'll extract data from the
# layer that is closest in date to the actual step.

# Note: I've been benchmarking multiple extraction techniques, including using
# the velox, raster, or terra packages. The speed of the different packages
# varies greatly depending on whether we are extracting along lines or points,
# and depending on whether the data is all loaded into memory or not. It appears
# that extracting along interpolated points using the raster package yields the
# best performance, yet only so if the data is all loaded into memory. The
# benefit of using raster or terra over using velox is also that in case of many
# covariates, we can easily extract only from the indexed layer instead of
# having to extract from all layers and only index later.

# Clear R's brain
rm(list = ls())

# Load packages
library(tidyverse)  # To wrangle data
library(terra)      # To handle spatial data
library(raster)     # To handle spatial data
library(lubridate)  # To handle dates
library(pbmcapply)  # To run stuff on multiple cores
library(hms)        # To work with times

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

# Specify name of the file for the extracted data. If it already exists, load
# it, otherwise, load the raw data
filename <- "03_Data/02_CleanData/SSFExtracted.rds"
if (file.exists(filename)) {
    ssf <- read_rds(filename)
  } else {
    ssf <- read_csv("03_Data/02_CleanData/SSF.csv", show_col_types = F)

    # Subset to some of the steps for now
    # keep <- sample(unique(ssf$step_id), size = 100)
    # ssf <- subset(ssf, step_id %in% keep)

    # Let's generate a set of interpolated points along each of the steps
    ssf$Points <- pbmclapply(
        X                  = 1:nrow(ssf)
      , ignore.interactive = T
      , FUN                = function(i) {
      p <- interpolatePoints(
          x1 = ssf$x[i]
        , x2 = ssf$x_to[i]
        , y1 = ssf$y[i]
        , y2 = ssf$y_to[i]
        , by = 250 / 111000
      )
      return(p)
    })
}

################################################################################
#### Seasons
################################################################################
# Identify the season into which each of the fixes fall
if (!all(c("SeasonClimate", "SeasonHerbivores") %in% names(ssf))) {
  ssf$SeasonClimate    <- getRainySeason(ssf$Timestamp)
  ssf$SeasonHerbivores <- getHerbivoreSeason(ssf$Timestamp)
}

################################################################################
#### Covariates
################################################################################
# Load covariate overview and the lookup table, then combine them
covs <- read_rds("03_Data/02_CleanData/Covariates.rds") %>% dplyr::select(-Dates)
look <- read_rds("03_Data/02_CleanData/LookupTable.rds")

# Check them
print(covs, n = 29)
head(look)

# Ignore aggregated data
covs <- subset(covs, Type != "DynamicAggregated")
look <- subset(look, Type != "DynamicAggregated")

# Combine them. This will give us a table with all covariates for both static
# and dynamic environments, as well as a lookup table for each covariate that
# tells us from which layer to extract for a given timestamp
covs <- look %>%
  nest(Lookup = -c(Covariate, Type)) %>%
  left_join(covs, ., by = c("Covariate", "Type"))

# We created the lookup table for way more dates than we actually have
# dispersers. We will therefore simplify it later. However, we need to know the
# timestamps during dispersal for this
dispdates <- ssf %>%
  pull(TimestampRounded) %>%
  unique()

# Now we loop through the different covariates (both static and dynamic ones)
# and extract values along each observed and random step
for (i in seq_len(nrow(covs))) {

  # We'll add the moonlight statistics later, as it's not in raster format and
  # can easily be joined by the timestamp
  if (covs$Covariate[i] == "Moonlight") {
    next
  }

  # Skip covariates that are already extracted (this helps when the extraction
  # should be rerun only for specific covariates)
  cov <- paste0(covs$Covariate[i], covs$Type[i])
  if (cov %in% names(ssf)) next
  cat("Extracting data from the following layer:", cov, "\n")

  # Load the covariate layer(s)
  covariate <- stack(covs$Filename[i])

  # We can now extract the lookup table for this specific covariate
  covs_sub <- covs %>%
    subset(Covariate == covs$Covariate[i] & Type == covs$Type[i]) %>%
    unnest(Lookup) %>%
    dplyr::select(-c(Type, Covariate, Filename))

  # Presently, the lookup table references a layer by its index. However, we
  # will later remove layers that represent dates that we actually don't need
  # (for performance purposes). Hence, the layerindex might change and we should
  # rather index layers by their names
  covs_sub$Layername <- names(covariate)[covs_sub$Layerindex]

  # Now we drop any layer that corresponds to a date which does not show up in
  # the dispersal dates. This allows us to reduce the number of layers we need
  # to carry with us
  covs_sub <- subset(covs_sub, Timestamp %in% dispdates)

  # If we did everything correctly, we should still have a layer for every
  # dispersal timestamp. Let's verify this.
  allgood <- all(ssf$TimestampRounded %in% covs_sub$Timestamp)
  if (!allgood) {
    stop("We are missing layers of a few timestamps\n")
  }

  # Drop all other layers from the covariate stack and read the remaining data
  # into memory
  covariate <- covariate[[sort(unique(covs_sub$Layerindex))]]
  covariate <- readAll(covariate)

  # Finally, we join the lookup table with the ssf data
  ssf <- covs_sub %>% left_join(ssf, ., by = c("TimestampRounded" = "Timestamp"))

  # Extract data
  extracted <- pbmclapply(
      X                  = 1:nrow(ssf)
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(i) {
    if (is.na(ssf$Layername[i])) {
      return(NA)
    }
    extr  <- raster::extract(covariate[[ssf$Layername[i]]], ssf$Points[[i]])
    extr  <- mean(extr, na.rm = T)
    return(extr)
  }) %>% do.call(c, .)

  # Remove undesired columns
  ssf$Layername  <- NULL
  ssf$Layerdate  <- NULL
  ssf$Layerindex <- NULL

  # Combine with ssf data
  ssf                   <- cbind(ssf, extracted)
  names(ssf)[ncol(ssf)] <- cov

  # Store the extracted info to file
  write_rds(ssf, file = filename)

  # Make space
  rm(covariate)
  invisible(gc())

}

# Remove moonlight stats if desired
print(names(ssf))
ssf$meanMoonPhase <- NULL
ssf$meanMoonAlt   <- NULL
ssf$meanSunAlt    <- NULL
ssf$meanMoonlight <- NULL
ssf$Night         <- NULL
ssf$NightPercent  <- NULL
ssf$LightType     <- NULL

# Add moonlight statistics if needed
if (!"meanMoonlight" %in% names(ssf)) {
  moon <- read_rds("03_Data/02_CleanData/Moonlight.rds")
  names(moon)[1] <- "TimestampRounded"
  ssf <- left_join(ssf, moon, by = "TimestampRounded")
}

# Store extracted data to file
# ssf <- dplyr::select(ssf, -Points)
write_rds(ssf, "03_Data/02_CleanData/SSFExtracted.rds")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/02_Analysis/01_CovariateExtraction.rds")
cat("Done :)\n")
