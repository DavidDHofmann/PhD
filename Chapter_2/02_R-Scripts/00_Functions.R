################################################################################
#### Function to Determine the Corrected Sentinel 2 Names
################################################################################
# Function to determine the name of the corrected Sentinel 2 product
correctedName <- function(x) {
  corr <- gsub(basename(x), pattern = "MSIL1C", replacement = "MSIL2A")
  return(corr)
}

################################################################################
#### Function to calculate Distances on a Raster Efficiently
################################################################################
# Compute distance to a certain value on a raster
distanceTo <- function(x, value = 1) {

  # Convert to regular raster if terra raster is provided
  if (inherits(x, "SpatRaster")) {
    x <- raster(x)
    returnterra <- T
  }

  # Function requires two packages
  require(spatstat)
  require(maptools)

  # Convert raster to ppp object
  ppp <- x %>%
    rasterToPoints(., fun = function(x){x == value}, spatial = TRUE) %>%
    spTransform(., CRS("+init=epsg:32734")) %>%
    as(., "ppp")

  # Create empty raster onto which the distances are calculated
  distance <- raster(x)
  values(distance) <- distance %>%
    as(., "SpatialPoints") %>%
    spTransform(., CRS("+init=epsg:32734")) %>%
    as(., "ppp") %>%
    nncross(., ppp) %>%
    .[["dist"]]

  # Return the distance raster
  return(distance)
}

################################################################################
#### Function to Resample GPS Data to a Desired Resolution
################################################################################
# Function to resample fixes to a desired resolution
.resFix <- function(data, hours, start, tol = 0.5) {

  # Identify the first date at which a fix was taken and set the time to the
  # desired start time
  first <- range(data$Timestamp)[1] %>%
    update(., hour = start, min = 0, sec = 0)

  # Identify the last date at which a fix was taken and set the end time to 24
  last <- range(data$Timestamp)[2] %>%
    update(., hour = 24, min = 0, sec = 0)

  # Prepare a range of dates for which we would expect to find data according
  # to the specified sampling scheme
  dates <- seq(first, last, by = paste0(hours, " hours")) %>%
    as.data.frame() %>%
    set_names("Timestamp")

  # For each Timestamp we now identify the closest fix
  closest <- sapply(1:nrow(dates), function(x){
    index <- which.min(abs(dates$Timestamp[x] - data$Timestamp))[1]
    close <- as.numeric(abs(dates$Timestamp[x] - data$Timestamp[index]), units = "hours") <= tol

    # In case the fix is close enough, return its index
    if (close) {
      return(index)
    } else {
      return(NA)
    }
  })

  # Return respective fixes
  closest <- na.omit(closest)
  return(data[closest, ])
}

# Wrapper to run the function on multiple individuals and multiple cores
resampleFixes <- function(data, ID = NULL, hours, start, tol = 0.5, cores = 1) {
  if (is.null(ID)) {
    data_resampled <- .resFix(data, hours, start, tol)
  } else {
    data_resampled <- data %>% nest(data = -ID)
    data_resampled$data <- pbmclapply(1:nrow(data_resampled), ignore.interactive = T, mc.cores = cores, function(x) {
      .resFix(data_resampled$data[[x]], hours = hours, start = start, tol = tol)
    })
    data_resampled <- unnest(data_resampled, data)
  }
  return(data_resampled)
}

################################################################################
#### Function to Compute Step Metrics
################################################################################
# Function to compute the absolute turning angle
absTA <- function(dx, dy) {
  absta <- atan2(dy, dx)
  absta <- (absta - pi / 2) * (-1)
  absta <- ifelse(absta < 0, 2 * pi + absta, absta)
  return(absta)
}

# Function to compute the relative turning angle
relTA <- function(absta) {
  relta <- (absta[-1] - absta[-length(absta)])
  relta <- c(NA, relta)
  relta <- ifelse(relta > +pi, relta - 2 * pi, relta)
  relta <- ifelse(relta < -pi, 2 * pi + relta, relta)
  return(relta)
}

# Function to compute step metrics
stepMet <- function(x, y, timestamp = NULL) {

    # Compute distances moved in x and y direction
    dx <- c(x[-1], NA) - x
    dy <- c(y[-1], NA) - y

    # Calculate step length
    sl <- sqrt(dx ** 2 + dy ** 2)

    # Compute absolute turn angle
    absta <- absTA(dx, dy)

    # Compute relative turn angle
    relta <- relTA(absta)

    # Compute step duration
    if (!is.null(timestamp)) {
      dt <- difftime(lead(timestamp), timestamp, units = "hours")
    }

    # Put metrics into data.frame
    metrics <- data.frame(sl = sl, absta = absta, relta = relta)
    if (!is.null(timestamp)) {
      metrics$dt <- dt
    }

    # Return the metrics
    return(metrics)
}
