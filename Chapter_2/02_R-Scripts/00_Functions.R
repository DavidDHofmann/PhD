################################################################################
#### Function to Determine the Corrected Sentinel 2 Names
################################################################################
# Function to determine the name of the corrected Sentinel 2 product
correctedName <- function(x) {
  corr <- gsub(basename(x), pattern = "MSIL1C", replacement = "MSIL2A")
  return(corr)
}

################################################################################
#### Function to Compute Normalized Difference Index
################################################################################
# Function to compute the normalized difference (nd) index of two bands
nd <- function(img, band_x, band_y) {
  x <- img[[band_x]]
  y <- img[[band_y]]
  nd <- (x - y) / (x + y)
  return(nd)
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
resampleFixes <- function(data, hours, start, tol = 0.5) {

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

################################################################################
#### Function to Compute Bursts
################################################################################
# Function to compute bursts
computeBursts <- function(data) {
  bursted <- data %>%
    mutate(dt_ = Timestamp - lag(Timestamp)) %>%
    mutate(dt_ = as.numeric(dt_, units = "hours")) %>%
    ungroup() %>%
    mutate(NewBurst = ifelse(
      dt_ > 8.25 |
      (dt_ > 4.25 & hour(TimestampRounded) != 15) |
      is.na(dt_), yes = 1, no = 0)
    ) %>%
    mutate(burst_id = cumsum(NewBurst)) %>%
    dplyr::select(-c(NewBurst, dt_))
  return(bursted)
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

# Function to compute step metrics (data needs to be projected)
stepMetrics <- function(data) {

    # Get relevant data
    x <- data$x
    y <- data$y
    timestamp <- data$Timestamp

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
    metrics <- data.frame(x_to = lead(x), y_to = lead(y), sl = sl, absta = absta, relta = relta)
    if (!is.null(timestamp)) {
      metrics$dt <- dt
    }

    # Combine with original data
    data <- cbind(data, metrics)
    data <- data[-nrow(data), ]

    # Return the metrics
    return(data)
}

################################################################################
#### Function to Generate Random Steps
################################################################################
# Function to generate random steps
randomSteps <- function(data, n_rsteps, scale, shape, slmax = NULL) {

  # Generate a new column that indicates that the steps are "observed" steps
  data$case <- 1

  # Generate step ids
  data$step_id <- 1:nrow(data)

  # Cannot work with steps that have no turning angle, so remove them
  data <- subset(data, !is.na(relta))

  # Create a new dataframe into which we can put alternative/random steps
  rand <- data[rep(1:nrow(data), each = n_rsteps), ]

  # Indicate that these steps are random steps (case = 0)
  rand$case <- 0

  # Sample random turning angles
  rand$sl <- rgamma(n = nrow(rand)
    , scale = scale
    , shape = shape
  )

  # If a maximum step length is provided, truncate the sampled step length if
  # necessary
  if (!is.null(slmax)) {
    rand$sl <- ifelse(rand$sl > slmax, slmax, rand$sl)
  }

  # Sample random step lengths
  rand$relta_new <- runif(n = nrow(rand)
    , min = -pi
    , max = +pi
  )

  # Calculate new "absolute" turning angle
  rand$relta_diff <- rand$relta_new - rand$relta
  rand$absta <- rand$absta + rand$relta_diff
  rand$absta <- ifelse(rand$absta < 0, 2 * pi + rand$absta, rand$absta)
  rand$absta <- ifelse(rand$absta > 2 * pi, rand$absta - 2 * pi, rand$absta)
  rand$relta <- rand$relta_new

  # Remove undesired stuff
  rand$relta_new <- NULL
  rand$relta_diff <- NULL

  # Put steps together
  all <- rbind(data, rand)
  all <- arrange(all, ID, step_id, desc(case))

  # Calculate new endpoints
  all$x_to <- all$x + sin(all$absta) * all$sl
  all$y_to <- all$y + cos(all$absta) * all$sl

  # Return the final dataframe
  return(all)

}

################################################################################
#### Function to Reproject Coordinates
################################################################################
# Reproject coordinates to another projection
reprojectCoords <- function(xy, from = NULL, to = NULL) {
  xy <- vect(xy, crs = from)
  xy <- terra::project(xy, to)
  return(crds(xy))
}

################################################################################
#### Function to Interpolate Between Points
################################################################################
# Function to interpolate coordinates between two points
interpolatePoints <- function(x1, x2, y1, y2, by = 1){

  # Calculate length of line between points
  length <- sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

  # Calculate how many segments we need
  nsegs <- max(ceiling(length / by), 1)

  # Interpolate between points
  x <- seq(x1, x2, length.out = nsegs + 1)
  y <- seq(y1, y2, length.out = nsegs + 1)
  return(cbind(x, y))
}
