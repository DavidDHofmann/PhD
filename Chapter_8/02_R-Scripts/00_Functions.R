################################################################################
#### Function to calculate Distances on a Raster Efficiently
################################################################################
#' Calculate Raster Distance
#'
#' Function to calculate the distance of a raster cell to the nearest cell
#' containing a certain value
#' @export
#' @param x \code{RasterLayer} on which distances should be calculated
#' @param value value to which the distance should be calculated
#' @return \code{RasterLayer}
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
#### Function to Extend Raster and Fill with Observed Values
################################################################################
#' Extend Raster and Fill with Observed Values
#'
#' Function to extend a \code{SpatRaster} an fill the extended area with values
#' sampled from the original \code{SpatRaster}
#' @export
#' @param x \code{SpatRaster}
#' @param y extent to which the raster should be extended
#' @return \code{SpatRaster}
extendRaster <- function(x, y) {

  # Create mask of NA values
  na_mask <- is.na(x)

  # Extract values from raster
  vals <- values(x, data.frame = F)

  # Remove na's
  vals <- na.omit(vals)

  # Extend raster to new extent
  r <- extend(x, y)

  # Also extend the mask
  na_mask <- extend(na_mask, y)

  # Identify NAs in the new extent
  indices <- which(is.na(values(r)))

  # Replace them with values sampled from old raster
  values(r)[indices] <- vals[runif(length(indices), min = 1, max = length(vals))]

  # Make sure original NA's are put back on the map
  r <- mask(r, na_mask, maskvalue = 1, updatevalue = NA)

  # Return the new raster
  return(r)
}


################################################################################
#### Function to Prepare Movement Model for Dispersal Simulation
################################################################################
#' Prepare Movement Model for Dispersal Simulation
#'
#' Function to prepare movement model for dispersal simulation
#' @export
#' @param model Movement model
#' @return \code{list}
prepareModel <- function(model = NULL){

  # Extract parameter estimates from model
  coeffs <- fixef(model)$cond

  # Remove intercept
  coeffs <- coeffs[names(coeffs) != "(Intercept)"]

  # Obtain original model formula and coerce it to a terms object
  form <- model
  form <- formula(form)
  form <- terms(form)

  # Identify and remove unnecessary terms
  newform <- c("Intercept", "step_id_", "0 +")
  newform <- paste(newform, collapse = "|")
  newform <- grep(newform, attr(form, "term.labels"))
  newform <- drop.terms(form, newform, keep.response = F)
  newform <- formula(newform)

  # Return the coefficients and the formula
  return(list(
      coefficients  = coeffs
    , formula       = newform
  ))
}

################################################################################
#### Function to Prepare Raster Layers for Dispersal Simulation
################################################################################
#' Prepare Covariates for Simulation
#'
#' Function to prepare covariate layers for dispersal simulation
#' @export
#' @param layers \code{RasterStack} of covariates
#' @return \code{list}
prepareCovars <- function(layers){

  # Extract layer names
  names <- names(layers)

  # If it is a terra raster, convert it
  if (inherits(layers, "SpatRaster")) {
    layers <- as(layers, "Raster")
  }

  # Extract layer extent
  extent <- extent(layers)
  extent <- list(
      xmin = extent@xmin
    , xmax = extent@xmax
    , ymin = extent@ymin
    , ymax = extent@ymax
  )

  # Convert layers to velox raster
  layers <- velox(layers)

  # Return a named list
  return(list(
      covars      = layers
    , covar_names = names
    , extent      = extent
  ))
}

################################################################################
#### Function to Check if Points are within or outside an extent
################################################################################
#' Check if Coordinates are within Extent
#'
#' Function to Check if Points are within or outside an extent
#' @export
#' @param \code{data.frame} or \code{matrix} containing columns representing x
#' and y coordinates for which you want to know if they lie within an extent or
#' not.
#' @param \code{list} containing xmin, xmax, ymin, and ymax
#' @return Logical (\code{TRUE} = within extent, \code{FALSE} = outside extent)
pointsInside <- function(xy = NULL, extent = NULL){

  # Check if x value is within boundaries
  xy[, 1] > extent$xmin & xy[, 1] < extent$xmax &

    # Check if y value is within boundaries
    xy[, 2] > extent$ymin & xy[, 2] < extent$ymax
}

################################################################################
#### Function to Extract Covariates (Fast)
################################################################################
#' Extract covariates quickly
#'
#' Function to extract raster values and calculate their average coverage along
#' a spatial feature. This function requires the velox package which allows a
#' much quicker value extraction.
#' @export
#' @param raster \code{velox raster} containing covariates
#' @param feature Spatial feature (\code{sp} package) below which covariates
#' should be extracted
#' @return \code{data.frame}
extrCov <- function(raster = NULL, feature = NULL) {

  # Copy velox raster
  raster_copy <- raster$copy()

  # Crop the raster
  raster_copy$crop(feature)

  # Extract values
  extracted <- raster_copy$extract(feature, small = TRUE)

  # Calculate average values
  extracted <- vapply(extracted, colMeans, numeric(raster_copy$nbands))

  # Coerce to dataframe
  extracted <- as.data.frame(t(extracted))

  # Return the final dataframe
  return(extracted)
}

################################################################################
#### Function to Scale a Covariate Using a Scaling Table
################################################################################
#' Scale Covariates
#'
#' Function to scale covariates using a scaling table
#' @export
#' @param covars Either a \code{data.frame} or a \code{RasterStack} containing
#' continuous covariates.
#' @param scaling \code{data.frame} containing a "Center" column and a "Scale"
#' column. The \code{data.frame} must contain rownames to indicating the
#' covariates to which the scaling parameters belong to.
#' @return \code{data.frame} or \code{RasterStack}
scaleCovars <- function (covars, scaling) {
  scaled <- sapply(1:ncol(covars), function(x) {
    scale(covars[, x]
      , center  = scaling$center[names(covars)[x]]
      , scale   = scaling$scale[names(covars)[x]]
    )
  })
  scaled <- as.data.frame(scaled)
  names(scaled) <- names(covars)
  return(scaled)
}

################################################################################
#### Function to Reproject Coordinates
################################################################################
#' Reproject coordinates to another projection
#'
#' Function to transform coordinates into coordinates of another projection
#' @export
#' @param xy \code{data.frame} or matrix of coordinates that need to be
#' transformed.
#' @param from CRS of the coordinates
#' @param to CRS to which the coordinates should be reprojected
#' @return numeric value or vector
#' @examples
#' degreesToMeters(0.01)
reprojCoords <- function(xy, from = NULL, to = NULL) {
  xy <- as.matrix(xy)
  xy <- suppressWarnings(vect(xy, crs = from))
  xy <- terra::project(xy, to)
  return(crds(xy))
}

################################################################################
#### Function to Coerce Single Simulated Trajectories to Proper Lines
################################################################################
#' Coerce multiple Simulated Trajectories to SpatialLines
#'
#' Function to coerce simulated coordinates to lines after a desired number of
#' steps. Basically a wrapper around \code{sim2tracks()}
#' @export
#' @param simulation \code{data.frame} resulting from the function
#' \code{disperse()}
#' @param keep.data logical Should the input dataframe be conserved?
#' @return \code{SpatialLinesDataFrame}
sim2tracks <- function(simulation = NULL, crs = NULL, keep.data = F) {

  # Unless keep.data is desired, remove all unnecessary columns
  if (!keep.data){
    simulation <- simulation[, c("x", "y")]
  }

  # Calculate number of steps
  steps <- nrow(simulation) - 1

  # If data does not need to be kept
  if (!keep.data){
    coordinates(simulation) <- c("x", "y")
    lines <- spLines(simulation)
    if (!is.null(crs)){
      crs(lines) <- crs
    }
    return(lines)

  # If data needs to be kept
  } else {
    pts <- simulation
    coordinates(pts) <- c("x", "y")
    line <- spLines(pts)
    lines <- createSegments(line)
    lines <- as(lines, "SpatialLinesDataFrame")
    lines@data <- simulation[1:steps, ]
    if (!is.null(crs)){
      crs(lines) <- crs
    }
    return(lines)
  }
}

################################################################################
#### Function to Coerce Multiple Simulated Trajectories to Proper Lines
################################################################################
#' Coerce multiple Simulated Trajectories to SpatialLines
#'
#' Function to coerce simulated coordinates to lines after a desired number of
#' steps. Basically a wrapper around \code{sim2tracks()}
#' @export
#' @param simulations \code{data.frame} resulting from the function
#' \code{disperse}
#' @param id character column name of the column containing the track id
#' @param keep.data logical Should the input dataframe be conserved?
#' @param mc.cores integer Number of cores used for the coercion. Set to
#' 1 by default.
#' @param messages logical Should a progress bar and messages be printed?
#' @return \code{SpatialLinesDataFrame}
sims2tracks <- function(
    simulations = NULL
  , id          = "TrackID"
  , crs         = NULL
  , keep.data   = F
  , mc.cores    = 1
  , messages    = T
  ) {

  # Keep only relevant data
  if (!keep.data){
    simulations <- simulations[, c("x", "y", "TrackID")]
  }

  # Nest data by id
  simulations <- nest(simulations, data = -all_of(id))

  # Coerce each trajectory to a spatial lines object
  if (messages){
    lines <- pbmclapply(1:nrow(simulations)
    , mc.cores           = mc.cores
    , ignore.interactive = T
    , FUN                = function(x){
      l <- sim2tracks(
          simulation = simulations$data[[x]]
        , crs        = crs
        , keep.data  = keep.data
      )
      l$TrackID <- simulations$TrackID[x]
      return(l)
    })
  } else {
    lines <- mclapply(1:nrow(simulations)
    , mc.cores           = mc.cores
    , FUN                = function(x){
      l <- sim2tracks(
          simulation = simulations$data[[x]]
        , keep.data  = keep.data
      )
      l$TrackID <- simulations$TrackID[x]
      return(l)
    })
  }

  # Bind them
  lines <- do.call(rbind, lines)

  # Return the lines
  return(lines)

}

################################################################################
#### Function to Rasterize Using SpatStat
################################################################################
#' Rasterize Shapes Using spatstat
#'
#' Function to rasterize using spatstat
#' @export
#' @param l \code{SpatialPoints}, \code{SpatialLines}, or \code{SpatialPolygons}
#' to be rasterized
#' @param r \code{RasterLayer} onto which the objects should be rasterized
#' @param mc.cores numeric How many cores should be used?
#' @return \code{RasterLayer}
rasterizeSpatstat <- function(l, r, mc.cores = 1){

  # In case we run the rasterization on a single core, run a loop
  if (mc.cores == 1){

    # Create im layer
    values(r) <- 0
    im <- as.im.RasterLayer(r)
    summed <- im

    # Prepare progress bar
    pb <- txtProgressBar(
        min     = 0
      , max     = length(l)
      , initial = 0
      , style   = 3
      , width   = 55
    )

    # Go through all lines and rasterize them
    for (y in 1:length(l)){
      line    <- as.psp(l[y, ], window = im)
      line    <- as.mask.psp(line)
      line_r  <- as.im.owin(line, na.replace = 0)
      summed  <- Reduce("+", list(summed, line_r))
      setTxtProgressBar(pb, y)
    }

    # Return heatmap as a raster
    return(raster(summed))

  # If multicore is desired
  } else {

    # Split lines into a package for each core
    lines <- splitShape(l, n = mc.cores)

    # Run in parallel
    heatmap <- mclapply(lines, mc.cores = mc.cores, function(x){

      # Create im layer
      im <- as.im.RasterLayer(r)
      summed <- im

      # Rasterize each line
      for (y in 1:length(x)){
        line    <- as.psp(x[y, ], window = im)
        line    <- as.mask.psp(line)
        line_r  <- as.im.owin(line, na.replace = 0)
        summed  <- Reduce("+", list(summed, line_r))
      }

      # Return the sum
      return(summed)

    })

    # Combine heatmaps
    heatmap <- Reduce("+", heatmap)

    # Return as raster
    return(raster(heatmap))

  }
}

################################################################################
#### Function to Rasterize Simulations
################################################################################
#' Rasterize Simulations
#'
#' Function to rasterize using spatstat
#' @export
#' @param simulations \code{data.frame} containing the simulated data
#' @param raster \code{RasterLayer} onto which the objects should be rasterized
#' @param steps numeric, how many steps should be rasterized
#' @param area numeric, ID of the source areas that should be considered
#' @param flood character, one of "Min", "Mean", "Max"
#' @param messages boolean, should messages be printed during the rasterization
#' @param mc.cores numeric How many cores should be used?
#' @return \code{RasterLayer}
rasterizeSims <- function(
      simulations = NULL      # Simulated trajectories
    , raster      = NULL      # Raster onto which we rasterize
    , steps       = 500       # How many steps should be considered
    , area        = NULL      # Simulations from which areas?
    , flood       = "Min"     # Which flood level?
    , messages    = T         # Print update messages?
    , mc.cores    = detectCores() - 1
  ) {

  # Subset to corresponding data
  sub <- simulations[simulations$StepNumber <= steps, ]
  sub <- simulations[simulations$FloodLevel == flood, ]
  if (!is.null(area)) {
    sub <- simulations[simulations$Area %in% area, ]
  }

  # Make sure raster values are all 0
  values(raster) <- 0

  # Create spatial lines
  sub_traj <- sims2tracks(
      simulations = sub
    , id          = "TrackID"
    , messages    = messages
    , mc.cores    = mc.cores
  )

  # Convert into spatial lines
  sub_traj <- as(sub_traj, "SpatialLines")
  crs(sub_traj) <- CRS("+init=epsg:32734")

  # Rasterize lines onto the cropped raster
  if (messages) {
    cat("Rasterizing spatial lines...\n")
  }
  heatmap <- rasterizeSpatstat(
      l        = sub_traj
    , r        = raster
    , mc.cores = 1
  )

  # Store heatmap to temporary file
  heatmap <- writeRaster(heatmap, tempfile())
  crs(heatmap) <- CRS("+init=epsg:32734")

  # Return the resulting heatmap
  return(heatmap)
}
